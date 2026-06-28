/*
 * test_scaling_cuda.c — CUDA-accelerated headless test for AND3tot scaling
 *
 * Build (Windows MSVC + CUDA):
 *   nvcc -c ca_cuda_sinc.cu -o ca_cuda_sinc.obj -DNO_SDL -DUSE_CUDA [-DL=161]
 *   cl /Ox /W3 /DNO_SDL /DUSE_CUDA [-DL=161] /Fe:test_cuda.exe ^
 *      test_scaling_cuda.c integrated.c ca_cuda_sinc.obj cudart.lib
 *
 * Build (Linux):
 *   nvcc -c ca_cuda_sinc.cu -o ca_cuda_sinc.o -DNO_SDL -DUSE_CUDA [-DL=161]
 *   gcc -O2 -DNO_SDL -DUSE_CUDA [-DL=161] -o test_cuda \
 *       test_scaling_cuda.c integrated.c ca_cuda_sinc.o -lcudart -lm
 *
 * The sinc_step() and pulse_step() in integrated.c are compiled out
 * when USE_CUDA is defined (#ifndef USE_CUDA). Instead, CUDA kernels
 * are called from this driver. init() and generate_spiral() still
 * run on CPU.
 */

#define NO_SDL
#include "integrated.h"
#include "ca_cuda_sinc.h"

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* sinc_converged and check_convergence() declared in integrated.h */

int main(void)
{
    /* ---- allocate host grids ---- */
    size_t grid_bytes = sizeof(Cell) * (size_t)L * L * L;
    grid      = (Cell (*)[L][L])malloc(grid_bytes);
    grid_next = (Cell (*)[L][L])malloc(grid_bytes);
    if (!grid || !grid_next) {
        printf("OOM for L=%d (need %lu MB host RAM)\n", L,
               (unsigned long)(2UL * grid_bytes / (1024*1024)));
        return 1;
    }

    /* ---- CPU init: seed wave, BFS, spiral geometry ---- */
    init();
    printf("L=%d R_MAX=%d spiral_n=%d PULSE_TOLERANCE=%d (CUDA)\n",
           L, R_MAX, spiral_n, PULSE_TOLERANCE);

    /* ---- allocate + upload to GPU ---- */
    cuda_alloc_grids();
    cuda_upload_grid(grid, grid_next);

    /* ---- cycle detection constants ---- */
    const unsigned int max_r2 =
        (unsigned int)((unsigned int)R_MAX * R_MAX * 92 / 100);
    const unsigned int span   = max_r2;
    const unsigned int period = span + span;
    unsigned int target_cycles = 5;
    unsigned int max_ticks =
        (period / PULSE_STEP) * (target_cycles + 2) + 2000;
    if (max_ticks > 2000000) max_ticks = 2000000;

    int cycles_seen = 0;
    unsigned int prev_cycle = 0;
    int converged = 0;

    /* running totals (accumulated on CPU from per-tick GPU counts) */
    int and2_total = 0;
    int and3_total = 0;
    int and2_last  = 0;
    int and3_last  = 0;

    /* convergence check interval (download grid every N ticks) */
    int conv_check_interval = 10;

    printf("Running up to %u ticks ...\n", max_ticks);
    fflush(stdout);

    for (unsigned int t = 0; t < max_ticks; t++) {
        /* ---- one step on GPU ---- */
        cuda_reset_and_counters();
        cuda_pulse_step(tick);
        cuda_sinc_step(tick);
        tick++;

        /* ---- download per-tick AND counters (cheap: 2 ints) ---- */
        int tick_and2 = 0, tick_and3 = 0;
        cuda_download_and_counters(&tick_and2, &tick_and3);
        and2_total += tick_and2;
        and3_total += tick_and3;

        /* ---- cycle detection (on CPU, using pulse arithmetic) ---- */
        {
            unsigned int cycle_now =
                ((unsigned int)tick * PULSE_STEP) / period;
            if (cycle_now != prev_cycle && prev_cycle > 0) {
                and2_last = and2_total;
                and3_last = and3_total;
                and2_total = 0;
                and3_total = 0;

                if (converged) {
                    cycles_seen++;
                    printf("  Cycle %d ended: AND2tot=%d AND3tot=%d\n",
                           cycles_seen, and2_last, and3_last);
                    fflush(stdout);
                    if (cycles_seen >= (int)target_cycles) {
                        prev_cycle = cycle_now;
                        break;
                    }
                }
            }
            prev_cycle = cycle_now;
        }

        /* ---- convergence detection (periodic, on CPU) ---- */
        if (!converged && (t % conv_check_interval == 0)) {
            cuda_download_grid(grid);
            check_convergence();
            if (sinc_converged) {
                converged = 1;
                printf("  Sinc converged at tick %d — AND triple active.\n",
                       tick);
                fflush(stdout);
                /* upload updated grid (sinc_p/sinc_q set by convergence) */
                cuda_upload_grid_full(grid);
            }
        }
    }

    printf("RESULT: L=%d spiral_n=%d AND3tot=%d\n\n",
           L, spiral_n, and3_last);

    cuda_free();
    free(grid);
    free(grid_next);
    return 0;
}
