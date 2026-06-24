/*
 * test_scaling.c — headless test to measure AND3 total per cycle vs L
 * Compile with: gcc -O2 -DNO_SDL -DL=N -o test_LN test_scaling.c integrated.c -lm
 */

#define NO_SDL
#include "integrated.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

/* externs now in integrated.h */

int main(void) {
    grid      = malloc(sizeof(Cell) * L * L * L);
    grid_next = malloc(sizeof(Cell) * L * L * L);
    if (!grid || !grid_next) {
        printf("OOM for L=%d (need %lu MB)\n", L,
               (unsigned long)(2UL * sizeof(Cell) * L * L * L / (1024*1024)));
        return 1;
    }

    init();
    printf("L=%d R_MAX=%d spiral_n=%d PULSE_TOLERANCE=%d\n",
           L, R_MAX, spiral_n, PULSE_TOLERANCE);

    /* Run until convergence + 5 full cycles */
    const unsigned int max_r2 =
        (unsigned int)((unsigned int)R_MAX * R_MAX * 92 / 100);
    const unsigned int span = max_r2;
    const unsigned int period = span + span;
    unsigned int target_cycles = 5;
    unsigned int max_ticks = (period / PULSE_STEP) * (target_cycles + 2) + 2000;

    /* cap at reasonable limit (PULSE_STEP=1 needs longer runs) */
    if (max_ticks > 2000000) max_ticks = 2000000;

    int cycles_seen = 0;
    unsigned int prev_cycle = 0;
    int converged = 0;

    for (unsigned int t = 0; t < max_ticks; t++) {
        step_all();

        /* detect convergence by first AND triple appearance */
        if (!converged) {
            if (and_triple_count > 0) {
                converged = 1;
                printf("  First AND3 at tick %d\n", tick);
            }
        }

        /* track cycles after convergence */
        if (converged) {
            unsigned int cycle_now = ((unsigned int)tick * PULSE_STEP) / period;
            if (cycle_now != prev_cycle) {
                cycles_seen++;
                printf("  Cycle %d ended: AND2tot=%d AND3tot=%d\n",
                       cycles_seen, and_double_last, and_triple_last);
                prev_cycle = cycle_now;
                if (cycles_seen >= (int)target_cycles) break;
            }
        }
    }

    printf("RESULT: L=%d spiral_n=%d AND3tot=%d\n\n",
           L, spiral_n, and_triple_last);

    free(grid);
    free(grid_next);
    return 0;
}
