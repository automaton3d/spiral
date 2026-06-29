/*
 * test_pulse_r.c — verify that sum-of-odds pulse_r matches isqrt(pulse_r2)
 *                  at every tick for several full cycles.
 */

#define NO_SDL
#include "integrated.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

int main(void) {
    grid      = malloc(sizeof(Cell) * L * L * L);
    grid_next = malloc(sizeof(Cell) * L * L * L);
    if (!grid || !grid_next) {
        printf("OOM for L=%d\n", L);
        return 1;
    }

    init();
    printf("L=%d R_MAX=%d PULSE_STEP=%d\n", L, R_MAX, PULSE_STEP);

    const unsigned int max_r2 =
        (unsigned int)((unsigned int)R_MAX * R_MAX * 92 / 100);
    /* run 3 full cycles */
    unsigned int ticks_per_cycle = max_r2 + max_r2;  /* PULSE_STEP=1 */
    unsigned int total_ticks = ticks_per_cycle * 3 + 100;

    int errors = 0;
    int checks = 0;

    for (unsigned int t = 0; t < total_ticks; t++) {
        step_all();
        checks++;

        /* verify: pulse_r_state == isqrt(pulse_r2_state) */
        int expected_r = isqrt((int)pulse_r2_state);
        if (pulse_r_state != expected_r) {
            printf("  MISMATCH at tick %d: pulse_r2=%u  pulse_r=%d  isqrt=%d\n",
                   tick, pulse_r2_state, pulse_r_state, expected_r);
            errors++;
            if (errors > 20) {
                printf("  Too many errors, stopping.\n");
                break;
            }
        }

        /* verify: pulse_r_sq == pulse_r_state * pulse_r_state */
        unsigned int expected_sq = (unsigned int)(pulse_r_state * pulse_r_state);
        if (pulse_r_sq != expected_sq) {
            printf("  SQ MISMATCH at tick %d: pulse_r=%d  pulse_r_sq=%u  expected=%u\n",
                   tick, pulse_r_state, pulse_r_sq, expected_sq);
            errors++;
        }

        /* verify: pulse_gap == 2*pulse_r_state + 1 */
        int expected_gap = 2 * pulse_r_state + 1;
        if (pulse_gap != expected_gap) {
            printf("  GAP MISMATCH at tick %d: pulse_r=%d  pulse_gap=%d  expected=%d\n",
                   tick, pulse_r_state, pulse_gap, expected_gap);
            errors++;
        }
    }

    printf("\nChecked %d ticks over 3 cycles.\n", checks);
    if (errors == 0)
        printf("PASS: pulse_r == isqrt(pulse_r2) at every tick.\n");
    else
        printf("FAIL: %d mismatches found.\n", errors);

    /* Also verify per-cell r values */
    int cell_errors = 0;
    int cells_checked = 0;
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        if (grid[x][y][z].r2 == INF_R2) continue;
        cells_checked++;
        int expected = isqrt((int)grid[x][y][z].r2);
        if (grid[x][y][z].r != expected) {
            if (cell_errors < 5)
                printf("  CELL r MISMATCH at (%d,%d,%d): r2=%u  r=%d  isqrt=%d\n",
                       x, y, z, grid[x][y][z].r2, grid[x][y][z].r, expected);
            cell_errors++;
        }
    }
    printf("Per-cell r: checked %d cells, %d mismatches.\n",
           cells_checked, cell_errors);

    free(grid);
    free(grid_next);
    return errors + cell_errors;
}
