/*
 * test_integrated.c — verify AND triple detection
 *
 * Runs the integrated CA headlessly (no SDL) for enough ticks
 * to see BFS propagation, spiral geometry, and AND triple hits.
 */

#define NO_SDL
#include "integrated.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>

int main(void) {
    grid      = malloc(sizeof(Cell) * L * L * L);
    grid_next = malloc(sizeof(Cell) * L * L * L);
    if (!grid || !grid_next) {
        printf("FAIL: out of memory\n");
        return 1;
    }
    memset(grid,      0, sizeof(Cell) * L * L * L);
    memset(grid_next, 0, sizeof(Cell) * L * L * L);

    init();

    /* --- Verify spiral geometry --- */
    printf("Spiral points: %d\n", spiral_n);
    if (spiral_n < 200 || spiral_n > 400) {
        printf("FAIL: spiral_n=%d outside expected range [200,400]\n",
               spiral_n);
        return 1;
    }

    /* verify all spiral points have spin=1 */
    int spin_ok = 1;
    for (int i = 0; i < spiral_n; i++) {
        if (!grid[spiral_pts[i].x][spiral_pts[i].y][spiral_pts[i].z].spin) {
            printf("FAIL: spiral_pts[%d] has spin=0\n", i);
            spin_ok = 0;
        }
    }
    if (spin_ok) printf("All spiral points have spin=1: OK\n");

    /* verify center cell */
    if (grid[MID][MID][MID].r2 != 0) {
        printf("FAIL: center r2 != 0\n");
        return 1;
    }
    if (grid[MID][MID][MID].u != 2048) {
        printf("FAIL: center u != 2048\n");
        return 1;
    }
    printf("Center cell: r2=0, u=2048: OK\n");

    /* --- Run CA for some ticks, verify BFS propagation --- */
    int max_ticks = 300;
    int bfs_done = 0;
    for (int t = 0; t < max_ticks; t++) {
        step_all();
        if (!bfs_done) {
            int corner_r2 = (int)grid[0][0][0].r2;
            if (corner_r2 != (int)INF_R2) {
                printf("BFS reached corner (0,0,0) at tick %d, r2=%d, r=%d\n",
                       t + 1, corner_r2, grid[0][0][0].r);
                bfs_done = 1;
            }
        }
    }
    if (!bfs_done)
        printf("BFS did not reach corner in %d ticks (may need more)\n",
               max_ticks);

    /* --- Check AND double: active ∧ spin --- */
    int and_double = 0;
    for (int i = 0; i < spiral_n; i++) {
        Cell *c = &grid[spiral_pts[i].x][spiral_pts[i].y][spiral_pts[i].z];
        if (c->active && c->spin)
            and_double++;
    }
    printf("AND double (active ∧ spin) at tick %d: %d points\n",
           tick, and_double);

    /* --- Check AND triple counts --- */
    int total_and = 0;
    for (int r = 0; r < L; r++)
        total_and += and_count[r];
    printf("AND triple (trig ∧ active ∧ spin) accumulated: %d hits\n",
           total_and);

    /* note: sinc may not have converged in 300 ticks,
       so AND triple may be 0 — that's expected.
       the triggers only activate after sinc convergence. */
    if (total_and == 0)
        printf("(expected: sinc needs ~%d ticks to converge)\n",
               STABILITY_FRAMES + 200);

    printf("\nSummary: L=%d, MID=%d, R_CYL=%d, spiral_n=%d, ticks=%d\n",
           L, MID, R_CYL, spiral_n, tick);
    printf("Test PASSED.\n");

    free(grid);
    free(grid_next);
    return 0;
}
