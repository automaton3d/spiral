/* Headless geometry check: build with -DNO_SDL */
#include "spiral.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#ifndef M_PI
#define M_PI 3.14159265358979323846
#endif

static void run(int ax, int ay, int az) {
    spiral_set_axis(ax, ay, az);
    memset(grid, 0, sizeof(Cell) * L * L * L);
    memset(grid_next, 0, sizeof(Cell) * L * L * L);
    for (int x = 0; x < L; x++) for (int y = 0; y < L; y++) for (int z = 0; z < L; z++)
        grid[x][y][z].spin = 0;
    tick = 0;
    spiral_init();

    int guard = 200000;
    while (!spiral_done && guard--) spiral_step();

    /* measure: distance from the cylinder axis line, r range, angle swept */
    double A[3] = { axis_x, axis_y, axis_z };
    double an = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
    double u0[3], last_ang = 0, total_ang = 0;
    double rmin = 1e9, rmax = -1e9, rr_max = 0;
    double e1[3], e2[3];
    int first = 1;

    /* basis perpendicular to A */
    {
        double w[3] = {0,0,0};
        int wi = (fabs(A[0]) <= fabs(A[1]) && fabs(A[0]) <= fabs(A[2])) ? 0 :
                 (fabs(A[1]) <= fabs(A[2]) ? 1 : 2);
        w[wi] = 1;
        e1[0] = A[1]*w[2]-A[2]*w[1]; e1[1] = A[2]*w[0]-A[0]*w[2]; e1[2] = A[0]*w[1]-A[1]*w[0];
        double n1 = sqrt(e1[0]*e1[0]+e1[1]*e1[1]+e1[2]*e1[2]);
        for (int i=0;i<3;i++) e1[i] /= n1;
        e2[0] = A[1]*e1[2]-A[2]*e1[1]; e2[1] = A[2]*e1[0]-A[0]*e1[2]; e2[2] = A[0]*e1[1]-A[1]*e1[0];
        double n2 = sqrt(e2[0]*e2[0]+e2[1]*e2[1]+e2[2]*e2[2]);
        for (int i=0;i<3;i++) e2[i] /= n2;
    }

    /* cylinder axis offset P estimated from the first point set:
       use the mean of the projected points as the circle centre */
    double mx = 0, my = 0;
    for (int i = 0; i < spiral_n; i++) {
        double v[3] = { spiral_pts[i].x - MID, spiral_pts[i].y - MID, spiral_pts[i].z - MID };
        mx += v[0]*e1[0]+v[1]*e1[1]+v[2]*e1[2];
        my += v[0]*e2[0]+v[1]*e2[1]+v[2]*e2[2];
    }
    mx /= spiral_n; my /= spiral_n;

    double hmin = 1e9, hmax = -1e9;
    for (int i = 0; i < spiral_n; i++) {
        double v[3] = { spiral_pts[i].x - MID, spiral_pts[i].y - MID, spiral_pts[i].z - MID };
        double p1 = v[0]*e1[0]+v[1]*e1[1]+v[2]*e1[2] - mx;
        double p2 = v[0]*e2[0]+v[1]*e2[1]+v[2]*e2[2] - my;
        double h  = (v[0]*A[0]+v[1]*A[1]+v[2]*A[2]) / an;
        double rad = sqrt(p1*p1 + p2*p2);
        if (rad < rmin) rmin = rad;
        if (rad > rmax) rmax = rad;
        double rr = sqrt(v[0]*v[0]+v[1]*v[1]+v[2]*v[2]);
        if (rr > rr_max) rr_max = rr;
        if (h < hmin) hmin = h;
        if (h > hmax) hmax = h;
        double ang = atan2(p2, p1);
        if (first) { last_ang = ang; first = 0; u0[0]=p1; u0[1]=p2; u0[2]=0; }
        else {
            double d = ang - last_ang;
            while (d >  M_PI) d -= 2*M_PI;
            while (d < -M_PI) d += 2*M_PI;
            total_ang += d;
            last_ang = ang;
        }
    }
    (void)u0;
    printf("axis(%4d,%4d,%4d)->(%4d,%4d,%4d) |A|=%.1f  pts=%4d  cyl r: %.2f..%.2f (R_CYL=%d)"
           "  r_max=%.1f (RADIUS=%d)  climb=%.1f  turn=%.2f rad\n",
           ax, ay, az, axis_x, axis_y, axis_z, an, spiral_n,
           rmin, rmax, R_CYL, rr_max, RADIUS, hmax - hmin, total_ang);
}

int main(void) {
    grid      = malloc(sizeof(Cell) * L * L * L);
    grid_next = malloc(sizeof(Cell) * L * L * L);
    if (!grid || !grid_next) { printf("oom\n"); return 1; }

    run(0, 0, 1);
    run(0, 0, -1);
    run(1, 0, 0);
    run(0, 1, 0);
    run(1, 1, 1);
    run(3, -2, 5);
    run(-7, 4, 1);
    run(1, 2, 0);

    free(grid); free(grid_next);
    return 0;
}
