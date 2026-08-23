/*
 * spiral.c — emergent spiral CA on pulsating wavefront
 *            around an ARBITRARY rotation axis (|AXIS| = L/2)
 *
 * Walker:
 *   - starts at the lattice centre (MID,MID,MID), which lies exactly on
 *     a cylinder of radius R_CYL whose axis is parallel to AXIS and
 *     offset by the perpendicular vector P (|P| = R_CYL);
 *   - each tick a Bresenham accumulator selects an orbital step (one of
 *     the 6 lattice moves, chosen to keep |u x AXIS|^2 closest to
 *     R_CYL^2 |AXIS|^2 while advancing along the tangent t = AXIS x u)
 *     or a climb step (6-neighbour DDA along AXIS);
 *   - it stops when the wave radius r reaches L/2 (theta = pi).
 *
 * Runtime operations: +, -, <<, comparison.  No multiplication,
 * no division, no floats, no tables indexed by angle.
 * (Init-time constants, pulse_from_time and rendering are host-only.)
 */

#include "spiral.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>

/* --- Grid allocation (heap, too large for stack) --- */
Cell (*grid)[L][L]      = NULL;
Cell (*grid_next)[L][L] = NULL;

int tick = 0;

/* --- Rotation axis (|AXIS| = R_MAX = L/2) --- */
int axis_x = AXIS_X, axis_y = AXIS_Y, axis_z = AXIS_Z;

/* --- Spiral walker state --- */
SpiralPt spiral_pts[MAX_SPIRAL_PTS];
int spiral_n    = 0;
int spiral_done = 0;

/* Lattice moves: 0:+x 1:-x 2:+y 3:-y 4:+z 5:-z */
static const int mvx[6] = { 1, -1, 0,  0, 0,  0 };
static const int mvy[6] = { 0,  0, 1, -1, 0,  0 };
static const int mvz[6] = { 0,  0, 0,  0, 1, -1 };

/* ---- init-time constants (axis dependent) ---- */
static long long Dx[6], Dy[6], Dz[6];  /* D[m] = e_m x AXIS            */
static long long K[6][6];              /* K[k][m] = 2 * D_k . D_m      */
static long long TGT;                  /* R_CYL^2 * |AXIS|^2           */
static long long TOL;                  /* radial band ~ +/-1 cell      */
static int absa[3];                    /* |ax|,|ay|,|az|               */
static int climb_total;                /* |ax|+|ay|+|az| (DDA period)  */
static int orbit_total;                /* CYL_CIRC                     */
static int period_total;               /* climb + orbit (Bresenham)    */

/* ---- runtime walker state (integers only) ---- */
static int tip_x, tip_y, tip_z;        /* lattice position             */
static int vx, vy, vz;                 /* offset from lattice centre   */
static long long q2;                   /* |v|^2 (sum of odds)          */
static long long cxv, cyv, czv;        /* c = u x AXIS, u = v - P      */
static long long E;                    /* |c|^2                        */
static long long G[6];                 /* E increment for move m       */
static int tip_acc;                    /* orbital/climb Bresenham      */
static int dda[3];                     /* climb DDA accumulators       */

/* ==========================================================
 * spiral_set_axis — install an arbitrary rotation axis.
 *
 * The vector is rescaled to |AXIS| = L/2 and every axis-dependent
 * constant is derived here (init-time arithmetic may multiply).
 * ========================================================== */
void spiral_set_axis(int ax, int ay, int az) {
    if (ax == 0 && ay == 0 && az == 0) { ax = 0; ay = 0; az = R_MAX; }

    /* normalise to length R_MAX = L/2 */
    {
        double n = sqrt((double)ax * ax + (double)ay * ay + (double)az * az);
        double s = (double)R_MAX / n;
        axis_x = (int)floor((double)ax * s + 0.5);
        axis_y = (int)floor((double)ay * s + 0.5);
        axis_z = (int)floor((double)az * s + 0.5);
        if (axis_x == 0 && axis_y == 0 && axis_z == 0) axis_z = R_MAX;
    }

    absa[0] = axis_x < 0 ? -axis_x : axis_x;
    absa[1] = axis_y < 0 ? -axis_y : axis_y;
    absa[2] = axis_z < 0 ? -axis_z : axis_z;
    climb_total = absa[0] + absa[1] + absa[2];
    /* Orbital step budget for one revolution.
     * A lattice circle of radius R_CYL in the plane perpendicular to AXIS
     * needs  R_CYL * integral(|t|_1 dtheta)  von-Neumann moves, where t is
     * the unit tangent.  For AXIS = z this integral is 8 (the classic 4/pi
     * correction, CYL_CIRC = 8R).  For a tilted axis it is larger, so it is
     * evaluated once here at init time — the CA runtime never sees it. */
    {
        double A[3] = { axis_x, axis_y, axis_z };
        double an = sqrt(A[0]*A[0] + A[1]*A[1] + A[2]*A[2]);
        double w[3] = { 0, 0, 0 };
        int wi = 0;
        double aa[3] = { fabs(A[0]), fabs(A[1]), fabs(A[2]) };
        if (aa[1] < aa[wi]) wi = 1;
        if (aa[2] < aa[wi]) wi = 2;
        w[wi] = 1.0;

        double e1[3], e2[3], n1, n2;
        e1[0] = A[1]*w[2] - A[2]*w[1];
        e1[1] = A[2]*w[0] - A[0]*w[2];
        e1[2] = A[0]*w[1] - A[1]*w[0];
        n1 = sqrt(e1[0]*e1[0] + e1[1]*e1[1] + e1[2]*e1[2]);
        for (int i = 0; i < 3; i++) e1[i] /= n1;
        e2[0] = A[1]*e1[2] - A[2]*e1[1];
        e2[1] = A[2]*e1[0] - A[0]*e1[2];
        e2[2] = A[0]*e1[1] - A[1]*e1[0];
        n2 = sqrt(e2[0]*e2[0] + e2[1]*e2[1] + e2[2]*e2[2]);
        for (int i = 0; i < 3; i++) e2[i] /= n2;
        (void)an;

        int N = 720;
        double acc = 0.0;
        for (int i = 0; i < N; i++) {
            double th = 6.283185307179586 * (double)i / (double)N;
            double s = sin(th), c = cos(th);
            /* tangent = -sin*e1 + cos*e2 */
            double tx = -s*e1[0] + c*e2[0];
            double ty = -s*e1[1] + c*e2[1];
            double tz = -s*e1[2] + c*e2[2];
            acc += (fabs(tx) + fabs(ty) + fabs(tz)) * 6.283185307179586 / (double)N;
        }
        orbit_total = (int)floor((double)R_CYL * acc + 0.5);
        if (orbit_total < CYL_CIRC) orbit_total = CYL_CIRC;
        period_total = orbit_total + climb_total;
    }


    /* D[m] = e_m x AXIS  (constant increment of c per lattice move) */
    Dx[0] =  0;               Dy[0] = -axis_z;  Dz[0] =  axis_y;   /* +x */
    Dx[1] =  0;               Dy[1] =  axis_z;  Dz[1] = -axis_y;   /* -x */
    Dx[2] =  axis_z;          Dy[2] =  0;       Dz[2] = -axis_x;   /* +y */
    Dx[3] = -axis_z;          Dy[3] =  0;       Dz[3] =  axis_x;   /* -y */
    Dx[4] = -axis_y;          Dy[4] =  axis_x;  Dz[4] =  0;        /* +z */
    Dx[5] =  axis_y;          Dy[5] = -axis_x;  Dz[5] =  0;        /* -z */

    /* K[k][m] = 2 * D_k . D_m  (how G[m] shifts after taking move k) */
    for (int k = 0; k < 6; k++)
        for (int m = 0; m < 6; m++)
            K[k][m] = 2 * (Dx[k] * Dx[m] + Dy[k] * Dy[m] + Dz[k] * Dz[m]);

    /* target: R_CYL^2 * |AXIS|^2 */
    {
        long long a2 = (long long)axis_x * axis_x
                     + (long long)axis_y * axis_y
                     + (long long)axis_z * axis_z;
        TGT = (long long)R_CYL_SQ * a2;
        TOL = (long long)(R_CYL + R_CYL + 1) * a2;   /* (R+1)^2 - R^2 */
    }
}

/* ==========================================================
 * Initialization
 * ========================================================== */
void init(void) {
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        Cell *c  = &grid[x][y][z];
        c->r     = 0;
        c->r2    = INF_R2;
        c->active = 0;
        c->spin  = 0;

        Cell *cn  = &grid_next[x][y][z];
        cn->r     = 0;
        cn->r2    = INF_R2;
        cn->active = 0;
        cn->spin  = 0;
    }

    /* seed center */
    grid[MID][MID][MID].r2 = 0;
    grid[MID][MID][MID].r  = 0;

    /* init walker */
    spiral_init();
}

/* ==========================================================
 * spiral_init — seed the helix walker at the lattice centre.
 *
 * P is a vector perpendicular to AXIS with |P| ~= R_CYL: the
 * cylinder axis is the line {P + t*AXIS}, so the centre (v = 0)
 * sits exactly on the cylinder surface, as in the +z version
 * (where P was (0, R_CYL, 0)).
 * ========================================================== */
void spiral_init(void) {
    int Px, Py, Pz;

    if (climb_total == 0) spiral_set_axis(axis_x, axis_y, axis_z);

    /* --- build P: perpendicular to AXIS, length R_CYL --- */
    {
        int wi = 0;                     /* least-aligned cardinal */
        if (absa[1] < absa[wi]) wi = 1;
        if (absa[2] < absa[wi]) wi = 2;

        /* P0 = AXIS x e_w */
        double p0x, p0y, p0z;
        if (wi == 0)      { p0x = 0;              p0y =  (double)axis_z; p0z = -(double)axis_y; }
        else if (wi == 1) { p0x = -(double)axis_z; p0y = 0;              p0z =  (double)axis_x; }
        else              { p0x =  (double)axis_y; p0y = -(double)axis_x; p0z = 0; }

        double n = sqrt(p0x * p0x + p0y * p0y + p0z * p0z);
        double s = (double)R_CYL / n;
        Px = (int)floor(p0x * s + 0.5);
        Py = (int)floor(p0y * s + 0.5);
        Pz = (int)floor(p0z * s + 0.5);
    }

    tip_x = MID; tip_y = MID; tip_z = MID;
    vx = 0; vy = 0; vz = 0;
    q2 = 0;

    /* u = v - P ; c = u x AXIS */
    {
        long long ux = -Px, uy = -Py, uz = -Pz;
        cxv = uy * axis_z - uz * axis_y;
        cyv = uz * axis_x - ux * axis_z;
        czv = ux * axis_y - uy * axis_x;
    }
    E = cxv * cxv + cyv * cyv + czv * czv;

    /* G[m] = 2 (c . D_m) + |D_m|^2  — the generalised "2dx+1" */
    for (int m = 0; m < 6; m++) {
        G[m] = 2 * (cxv * Dx[m] + cyv * Dy[m] + czv * Dz[m])
             + (Dx[m] * Dx[m] + Dy[m] * Dy[m] + Dz[m] * Dz[m]);
    }

    tip_acc = 0;
    dda[0] = dda[1] = dda[2] = 0;
    spiral_done = 0;

    grid[MID][MID][MID].spin = 1;
    spiral_pts[0].x = MID;
    spiral_pts[0].y = MID;
    spiral_pts[0].z = MID;
    spiral_n = 1;
}

/* ==========================================================
 * BFS wavefront propagation (sum-of-odds for r2)
 * ========================================================== */
static const int ddx[6] = {1, -1, 0,  0, 0,  0};
static const int ddy[6] = {0,  0, 1, -1, 0,  0};
static const int ddz[6] = {0,  0, 0,  0, 1, -1};

/* ==========================================================
 * Pure CA pulse update: each cell computes its next r2 from the
 * minimum among itself and its 6 neighbours (sum-of-odds).
 * ========================================================== */
static void pulse_update_wavefront(void) {
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        unsigned int best = grid[x][y][z].r2;

        for (int d = 0; d < 6; d++) {
            int nx = x + ddx[d];
            int ny = y + ddy[d];
            int nz = z + ddz[d];
            if (nx < 0 || nx >= L || ny < 0 || ny >= L ||
                nz < 0 || nz >= L)
                continue;

            unsigned int nr2 = grid[nx][ny][nz].r2;
            if (nr2 == INF_R2) continue;

            unsigned int diff;
            if (nx != x) {
                diff = (unsigned int)((nx > MID) ? (nx - MID) : (MID - nx));
                diff = (diff << 1) + 1;
            } else if (ny != y) {
                diff = (unsigned int)((ny > MID) ? (ny - MID) : (MID - ny));
                diff = (diff << 1) + 1;
            } else {
                diff = (unsigned int)((nz > MID) ? (nz - MID) : (MID - nz));
                diff = (diff << 1) + 1;
            }

            unsigned int cand = nr2 + diff;
            if (cand < best) best = cand;
        }

        grid_next[x][y][z].r2 = best;
        grid_next[x][y][z].r  = (best != INF_R2) ? isqrt((int)best) : 0;
    }

    grid_next[MID][MID][MID].r2 = 0;
    grid_next[MID][MID][MID].r  = 0;
}

/* ==========================================================
 * pulse_step — synchronous CA update + active shell
 * ========================================================== */
void pulse_step(void) {
    pulse_update_wavefront();

    unsigned int pulse_r2 = pulse_from_time((unsigned int)tick);

    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        unsigned int r2 = grid_next[x][y][z].r2;
        grid[x][y][z].r2 = r2;
        grid[x][y][z].r  = grid_next[x][y][z].r;

        if (r2 == INF_R2) {
            grid[x][y][z].active = 0;
        } else {
            unsigned int delta = (r2 > pulse_r2) ? (r2 - pulse_r2)
                                               : (pulse_r2 - r2);
            grid[x][y][z].active = (delta <= PULSE_TOLERANCE) ? 1 : 0;
        }
    }
}

/* ==========================================================
 * apply_move — commit lattice move m (additions/shifts only)
 *
 *   E   += G[m]                       (|c|^2 update)
 *   G[k]+= K[m][k]  for all k         (re-centre the increments)
 *   c   += D[m]                       (cross product update)
 *   q2  += (2 v_i + 1) or -(2 v_i - 1)   (sum of odds, |v|^2)
 * ========================================================== */
static void apply_move(int m) {
    E += G[m];
    for (int k = 0; k < 6; k++) G[k] += K[m][k];
    cxv += Dx[m];
    cyv += Dy[m];
    czv += Dz[m];

    switch (m) {
    case 0: q2 += ((long long)vx << 1) + 1; vx++; tip_x++; break;
    case 1: q2 -= ((long long)vx << 1) - 1; vx--; tip_x--; break;
    case 2: q2 += ((long long)vy << 1) + 1; vy++; tip_y++; break;
    case 3: q2 -= ((long long)vy << 1) - 1; vy--; tip_y--; break;
    case 4: q2 += ((long long)vz << 1) + 1; vz++; tip_z++; break;
    case 5: q2 -= ((long long)vz << 1) - 1; vz--; tip_z--; break;
    }
}

/* in-bounds and not already part of the arm */
static int move_ok(int m) {
    int nx = tip_x + mvx[m];
    int ny = tip_y + mvy[m];
    int nz = tip_z + mvz[m];
    if (nx < 0 || nx >= L || ny < 0 || ny >= L || nz < 0 || nz >= L) return 0;
    return !grid[nx][ny][nz].spin;
}

/* ==========================================================
 * spiral_step — advance the helix walker by one cell
 *
 * Bresenham: climb_total climb steps distributed among
 * orbit_total orbital steps (period = climb_total + orbit_total).
 *
 * Orbital step: among the 6 lattice moves pick the one minimising
 *   | (E + G[m]) - TGT |      (radial error on the cylinder)
 * restricted to moves with positive tangent projection
 *   t . e_m,  t = AXIS x u = -c    =>   t_+x = -c_x, t_-x = +c_x, ...
 * which enforces a single rotation sense (CCW around AXIS).
 *
 * Climb step: multi-dimensional DDA along AXIS —
 *   dda[i] += |a_i| ; take the largest ; dda[i] -= |a_x|+|a_y|+|a_z|
 * so lattice moves are distributed exactly in proportion to the
 * axis components, tracing the tilted axis direction.
 * ========================================================== */
void spiral_step(void) {
    if (spiral_done) return;

    int new_acc = tip_acc + climb_total;
    int do_climb = 0;
    if (new_acc >= period_total) { new_acc -= period_total; do_climb = 1; }
    tip_acc = new_acc;

    /* precompute legality once per step */
    int ok[6];
    for (int m = 0; m < 6; m++) ok[m] = move_ok(m);

    int chosen = -1;

    if (do_climb) {
        /* ---- climb step: DDA along AXIS ---- */
        dda[0] += absa[0];
        dda[1] += absa[1];
        dda[2] += absa[2];

        int i = 0;
        if (dda[1] > dda[i]) i = 1;
        if (dda[2] > dda[i]) i = 2;
        dda[i] -= climb_total;

        int sgn_neg = (i == 0) ? (axis_x < 0) : (i == 1) ? (axis_y < 0)
                                                          : (axis_z < 0);
        int m = (i << 1) + (sgn_neg ? 1 : 0);
        if (ok[m]) chosen = m;
        /* if blocked, fall through to an orbital choice below */
    }

    if (chosen < 0) {
        /* ---- orbital step: stay on the cylinder, keep turning ---- */
        long long t[6];
        t[0] = -cxv; t[1] =  cxv;
        t[2] = -cyv; t[3] =  cyv;
        t[4] = -czv; t[5] =  czv;

        long long best_err = 0, best_tan = 0;
        int best = -1;

        /* Preferred rule: among the moves that keep the cylinder radius
         * inside the tolerance band (|E' - TGT| <= TOL), take the one with
         * the largest tangential projection — maximal angular progress.
         * Ties and out-of-band situations fall back to minimum radial
         * error.  Only comparisons and additions are involved. */
        for (int m = 0; m < 6; m++) {
            if (!ok[m]) continue;
            if (t[m] <= 0) continue;                 /* wrong rotation sense */
            long long e = E + G[m] - TGT;
            if (e < 0) e = -e;
            if (e <= TOL) {
                if (best < 0 || best_err > TOL || t[m] > best_tan) {
                    best = m; best_err = e; best_tan = t[m];
                }
            } else if (best < 0 || (best_err > TOL && e < best_err)) {
                best = m; best_err = e; best_tan = t[m];
            }
        }

        if (best < 0) {                              /* fallback: any move */
            for (int m = 0; m < 6; m++) {
                if (!ok[m]) continue;
                long long e = E + G[m] - TGT;
                if (e < 0) e = -e;
                if (best < 0 || e < best_err) { best = m; best_err = e; }
            }
        }
        if (best < 0) { spiral_done = 1; return; }   /* fully boxed in */
        chosen = best;
    }

    apply_move(chosen);

    grid[tip_x][tip_y][tip_z].spin = 1;
    if (spiral_n < MAX_SPIRAL_PTS) {
        spiral_pts[spiral_n].x = tip_x;
        spiral_pts[spiral_n].y = tip_y;
        spiral_pts[spiral_n].z = tip_z;
        spiral_n++;
    }

    /* theta = pi reached: the arm touched the outer shell r = L/2 */
    if (q2 >= RADIUS_SQ || spiral_n >= MAX_SPIRAL_PTS) spiral_done = 1;
}

/* ==========================================================
 * Unified step — pulsating wavefront, then spiral advance
 * ========================================================== */
void step_all(void) {
    pulse_step();
    spiral_step();
    tick++;
}

/* ==========================================================
 * Rendering (SDL3)
 * ========================================================== */
#ifndef NO_SDL

void render_frame(SDL_Renderer *ren) {
    SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
    SDL_RenderClear(ren);

    unsigned int pulse_r2 = pulse_from_time((unsigned int)tick);
    int cur_r = isqrt((int)pulse_r2);

    /* -------------------------------------------------------
     * Panel 1 (top-left): wavefront z=MID slice
     *   green = distance gradient, yellow = active shell
     * ------------------------------------------------------- */
    {
        int ox = 30;
        for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            Cell *c = &grid[x][y][MID];
            unsigned int pr = 0, pg = 0, pb = 0;

            if (c->r2 != INF_R2) {
                int rr = isqrt((int)c->r2);
                int g = (rr < RADIUS) ? 255 - (rr * 255 / RADIUS) : 0;
                pg = (unsigned int)g;
            }

            {
                int delta = (int)c->r2 - (int)pulse_r2;
                if (delta < 0) delta = -delta;
                if (delta <= YELLOW_VIS_TOL) {
                    pr = 255; pg = 255; pb = 0;
                }
            }

            SDL_SetRenderDrawColor(ren,
                (Uint8)pr, (Uint8)pg, (Uint8)pb, 255);
            SDL_RenderPoint(ren, (float)(x + ox), (float)(y + 10));
        }
    }

    /* -------------------------------------------------------
     * Panel 2 (top-right, large): 3D isometric view of spiral
     *   Uses stored spiral_pts[] — no Bresenham re-trace.
     *   cyan = spin, white = spin+active.
     * ------------------------------------------------------- */
    {
        int panel_x = 30 + L + 40;
        int panel_y = 10;
        int panel_w = 480;
        int panel_h = L;
        float cx = (float)(panel_x + panel_w / 2);
        float cy = (float)(panel_y + panel_h / 2);
        float scale = (float)panel_h / (3.2f * (float)R_MAX);

        float ax = 0.866f;   /* cos(30°) */
        float ay = 0.5f;     /* sin(30°) */
        float ez = 1.0f;

        /* Render spiral points */
        for (int i = 0; i < spiral_n; i++) {
            int sx = spiral_pts[i].x;
            int sy = spiral_pts[i].y;
            int sz = spiral_pts[i].z;
            int ldx = sx - MID;
            int ldy = sy - MID;
            int dz  = sz - MID;

            float px = cx + ((float)ldx - (float)ldy) * ax * scale;
            float py = cy - (float)dz * ez * scale
                     + ((float)ldx + (float)ldy) * ay * scale;

            if (grid[sx][sy][sz].active) {
                SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
            } else {
                int bright = 150 + (ldx + ldy) / 4;
                if (bright > 255) bright = 255;
                if (bright < 80) bright = 80;
                SDL_SetRenderDrawColor(ren,
                    0, (Uint8)bright, (Uint8)bright, 255);
            }
            SDL_FRect rc = { px - 1, py - 1, 3, 3 };
            SDL_RenderFillRect(ren, &rc);
        }

        /* --- Coordinate axes (RGB = XYZ) --- */
        {
            float axis_len = (float)RADIUS * scale;

            /* X axis (red) */
            float x_end_px = cx + axis_len * ax;
            float x_end_py = cy + axis_len * ay;
            SDL_SetRenderDrawColor(ren, 180, 50, 50, 255);
            SDL_RenderLine(ren, cx, cy, x_end_px, x_end_py);

            /* Y axis (green) */
            float y_end_px = cx - axis_len * ax;
            float y_end_py = cy + axis_len * ay;
            SDL_SetRenderDrawColor(ren, 50, 180, 50, 255);
            SDL_RenderLine(ren, cx, cy, y_end_px, y_end_py);

            /* Z axis (blue up, dark blue down) */
            float z_end_py = cy - (float)R_MAX * ez * scale;
            SDL_SetRenderDrawColor(ren, 80, 80, 255, 255);
            SDL_RenderLine(ren, cx, cy, cx, z_end_py);
            float z_neg_py = cy + (float)R_MAX * ez * scale;
            SDL_SetRenderDrawColor(ren, 40, 40, 120, 255);
            SDL_RenderLine(ren, cx, cy, cx, z_neg_py);
        }

        /* --- Cosmetic sphere silhouette --- */
        {
            int n_seg = 64;
            float sil_r = (float)RADIUS * scale;
            float prev_px2 = 0.0f, prev_py2 = 0.0f;

            SDL_SetRenderDrawColor(ren, 35, 35, 35, 255);
            for (int i = 0; i <= n_seg; i++) {
                float angle = (float)i * 6.2832f / (float)n_seg;
                float ppx = cx + sil_r * cosf(angle);
                float ppy = cy + sil_r * sinf(angle);
                if (i > 0)
                    SDL_RenderLine(ren, prev_px2, prev_py2, ppx, ppy);
                prev_px2 = ppx;
                prev_py2 = ppy;
            }
        }
    }

    /* -------------------------------------------------------
     * Panels 3-5: Orthographic projections (XY, XZ, YZ)
     * ------------------------------------------------------- */
    {
        int tp_sz = 155;
        int tp_base_x = 30 + L + 40;
        int tp_base_y = L + 30;
        float tp_scale = (float)tp_sz / (2.2f * (float)RADIUS);

        int panel_idx;
        for (panel_idx = 0; panel_idx < 3; panel_idx++) {
            int tp_x = tp_base_x + panel_idx * (tp_sz + 5);
            float tp_cx = (float)tp_x + (float)tp_sz * 0.5f;
            float tp_cy = (float)tp_base_y + (float)tp_sz * 0.5f;

            /* Outer boundary circle */
            {
                int n_seg = 48;
                float prev_px2 = 0.0f, prev_py2 = 0.0f;
                SDL_SetRenderDrawColor(ren, 35, 35, 35, 255);
                for (int i = 0; i <= n_seg; i++) {
                    float angle = (float)i * 6.2832f / (float)n_seg;
                    float px2 = tp_cx + (float)RADIUS * cosf(angle) * tp_scale;
                    float py2 = tp_cy + (float)RADIUS * sinf(angle) * tp_scale;
                    if (i > 0)
                        SDL_RenderLine(ren, prev_px2, prev_py2, px2, py2);
                    prev_px2 = px2;
                    prev_py2 = py2;
                }
            }

            /* Wavefront circle (current radius) */
            {
                int n_seg = 48;
                float prev_px2 = 0.0f, prev_py2 = 0.0f;
                SDL_SetRenderDrawColor(ren, 50, 50, 0, 255);
                for (int i = 0; i <= n_seg; i++) {
                    float angle = (float)i * 6.2832f / (float)n_seg;
                    float px2 = tp_cx + (float)cur_r * cosf(angle) * tp_scale;
                    float py2 = tp_cy + (float)cur_r * sinf(angle) * tp_scale;
                    if (i > 0)
                        SDL_RenderLine(ren, prev_px2, prev_py2, px2, py2);
                    prev_px2 = px2;
                    prev_py2 = py2;
                }
            }

            /* Axes */
            {
                float axis_len = (float)RADIUS * tp_scale;
                if (panel_idx == 0) {
                    SDL_SetRenderDrawColor(ren, 140, 40, 40, 255);
                    SDL_RenderLine(ren, tp_cx, tp_cy, tp_cx + axis_len, tp_cy);
                    SDL_SetRenderDrawColor(ren, 40, 140, 40, 255);
                    SDL_RenderLine(ren, tp_cx, tp_cy, tp_cx, tp_cy + axis_len);
                } else if (panel_idx == 1) {
                    SDL_SetRenderDrawColor(ren, 140, 40, 40, 255);
                    SDL_RenderLine(ren, tp_cx, tp_cy, tp_cx + axis_len, tp_cy);
                    SDL_SetRenderDrawColor(ren, 40, 40, 180, 255);
                    SDL_RenderLine(ren, tp_cx, tp_cy, tp_cx, tp_cy - axis_len);
                } else {
                    SDL_SetRenderDrawColor(ren, 40, 140, 40, 255);
                    SDL_RenderLine(ren, tp_cx, tp_cy, tp_cx + axis_len, tp_cy);
                    SDL_SetRenderDrawColor(ren, 40, 40, 180, 255);
                    SDL_RenderLine(ren, tp_cx, tp_cy, tp_cx, tp_cy - axis_len);
                }
            }

            /* Spiral points from stored array */
            for (int i = 0; i < spiral_n; i++) {
                int spx = spiral_pts[i].x - MID;
                int spy = spiral_pts[i].y - MID;
                int spz = spiral_pts[i].z - MID;
                float scr_x, scr_y;

                if (panel_idx == 0) {
                    scr_x = tp_cx + (float)spx * tp_scale;
                    scr_y = tp_cy + (float)spy * tp_scale;
                } else if (panel_idx == 1) {
                    scr_x = tp_cx + (float)spx * tp_scale;
                    scr_y = tp_cy - (float)spz * tp_scale;
                } else {
                    scr_x = tp_cx + (float)spy * tp_scale;
                    scr_y = tp_cy - (float)spz * tp_scale;
                }

                if (grid[spiral_pts[i].x][spiral_pts[i].y][spiral_pts[i].z].active) {
                    SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
                    SDL_FRect rc = { scr_x - 2, scr_y - 2, 5, 5 };
                    SDL_RenderFillRect(ren, &rc);
                } else {
                    SDL_SetRenderDrawColor(ren, 0, 200, 200, 255);
                    SDL_FRect rc = { scr_x - 1, scr_y - 1, 3, 3 };
                    SDL_RenderFillRect(ren, &rc);
                }
            }
        }
    }

    /* -------------------------------------------------------
     * Bottom graph: spiral XY trajectory vs z
     *   red = (x - MID), green = (y - MID) for each spiral point
     * ------------------------------------------------------- */
    {
        int gx0 = 30;
        int gy0 = L + 50;
        int gw  = L;
        int gh  = 180;

        /* axis */
        SDL_SetRenderDrawColor(ren, 80, 80, 80, 255);
        SDL_RenderLine(ren, (float)gx0, (float)(gy0 + gh),
                       (float)(gx0 + gw), (float)(gy0 + gh));
        SDL_RenderLine(ren, (float)gx0, (float)gy0,
                       (float)gx0, (float)(gy0 + gh));
        /* zero line */
        SDL_SetRenderDrawColor(ren, 40, 40, 40, 255);
        SDL_RenderLine(ren, (float)gx0, (float)(gy0 + gh / 2),
                       (float)(gx0 + gw), (float)(gy0 + gh / 2));

        /* plot spiral points: x-displacement (red), y-displacement (green) */
        for (int i = 0; i < spiral_n; i++) {
            int z_idx = spiral_pts[i].z;
            int sx_off = spiral_pts[i].x - MID;
            int sy_off = spiral_pts[i].y - MID;

            int px = gx0 + z_idx;
            if (px < gx0 || px > gx0 + gw) continue;

            /* x displacement (red) */
            int py_x = gy0 + gh / 2 - (sx_off * gh / (R_MAX + R_MAX + 1));
            if (py_x < gy0) py_x = gy0;
            if (py_x > gy0 + gh) py_x = gy0 + gh;
            SDL_SetRenderDrawColor(ren, 255, 80, 80, 255);
            SDL_RenderPoint(ren, (float)px, (float)py_x);

            /* y displacement (green) */
            int py_y = gy0 + gh / 2 - (sy_off * gh / (R_MAX + R_MAX + 1));
            if (py_y < gy0) py_y = gy0;
            if (py_y > gy0 + gh) py_y = gy0 + gh;
            SDL_SetRenderDrawColor(ren, 80, 255, 80, 255);
            SDL_RenderPoint(ren, (float)px, (float)py_y);
        }

        /* vertical marker at current sweep z-levels */
        if (cur_r < R_MAX) {
            SDL_SetRenderDrawColor(ren, 60, 60, 60, 255);
            int z1 = MID + cur_r;
            int z2 = MID - cur_r;
            if (z1 < L) {
                SDL_RenderLine(ren,
                    (float)(gx0 + z1), (float)gy0,
                    (float)(gx0 + z1), (float)(gy0 + gh));
            }
            if (z2 >= 0) {
                SDL_RenderLine(ren,
                    (float)(gx0 + z2), (float)gy0,
                    (float)(gx0 + z2), (float)(gy0 + gh));
            }
        }
    }

    /* status line */
    {
        int spin_count = 0;
        for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            if (grid[x][y][MID].spin) spin_count++;
        }
        printf("\r[tick %4d] r=%3d  spiral=%3d  axis=(%d,%d,%d)  R_CYL=%d  ",
               tick, cur_r, spiral_n, axis_x, axis_y, axis_z, R_CYL);
        fflush(stdout);
    }

    SDL_RenderPresent(ren);
}

#endif /* !NO_SDL */

/* ==========================================================
 * Main (requires SDL)
 * ========================================================== */
#ifndef NO_SDL
int main(int argc, char **argv) {
    if (argc >= 4) {
        spiral_set_axis(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    } else {
        spiral_set_axis(AXIS_X, AXIS_Y, AXIS_Z);
    }

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        printf("SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow(
        "CA — cylindrical helix on pulsating wavefront",
        WINDOW_W, WINDOW_H, 0);
    if (!window) {
        printf("Window error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);
    if (!renderer) {
        printf("Renderer error: %s\n", SDL_GetError());
        SDL_DestroyWindow(window);
        SDL_Quit();
        return 1;
    }

    /* allocate grids on heap */
    grid      = malloc(sizeof(Cell) * L * L * L);
    grid_next = malloc(sizeof(Cell) * L * L * L);
    if (!grid || !grid_next) {
        printf("Out of memory (need ~%lu MB)\n",
               (unsigned long)(2 * sizeof(Cell) * L * L * L / (1024*1024)));
        return 1;
    }
    memset(grid,      0, sizeof(Cell) * L * L * L);
    memset(grid_next, 0, sizeof(Cell) * L * L * L);

    init();

    int running = 1;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            if (e.type == SDL_EVENT_QUIT) running = 0;
        }

        step_all();
        render_frame(renderer);
        SDL_Delay(16);
    }

    free(grid);
    free(grid_next);
    SDL_DestroyRenderer(renderer);
    SDL_DestroyWindow(window);
    SDL_Quit();
    return 0;
}
#endif /* !NO_SDL */