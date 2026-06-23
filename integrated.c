/*
 * integrated.c — sinc wave CA + cylindrical helix + pulsating wavefront
 *
 * Three mechanisms, one grid, AND triple detection:
 *   trig ∧ active ∧ spin
 *
 * sinc_step():  integer-only wave equation (no *, no /, no LUTs)
 * pulse_step(): BFS sum-of-odds + pulsating shell
 * generate_spiral(): cylindrical helix walker (init-time only)
 */

#include "integrated.h"
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <limits.h>

/* =================================================================
 * Global state
 * ================================================================= */
Cell (*grid)[L][L]      = NULL;
Cell (*grid_next)[L][L] = NULL;
int tick = 0;

SpiralPt spiral_pts[MAX_SPIRAL_PTS];
int      spiral_n = 0;

int and_count[L];   /* accumulated AND triple hits per shell radius */

/* =================================================================
 * Spiral geometry generator (walker, init-time only)
 *
 * Traces a cylindrical helix from origin (MID,MID,MID) to pole
 * (MID,MID,MID+R_MAX), marking spin=1. One full CCW revolution.
 *
 * Operations: addition, subtraction, shift, comparison.
 * ================================================================= */
void generate_spiral(void) {
    int tip_x  = MID;
    int tip_y  = MID;
    int tip_z  = MID;
    int tip_dx = MID - CYL_X0;
    int tip_dy = MID - CYL_Y0;
    int tip_d2 = tip_dx * tip_dx + tip_dy * tip_dy;
    int tip_acc = 0;

    grid[MID][MID][MID].spin = 1;
    spiral_pts[0].x = MID;
    spiral_pts[0].y = MID;
    spiral_pts[0].z = MID;
    spiral_n = 1;

    for (int step = 0; step < CYL_TOTAL + 50; step++) {
        int new_acc = tip_acc + R_MAX;
        int do_z;
        if (new_acc >= CYL_TOTAL) {
            new_acc -= CYL_TOTAL;
            do_z = 1;
        } else {
            do_z = 0;
        }
        tip_acc = new_acc;

        if (do_z) {
            tip_z++;
            if (tip_z >= L) break;
        } else {
            int d2_px = tip_d2 + (tip_dx << 1) + 1;
            int d2_mx = tip_d2 - (tip_dx << 1) + 1;
            int d2_py = tip_d2 + (tip_dy << 1) + 1;
            int d2_my = tip_d2 - (tip_dy << 1) + 1;

            int e_px = d2_px > R_CYL_SQ ? d2_px - R_CYL_SQ
                                        : R_CYL_SQ - d2_px;
            int e_mx = d2_mx > R_CYL_SQ ? d2_mx - R_CYL_SQ
                                        : R_CYL_SQ - d2_mx;
            int e_py = d2_py > R_CYL_SQ ? d2_py - R_CYL_SQ
                                        : R_CYL_SQ - d2_py;
            int e_my = d2_my > R_CYL_SQ ? d2_my - R_CYL_SQ
                                        : R_CYL_SQ - d2_my;

            int v_px = (tip_x + 1 < L)
                    && !grid[tip_x + 1][tip_y][tip_z].spin;
            int v_mx = (tip_x - 1 >= 0)
                    && !grid[tip_x - 1][tip_y][tip_z].spin;
            int v_py = (tip_y + 1 < L)
                    && !grid[tip_x][tip_y + 1][tip_z].spin;
            int v_my = (tip_y - 1 >= 0)
                    && !grid[tip_x][tip_y - 1][tip_z].spin;

            int t_px = -tip_dy;
            int t_mx =  tip_dy;
            int t_py =  tip_dx;
            int t_my = -tip_dx;

            int best = -1;
            int best_err = 0x7FFFFFFF;

            if (v_px && t_px > 0 && e_px < best_err)
                { best = 0; best_err = e_px; }
            if (v_mx && t_mx > 0 && e_mx < best_err)
                { best = 1; best_err = e_mx; }
            if (v_py && t_py > 0 && e_py < best_err)
                { best = 2; best_err = e_py; }
            if (v_my && t_my > 0 && e_my < best_err)
                { best = 3; best_err = e_my; }

            if (best < 0) {
                if (v_px && e_px < best_err)
                    { best = 0; best_err = e_px; }
                if (v_mx && e_mx < best_err)
                    { best = 1; best_err = e_mx; }
                if (v_py && e_py < best_err)
                    { best = 2; best_err = e_py; }
                if (v_my && e_my < best_err)
                    { best = 3; best_err = e_my; }
            }

            if (best < 0) {
                tip_z++;
                if (tip_z >= L) break;
            } else {
                switch (best) {
                case 0: tip_d2 += (tip_dx << 1) + 1;
                        tip_dx++; tip_x++; break;
                case 1: tip_d2 -= (tip_dx << 1) - 1;
                        tip_dx--; tip_x--; break;
                case 2: tip_d2 += (tip_dy << 1) + 1;
                        tip_dy++; tip_y++; break;
                case 3: tip_d2 -= (tip_dy << 1) - 1;
                        tip_dy--; tip_y--; break;
                }
            }
        }

        if (tip_x >= 0 && tip_x < L &&
            tip_y >= 0 && tip_y < L &&
            tip_z >= 0 && tip_z < L) {
            grid[tip_x][tip_y][tip_z].spin = 1;
            if (spiral_n < MAX_SPIRAL_PTS) {
                spiral_pts[spiral_n].x = tip_x;
                spiral_pts[spiral_n].y = tip_y;
                spiral_pts[spiral_n].z = tip_z;
                spiral_n++;
            }
        }

        if (tip_z >= MID + R_MAX) break;
    }
}

/* =================================================================
 * Initialization — sinc wave + wavefront seed + spiral geometry
 * ================================================================= */
void init(void) {
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        Cell *c   = &grid[x][y][z];
        c->u      = 0;
        c->v      = 0;
        c->acc    = 0;
        c->sinc_p = 0;
        c->sinc_q = 1;
        c->r      = 0;
        c->r2     = INF_R2;
        c->active = 0;
        c->ttl    = 0;
        c->trig   = 0;
        c->spin   = 0;
    }

    /* sinc wave seed */
    grid[MID][MID][MID].u  = 2048;
    /* BFS wavefront seed */
    grid[MID][MID][MID].r2 = 0;

    /* paint spiral geometry (static, init-time only) */
    generate_spiral();

    memset(and_count, 0, sizeof(and_count));
}

/* =================================================================
 * Sinc wave CA — one tick (integer-only in the hot loop)
 *
 * AND triple: triggered ∧ active ∧ spin
 * ================================================================= */
#ifndef USE_CUDA
void sinc_step(void) {
    int cur_sweep_r = isqrt((int)pulse_from_time((unsigned int)tick));

    for (int x = 1; x < L - 1; x++)
    for (int y = 1; y < L - 1; y++)
    for (int z = 1; z < L - 1; z++) {
        int u = grid[x][y][z].u;
        int v = grid[x][y][z].v;

        int neighbors =
            grid[x+1][y][z].u + grid[x-1][y][z].u +
            grid[x][y+1][z].u + grid[x][y-1][z].u +
            grid[x][y][z+1].u + grid[x][y][z-1].u;

        int lap = neighbors - (u << 2) - (u << 1);

        int r = grid[x][y][z].r;

        int diff_shift = DIFF_SHIFT + 1 - (r >> DIFF_DIV_SHIFT);
        if (diff_shift < DIFF_SHIFT - 1)
            diff_shift = DIFF_SHIFT - 1;

        int v_new = v + (lap >> diff_shift);
        int u_new = u + v_new;

        v_new -= (v_new >> VEL_DAMP_SHIFT);

        /* shell forcing */
        int dr = r - SHELL_R;
        if (dr < 0) dr = -dr;

        if (dr <= SHELL_W) {
            if (u > SHELL_TARGET) {
                int excess = u - SHELL_TARGET;
                v_new -= (excess >> 4);
            } else if ((tick & 3) == 0) {
                int deficit = SHELL_TARGET - u;
                v_new += (deficit >> 10) + 1;
            }
        }

        /* boundary absorption */
        if (r > RADIUS - ABSORB_W) {
            int dist = r - (RADIUS - ABSORB_W);
            if (dist >= ABSORB_W)
                u_new = 0;
            else
                u_new >>= dist;
        }
        if (r >= RADIUS)
            u_new = 0;
        if (u_new < 0)
            u_new = 0;

        /* Bresenham trigger */
        int acc = grid[x][y][z].acc + grid[x][y][z].sinc_p;
        int triggered = 0;
        if (acc >= grid[x][y][z].sinc_q &&
            grid[x][y][z].sinc_q > 0) {
            acc -= grid[x][y][z].sinc_q;
            triggered = 1;
        }

        /* TTL persistence */
        unsigned char ttl = grid[x][y][z].ttl;
        if ((tick & TTL_DECAY_MASK) == 0 && ttl > 0)
            ttl--;

        /* AND interaction: Bresenham trigger × pulsating active */
        if (triggered && grid[x][y][z].active)
        {
            ttl = 32 + ((223 * grid[x][y][z].sinc_p) /
                        grid[x][y][z].sinc_q);
            int rr = grid[x][y][z].r;
            if (rr >= 0 && rr < L && rr == cur_sweep_r)
                and_count[rr]++;
        }

        /* AND triple: trig ∧ active ∧ spin — separate visual marker */
        unsigned char ttl3 = grid[x][y][z].ttl_triple;
        if ((tick & TTL_DECAY_MASK) == 0 && ttl3 > 0)
            ttl3--;
        if (triggered && grid[x][y][z].active && grid[x][y][z].spin)
            ttl3 = 255;

        grid_next[x][y][z].u      = u_new;
        grid_next[x][y][z].v      = v_new;
        grid_next[x][y][z].acc    = acc;
        grid_next[x][y][z].sinc_p = grid[x][y][z].sinc_p;
        grid_next[x][y][z].sinc_q = grid[x][y][z].sinc_q;
        grid_next[x][y][z].ttl    = ttl;
        grid_next[x][y][z].ttl_triple = ttl3;
        grid_next[x][y][z].trig   = (unsigned char)triggered;
    }

    /* copy back with global damping */
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        grid[x][y][z].u =
            grid_next[x][y][z].u -
            (grid_next[x][y][z].u >> 12);
        grid[x][y][z].v =
            grid_next[x][y][z].v -
            (grid_next[x][y][z].v >> 12);
        grid[x][y][z].acc    = grid_next[x][y][z].acc;
        grid[x][y][z].sinc_p = grid_next[x][y][z].sinc_p;
        grid[x][y][z].sinc_q = grid_next[x][y][z].sinc_q;
        grid[x][y][z].ttl    = grid_next[x][y][z].ttl;
        grid[x][y][z].ttl_triple = grid_next[x][y][z].ttl_triple;
        grid[x][y][z].trig   = grid_next[x][y][z].trig;
    }
}

/* =================================================================
 * Pulsating wavefront CA — one tick
 * ================================================================= */
static const int ddx[6] = {1, -1, 0,  0, 0,  0};
static const int ddy[6] = {0,  0, 1, -1, 0,  0};
static const int ddz[6] = {0,  0, 0,  0, 1, -1};

static void pulse_update_wavefront(void) {
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++)
        grid_next[x][y][z].r2 = grid[x][y][z].r2;

    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        if (grid[x][y][z].r2 == INF_R2) continue;

        unsigned int ax = (x > MID) ? (unsigned int)(x - MID)
                                     : (unsigned int)(MID - x);
        unsigned int ay = (y > MID) ? (unsigned int)(y - MID)
                                     : (unsigned int)(MID - y);
        unsigned int az = (z > MID) ? (unsigned int)(z - MID)
                                     : (unsigned int)(MID - z);

        for (int d = 0; d < 6; d++) {
            int nx = x + ddx[d];
            int ny = y + ddy[d];
            int nz = z + ddz[d];
            if (nx < 0 || nx >= L || ny < 0 || ny >= L ||
                nz < 0 || nz >= L)
                continue;

            unsigned int diff;
            if (d < 2)      diff = (ax << 1) + 1;
            else if (d < 4) diff = (ay << 1) + 1;
            else            diff = (az << 1) + 1;

            unsigned int new_r2 = grid[x][y][z].r2 + diff;
            if (new_r2 < grid_next[nx][ny][nz].r2)
                grid_next[nx][ny][nz].r2 = new_r2;
        }
    }
}

void pulse_step(void) {
    pulse_update_wavefront();
    grid_next[MID][MID][MID].r2 = 0;

    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        unsigned int old_r2 = grid[x][y][z].r2;
        unsigned int new_r2 = grid_next[x][y][z].r2;
        grid[x][y][z].r2 = new_r2;
        if (new_r2 != INF_R2 && old_r2 == INF_R2)
            grid[x][y][z].r = isqrt((int)new_r2);
    }

    unsigned int pulse_r2 = pulse_from_time((unsigned int)tick);
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        unsigned int r2 = grid[x][y][z].r2;
        if (r2 == INF_R2) {
            grid[x][y][z].active = 0;
        } else {
            unsigned int delta = (r2 > pulse_r2)
                               ? (r2 - pulse_r2)
                               : (pulse_r2 - r2);
            grid[x][y][z].active = (delta <= PULSE_TOLERANCE) ? 1 : 0;
        }
    }
}
#endif /* !USE_CUDA */

/* =================================================================
 * Unified step — wavefront first, then sinc wave
 * ================================================================= */
void step_all(void) {
    pulse_step();
    sinc_step();
    tick++;
}

/* =================================================================
 * Rendering (SDL3) — 4 panels + sinc profile graph + AND triple dots
 * ================================================================= */
#ifndef NO_SDL

#define PEAK_HIST_W 600

static int64_t profile[L];
static int64_t prev_profile[L];
static int     rcount[L];
static int     peak_history[PEAK_HIST_W];
static int     peak_idx = 0;
static int     sinc_stable_frames = 0;
static int     sinc_converged = 0;
static int64_t u_peak = 0;

static void compute_profile(void) {
    for (int i = 0; i < L; i++) {
        prev_profile[i] = profile[i];
        profile[i] = 0;
        rcount[i] = 0;
    }
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        if (grid[x][y][z].r2 == INF_R2) continue;
        int r = grid[x][y][z].r;
        if (r < L) {
            profile[r] += grid[x][y][z].u;
            rcount[r]++;
        }
    }
    for (int i = 0; i < L; i++)
        if (rcount[i] > 0)
            profile[i] /= rcount[i];
}

static int64_t profile_max_change(void) {
    int64_t maxd = 0;
    for (int i = 0; i < RADIUS; i++) {
        int64_t d = profile[i] - prev_profile[i];
        if (d < 0) d = -d;
        if (d > maxd) maxd = d;
    }
    return maxd;
}

void render_frame(SDL_Renderer *ren) {
    SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
    SDL_RenderClear(ren);

    compute_profile();
    int64_t max_change = profile_max_change();

    /* --- Detect sinc convergence → set rationals --- */
    if (!sinc_converged) {
        if (max_change < STABILITY_THRESHOLD) sinc_stable_frames++;
        else sinc_stable_frames = 0;
        if (sinc_stable_frames >= STABILITY_FRAMES) {
            sinc_converged = 1;
            u_peak = 0;
            for (int r = 0; r < RADIUS; r++)
                if (profile[r] > u_peak) u_peak = profile[r];
            if (u_peak > 0) {
                for (int x = 0; x < L; x++)
                for (int y = 0; y < L; y++)
                for (int z = 0; z < L; z++) {
                    grid[x][y][z].sinc_p = grid[x][y][z].u;
                    grid[x][y][z].sinc_q = (int)u_peak;
                    grid[x][y][z].acc = 0;
                }
            }
            printf("\nSinc converged at tick %d — AND triple active.\n", tick);
        }
    }

    int64_t peak = 1;
    for (int r = 0; r < RADIUS; r++)
        if (rcount[r] > 0 && profile[r] > peak)
            peak = profile[r];
    peak_history[peak_idx] = (int)(peak > INT32_MAX ? INT32_MAX : peak);
    peak_idx = (peak_idx + 1) % PEAK_HIST_W;

    /* -------------------------------------------------------
     * Panel 1 (top-left): sinc displacement slice z=MID
     * ------------------------------------------------------- */
    {
        int slice_max = 1;
        for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            int val = grid[x][y][MID].u;
            if (val > slice_max) slice_max = val;
        }
        for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            int c = (grid[x][y][MID].u * 255) / slice_max;
            if (c > 255) c = 255;
            if (c < 0) c = 0;
            SDL_SetRenderDrawColor(ren, (Uint8)c, (Uint8)c, (Uint8)c, 255);
            SDL_RenderPoint(ren, (float)(x + 30), (float)(y + 10));
        }
    }

    /* -------------------------------------------------------
     * Panel 2 (top-middle): wavefront + spiral overlay z=MID
     *   green = distance, yellow = shell, cyan = spiral,
     *   white = AND triple hit
     * ------------------------------------------------------- */
    {
        int ox = 30 + L + 20;
        unsigned int pulse_thr = pulse_from_time((unsigned int)tick);

        for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            Cell *c = &grid[x][y][MID];
            uint32_t pix_r = 0, pix_g = 0, pix_b = 0;

            /* wavefront distance (green gradient) */
            if (c->r2 != INF_R2) {
                int rr = isqrt((int)c->r2);
                int g = (rr < RADIUS) ? 255 - (rr * 255 / RADIUS) : 0;
                pix_g = (uint32_t)g;
            }

            /* spiral overlay (cyan) */
            if (c->spin) {
                pix_r = 0; pix_g = 200; pix_b = 200;
            }

            /* wavefront shell (yellow ring) */
            {
                int delta = (int)c->r2 - (int)pulse_thr;
                if (delta < 0) delta = -delta;
                if (delta <= YELLOW_VIS_TOL) {
                    pix_r = 255; pix_g = 255; pix_b = 0;
                }
            }

            /* AND triple: trig ∧ active ∧ spin → white */
            if (c->trig && c->active && c->spin) {
                pix_r = 255; pix_g = 255; pix_b = 255;
            }

            SDL_SetRenderDrawColor(ren,
                (Uint8)pix_r, (Uint8)pix_g, (Uint8)pix_b, 255);
            SDL_RenderPoint(ren, (float)(x + ox), (float)(y + 10));
        }
    }

    /* -------------------------------------------------------
     * Panel 3 (top-right): TTL pattern z=MID
     *   orange/yellow = AND double (trig∧active)
     *   magenta       = AND triple (trig∧active∧spin)
     * ------------------------------------------------------- */
    {
        int ox = 30 + L + 20 + L + 20;
        for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            Cell *c = &grid[x][y][MID];
            uint32_t pix_r, pix_g, pix_b;
            if (c->ttl_triple > 0) {
                /* magenta for AND triple */
                pix_r = c->ttl_triple;
                pix_g = 0;
                pix_b = c->ttl_triple;
            } else {
                /* orange/yellow for AND double */
                pix_r = c->ttl;
                pix_g = c->ttl >> 1;
                pix_b = 0;
            }
            SDL_SetRenderDrawColor(ren,
                (Uint8)pix_r, (Uint8)pix_g, (Uint8)pix_b, 255);
            SDL_RenderPoint(ren, (float)(x + ox), (float)(y + 10));
        }
    }

    /* -------------------------------------------------------
     * Panel 4 (top-far-right): triggering cut z=MID
     *   cyan = trig, magenta = trig ∧ spin
     * ------------------------------------------------------- */
    {
        int ox = 30 + L + 20 + L + 20 + L + 20;
        for (int x = 0; x < L; x++)
        for (int y = 0; y < L; y++) {
            Cell *c = &grid[x][y][MID];
            uint32_t pix_r = 0, pix_g = 0, pix_b = 0;

            if (c->trig && c->spin) {
                pix_r = 255; pix_g = 0; pix_b = 255;  /* magenta */
            } else if (c->trig) {
                pix_r = 0; pix_g = 255; pix_b = 255;  /* cyan */
            }

            SDL_SetRenderDrawColor(ren,
                (Uint8)pix_r, (Uint8)pix_g, (Uint8)pix_b, 255);
            SDL_RenderPoint(ren, (float)(x + ox), (float)(y + 10));
        }
    }

    /* -------------------------------------------------------
     * Panel 5 (top, after panel 4): spiral XY projection
     *   All z-levels collapsed onto XY plane.
     *   cyan = spiral, white = spin∧active, yellow circle = pulse
     * ------------------------------------------------------- */
    {
        int ox = 30 + L + 20 + L + 20 + L + 20 + L + 20;
        int oy = 10;
        int half = L / 2;

        /* dark gray sphere boundary */
        SDL_SetRenderDrawColor(ren, 30, 30, 30, 255);
        for (int a = 0; a < 360; a++) {
            float rad = (float)a * 3.14159265f / 180.0f;
            float cx = half + RADIUS * cosf(rad);
            float cy = half + RADIUS * sinf(rad);
            SDL_RenderPoint(ren, (float)ox + cx, (float)oy + cy);
        }

        /* dark yellow: current pulse radius */
        {
            unsigned int pr2 = pulse_from_time((unsigned int)tick);
            int pr = isqrt((int)pr2);
            SDL_SetRenderDrawColor(ren, 80, 80, 0, 255);
            for (int a = 0; a < 360; a++) {
                float rad = (float)a * 3.14159265f / 180.0f;
                float cx = half + pr * cosf(rad);
                float cy = half + pr * sinf(rad);
                SDL_RenderPoint(ren, (float)ox + cx, (float)oy + cy);
            }
        }

        /* axes */
        SDL_SetRenderDrawColor(ren, 60, 0, 0, 255);  /* X red */
        SDL_RenderLine(ren, (float)(ox + half), (float)(oy + half),
                       (float)(ox + half + RADIUS), (float)(oy + half));
        SDL_SetRenderDrawColor(ren, 0, 60, 0, 255);  /* Y green */
        SDL_RenderLine(ren, (float)(ox + half), (float)(oy + half),
                       (float)(ox + half), (float)(oy + half + RADIUS));

        /* cyan: spiral points projected to XY */
        SDL_SetRenderDrawColor(ren, 0, 200, 200, 255);
        for (int i = 0; i < spiral_n; i++) {
            int px = spiral_pts[i].x - MID + half;
            int py = spiral_pts[i].y - MID + half;
            /* 3×3 for visibility */
            for (int dy = -1; dy <= 1; dy++)
            for (int dx = -1; dx <= 1; dx++)
                SDL_RenderPoint(ren,
                    (float)(ox + px + dx), (float)(oy + py + dy));
        }

        /* white: AND triple points (spin ∧ active ∧ trig) */
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        for (int i = 0; i < spiral_n; i++) {
            int sx = spiral_pts[i].x;
            int sy = spiral_pts[i].y;
            int sz = spiral_pts[i].z;
            Cell *c = &grid[sx][sy][sz];
            if (c->active && c->trig) {
                int px = sx - MID + half;
                int py = sy - MID + half;
                for (int dy = -2; dy <= 2; dy++)
                for (int dx = -2; dx <= 2; dx++)
                    SDL_RenderPoint(ren,
                        (float)(ox + px + dx), (float)(oy + py + dy));
            }
        }

        /* green: spin ∧ active (AND double) — slightly smaller */
        SDL_SetRenderDrawColor(ren, 0, 255, 0, 255);
        for (int i = 0; i < spiral_n; i++) {
            int sx = spiral_pts[i].x;
            int sy = spiral_pts[i].y;
            int sz = spiral_pts[i].z;
            Cell *c = &grid[sx][sy][sz];
            if (c->active && !c->trig) {
                int px = sx - MID + half;
                int py = sy - MID + half;
                for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++)
                    SDL_RenderPoint(ren,
                        (float)(ox + px + dx), (float)(oy + py + dy));
            }
        }

        /* cross at center */
        SDL_SetRenderDrawColor(ren, 100, 100, 100, 255);
        SDL_RenderLine(ren, (float)(ox + half - 3), (float)(oy + half),
                       (float)(ox + half + 3), (float)(oy + half));
        SDL_RenderLine(ren, (float)(ox + half), (float)(oy + half - 3),
                       (float)(ox + half), (float)(oy + half + 3));
    }

    /* -------------------------------------------------------
     * Bottom: sinc(r) profile + AND triple scatter
     * ------------------------------------------------------- */
    {
        int px0 = 80;
        int py0 = WINDOW_H - 40;
        int graph_w = RADIUS * GRAPH_SCALE_X;

        /* axis */
        SDL_SetRenderDrawColor(ren, 80, 80, 80, 255);
        SDL_RenderLine(ren, (float)px0, (float)py0,
                       (float)(px0 + graph_w), (float)py0);

        /* green: sinc(r) profile */
        SDL_SetRenderDrawColor(ren, 0, 255, 0, 255);
        {
            float gpx = -1, gpy = -1;
            for (int r = 0; r < RADIUS; r++) {
                if (rcount[r] > 0) {
                    float yf = (float)py0 -
                        ((float)profile[r] * GRAPH_HEIGHT) /
                        (float)PROFILE_PEAK_REF;
                    float xf = (float)px0 + (float)(r * GRAPH_SCALE_X);
                    if (gpx >= 0)
                        SDL_RenderLine(ren, gpx, gpy, xf, yf);
                    gpx = xf;
                    gpy = yf;
                }
            }
        }

        /* yellow: peak history */
        SDL_SetRenderDrawColor(ren, 255, 255, 0, 255);
        {
            int max_val = 1;
            for (int i = 0; i < PEAK_HIST_W; i++)
                if (peak_history[i] > max_val) max_val = peak_history[i];
            max_val += (max_val >> 3);

            float last_x = -1, last_y = -1;
            for (int i = 0; i < PEAK_HIST_W; i++) {
                int val = peak_history[i];
                if (val == 0) continue;
                float xf = (float)px0 +
                    ((float)i * (float)graph_w) / (float)PEAK_HIST_W;
                float yf = (float)py0 -
                    ((float)val * (float)GRAPH_HEIGHT) / (float)max_val;
                if (i == peak_idx) last_x = -1;
                if (last_x >= 0)
                    SDL_RenderLine(ren, last_x, last_y, xf, yf);
                last_x = xf;
                last_y = yf;
            }
        }

        /* cyan: trigger rate (after convergence) */
        if (sinc_converged) {
            SDL_SetRenderDrawColor(ren, 0, 200, 255, 255);
            float cpx = -1, cpy = -1;
            for (int r = 0; r < RADIUS; r++) {
                float rate = (u_peak > 0)
                    ? (float)profile[r] / (float)u_peak : 0;
                float yf = (float)py0 - rate * (float)GRAPH_HEIGHT;
                float xf = (float)px0 + (float)(r * GRAPH_SCALE_X);
                if (cpx >= 0)
                    SDL_RenderLine(ren, cpx, cpy, xf, yf);
                cpx = xf;
                cpy = yf;
            }
        }

        /* dark gray: current pulse radius */
        {
            unsigned int pr2 = pulse_from_time((unsigned int)tick);
            int cur_r = isqrt((int)pr2);
            if (cur_r < RADIUS) {
                int cx = px0 + cur_r * GRAPH_SCALE_X;
                SDL_SetRenderDrawColor(ren, 60, 60, 60, 255);
                SDL_RenderLine(ren, (float)cx,
                    (float)(py0 - GRAPH_HEIGHT), (float)cx, (float)py0);
            }
        }

        /* red: AND triple count per shell radius (scatter) */
        {
            int max_cnt = 0;
            for (int r = 0; r < RADIUS; r++)
                if (and_count[r] > max_cnt) max_cnt = and_count[r];
            if (max_cnt < 1) max_cnt = 1;

            SDL_SetRenderDrawColor(ren, 255, 60, 60, 255);
            for (int r = 0; r < RADIUS; r++) {
                if (and_count[r] == 0) continue;
                float yf = (float)py0 -
                    ((float)and_count[r] / (float)max_cnt) *
                    (float)GRAPH_HEIGHT;
                float xf = (float)px0 + (float)(r * GRAPH_SCALE_X);
                for (int dy = -1; dy <= 1; dy++)
                for (int dx = -1; dx <= 1; dx++)
                    SDL_RenderPoint(ren, xf + dx, yf + dy);
            }
        }

        /* status */
        {
            unsigned int pr2 = pulse_from_time((unsigned int)tick);
            int cr = isqrt((int)pr2);
            int total_and = 0;
            for (int r = 0; r < L; r++) total_and += and_count[r];
            printf("\r[tick %4d] peak=%lld stable=%d conv=%d r=%d spiral=%d AND=%d  ",
                   tick, (long long)peak, sinc_stable_frames,
                   sinc_converged, cr, spiral_n, total_and);
            fflush(stdout);
        }
    }

    SDL_RenderPresent(ren);
}
#endif /* !NO_SDL */

/* =================================================================
 * Main (SDL3)
 * ================================================================= */
#ifndef NO_SDL
int main(void) {
    if (!SDL_Init(SDL_INIT_VIDEO)) {
        printf("SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *window = SDL_CreateWindow(
        "CA — sinc(r) + spiral + pulsating wavefront — AND triple",
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
