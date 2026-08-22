/*
 *  spiral_auto.c — THE FUSED AUTOMATON: a spiral CA that defines its
 *                  own direction from the w index.
 *
 *  One cellular automaton, two coupled layers, zero prescribed
 *  geometry:
 *
 *  1) ISOTROPIC BREATHING (no direction needed)
 *     An L^3 lattice diffuses r^2 by BFS sum-of-odds while a
 *     triangular pulse p(t) drives an active spherical shell.  This
 *     layer never consults any axis.
 *
 *  2) DIRECTION ELECTION AT THE EXPANSION LIMIT (ant view)
 *     Every cell carries a predefined value w in [0, 9L-1].  When
 *     the pulse hits its expansion limit (turn-around of the wave),
 *     exactly ONE cell of the outermost active shell is elected by a
 *     neighbour-to-neighbour tournament: payload = (score << 24) |
 *     code, score = hash((w + seed) mod 9L, code); each cell adopts
 *     the greatest payload among its six neighbours.  The strict
 *     total order makes the winner unique and scan-order
 *     independent; changing any w (or the seed, key W) elects a
 *     different surface point.  The winner BECOMES the momentum
 *     vector m — the rotation axis of layer 3.
 *
 *  3) THE HELIX (direction consumer)
 *     As soon as m exists, the classic Bresenham helix walker starts
 *     around the axis m (rescaled to |m| = L/2 by spiral_set_axis),
 *     tracing the spiral arm on the pulsating bubble.  The axis is
 *     frozen while the walker lives; the next election (next
 *     expansion limit) prepares a fresh direction for the next run.
 *
 *  Nothing is prescribed: the direction is harvested from the CA's
 *  own state at the exact instant the wave reaches its limit.
 *
 *  4) BROADCAST (ant view, whole-grid reach)
 *     Every newly computed spiral point stamps the current tick into
 *     its own cell of a second value grid; a copy-the-greater-
 *     neighbour diffusion then carries these monotonic values across
 *     the ENTIRE lattice through purely local moves.  Rendered as a
 *     faint 1-px dust: hue = which wave arrived last, brightness =
 *     how recently.  Lowest-priority layer, drawn first.
 *
 *  Controls:
 *    LMB drag   — orbit camera (yaw / pitch)
 *    RMB drag   — pan
 *    WHEEL      — zoom
 *    1 .. 5     — simulation steps per frame
 *    SPACE      — pause / resume
 *    A          — toggle auto-rotation
 *    W          — bump w seed and re-elect m immediately
 *    X          — toggle broadcast dust
 *    S          — toggle analytic wavefront sphere
 *    C          — toggle sampled active cells
 *    B          — toggle spiral arm voxels
 *    G          — toggle axes & bounding sphere
 *    R          — reset camera
 *    ESC        — quit
 */

#include "spiral.h"          /* CA core types, constants, pulse, isqrt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define TAU (6.283185307179586476925286766559f)

/* ============================================================
 * Camera state
 * ============================================================ */
static float cam_yaw   = -0.9f;
static float cam_pitch = -0.5f;
static float cam_dist  = 4.4f * (float)RADIUS;
static float pan_x = 0.0f, pan_y = 0.0f;

static bool auto_rotate  = false;
static bool show_axes    = true;
static bool show_sphere  = true;
static bool show_cells   = true;
static bool show_arm     = true;
static bool show_cast    = true;
static bool paused       = false;
static int  steps_pp     = 1;

/* ================================================================
 * Tiny 5x7 bitmap font (subset: A-Z, 0-9, space - . : /)
 * ================================================================ */
typedef struct { unsigned char ch; unsigned char b[7]; } Glyph;

static const Glyph FONT[] = {
    {'A',{0x0E,0x11,0x11,0x1F,0x11,0x11,0x11}},
    {'B',{0x1E,0x11,0x11,0x1E,0x11,0x11,0x1E}},
    {'C',{0x0E,0x11,0x10,0x10,0x10,0x11,0x0E}},
    {'D',{0x1E,0x11,0x11,0x11,0x11,0x11,0x1E}},
    {'E',{0x1F,0x10,0x10,0x1E,0x10,0x10,0x1F}},
    {'F',{0x1F,0x10,0x10,0x1E,0x10,0x10,0x10}},
    {'G',{0x0E,0x11,0x10,0x17,0x11,0x11,0x0F}},
    {'H',{0x11,0x11,0x11,0x1F,0x11,0x11,0x11}},
    {'I',{0x07,0x02,0x02,0x02,0x02,0x02,0x07}},
    {'J',{0x07,0x02,0x02,0x02,0x02,0x12,0x0C}},
    {'K',{0x11,0x12,0x14,0x18,0x14,0x12,0x11}},
    {'L',{0x10,0x10,0x10,0x10,0x10,0x10,0x1F}},
    {'M',{0x11,0x1B,0x15,0x15,0x11,0x11,0x11}},
    {'N',{0x11,0x19,0x15,0x13,0x11,0x11,0x11}},
    {'O',{0x0E,0x11,0x11,0x11,0x11,0x11,0x0E}},
    {'P',{0x1E,0x11,0x11,0x11,0x1E,0x10,0x10}},
    {'Q',{0x0E,0x11,0x11,0x11,0x15,0x12,0x0D}},
    {'R',{0x1E,0x11,0x11,0x11,0x1E,0x14,0x12}},
    {'S',{0x0F,0x10,0x10,0x0E,0x01,0x01,0x1E}},
    {'T',{0x1F,0x02,0x02,0x02,0x02,0x02,0x02}},
    {'U',{0x11,0x11,0x11,0x11,0x11,0x11,0x0E}},
    {'V',{0x11,0x11,0x11,0x11,0x11,0x0A,0x04}},
    {'W',{0x11,0x11,0x11,0x11,0x15,0x1B,0x11}},
    {'X',{0x11,0x11,0x0A,0x04,0x0A,0x11,0x11}},
    {'Y',{0x11,0x11,0x0A,0x04,0x04,0x04,0x04}},
    {'Z',{0x1F,0x01,0x02,0x04,0x08,0x10,0x1F}},
    {'0',{0x0E,0x11,0x13,0x15,0x19,0x11,0x0E}},
    {'1',{0x04,0x0C,0x04,0x04,0x04,0x04,0x0E}},
    {'2',{0x0E,0x11,0x01,0x02,0x04,0x08,0x1F}},
    {'3',{0x0E,0x11,0x01,0x06,0x01,0x11,0x0E}},
    {'4',{0x02,0x06,0x0A,0x12,0x1F,0x02,0x02}},
    {'5',{0x1F,0x10,0x1E,0x01,0x01,0x01,0x1E}},
    {'6',{0x06,0x08,0x10,0x1E,0x11,0x11,0x0E}},
    {'7',{0x1F,0x01,0x02,0x04,0x08,0x08,0x08}},
    {'8',{0x0E,0x11,0x11,0x0E,0x11,0x11,0x0E}},
    {'9',{0x0E,0x11,0x11,0x0F,0x01,0x02,0x0C}},
    {'.',{0x00,0x00,0x00,0x00,0x00,0x0C,0x0C}},
    {'-',{0x00,0x00,0x00,0x0F,0x00,0x00,0x00}},
    {':',{0x00,0x0C,0x0C,0x00,0x0C,0x0C,0x00}},
    {'/',{0x01,0x02,0x02,0x04,0x08,0x08,0x10}},
    {' ',{0x00,0x00,0x00,0x00,0x00,0x00,0x00}},
};

static void draw_char(SDL_Renderer *ren, char ch, int x, int y) {
    const Glyph *g = NULL;
    size_t i;
    for (i = 0; i < sizeof(FONT)/sizeof(FONT[0]); i++) {
        if (FONT[i].ch == (unsigned char)ch) { g = &FONT[i]; break; }
    }
    if (!g) { g = &FONT[sizeof(FONT)/sizeof(FONT[0]) - 1]; } /* space */

    for (int row = 0; row < 7; row++) {
        unsigned char b = g->b[row];
        for (int col = 0; col < 5; col++) {
            if (b & (0x10 >> col)) {
                SDL_FRect rc = { (float)(x + col * 2), (float)(y + row * 2), 2, 2 };
                SDL_RenderFillRect(ren, &rc);
            }
        }
    }
}

static void draw_text(SDL_Renderer *ren, const char *s, int x, int y) {
    while (*s) {
        draw_char(ren, *s, x, y);
        x += 12;
        s++;
    }
}

/* ================================================================
 * Projection helper (world coords centred on the lattice origin)
 * ================================================================ */
static void proj_point(float wx, float wy, float wz, int W, int H,
                       float *px, float *py, float *pd) {
    float cy = cosf(cam_yaw),   sy = sinf(cam_yaw);
    float cp = cosf(cam_pitch), sp = sinf(cam_pitch);

    float x1 =  cy * wx - sy * wz;
    float z1 =  sy * wx + cy * wz;
    float y2 =  cp * wy - sp * z1;
    float z2 =  sp * wy + cp * z1;

    float dz = cam_dist + z2;
    if (dz < 0.5f) dz = 0.5f;           /* near clip */
    float f = (float)H * 1.5f;
    float s = f / dz;

    *px = (float)W * 0.5f + pan_x + x1 * s;
    *py = (float)H * 0.5f - pan_y - y2 * s;
    *pd = dz;
}

/* ================================================================
 * Wireframe helpers
 * ================================================================ */
static void draw_circle(SDL_Renderer *ren, int W, int H,
                        float cx, float cy, float cz,
                        float r, Uint8 cr, Uint8 cg, Uint8 cb) {
    int nseg = 80;
    float prev_x = 0, prev_y = 0, prev_d = 0;
    SDL_SetRenderDrawColor(ren, cr, cg, cb, 180);

    for (int i = 0; i <= nseg; i++) {
        float a = TAU * (float)i / (float)nseg;
        float wx = cx + r * cosf(a);
        float wy = cy + r * sinf(a);
        float px, py, d;
        proj_point(wx, wy, cz, W, H, &px, &py, &d);
        if (i > 0 && d < cam_dist && prev_d < cam_dist) {
            SDL_RenderLine(ren, prev_x, prev_y, px, py);
        }
        prev_x = px; prev_y = py; prev_d = d;
    }
}

static void draw_sphere_frames(SDL_Renderer *ren, int W, int H, float R) {
    draw_circle(ren, W, H, 0, 0, 0, R, 60, 65, 80);
    draw_circle(ren, W, H, 0, 0, 0, R, 50, 55, 70);
    draw_circle(ren, W, H, 0, 0, 0, R, 55, 60, 75);
}

static void draw_world_axes(SDL_Renderer *ren, int W, int H, float R) {
    const int SEG = 64;
    Uint8 col[3][3] = {
        {220, 70, 70},   /* X red    */
        {70, 220, 70},   /* Y green  */
        {70, 120, 250},  /* Z blue   */
    };
    const float dirs[3][3] = {
        {1, 0, 0},
        {0, 1, 0},
        {0, 0, 1},
    };

    for (int a = 0; a < 3; a++) {
        float ex = dirs[a][0] * R, ey = dirs[a][1] * R, ez = dirs[a][2] * R;
        float px0, py0, pd0, px1, py1, pd1;
        proj_point(-ex, -ey, -ez, W, H, &px0, &py0, &pd0);
        proj_point( ex,  ey,  ez, W, H, &px1, &py1, &pd1);

        SDL_SetRenderDrawColor(ren, col[a][0], col[a][1], col[a][2], 200);
        SDL_RenderLine(ren, px0, py0, px1, py1);

        for (int i = 1; i < SEG; i++) {
            float s = -R + 2.0f * R * (float)i / (float)SEG;
            float wx = dirs[a][0] * s, wy = dirs[a][1] * s, wz = dirs[a][2] * s;
            float bx, by, bd;
            proj_point(wx, wy, wz, W, H, &bx, &by, &bd);
            SDL_FRect rc = { bx - 1.0f, by - 1.0f, 2.5f, 2.5f };
            SDL_RenderFillRect(ren, &rc);
        }
    }
}

/* ================================================================
 * Analytic pulsating sphere — Fibonacci-sphere point cloud drawn
 * at the nominal wavefront radius cur_r = isqrt(pulse(t)).
 * ================================================================ */
#define WAVE_PTS 3600

static float g_wav_depth[WAVE_PTS];

static int cmp_wav_depth(const void *a, const void *b) {
    int ia = *(const int *)a, ib = *(const int *)b;
    return (g_wav_depth[ia] < g_wav_depth[ib]) -
           (g_wav_depth[ia] > g_wav_depth[ib]);   /* far -> near */
}

static void draw_wave_sphere(SDL_Renderer *ren, int W, int H, float R) {
    static bool  seeded = false;
    static float ux[WAVE_PTS], uy[WAVE_PTS], uz[WAVE_PTS];
    static float sx[WAVE_PTS], sy[WAVE_PTS];
    static int   order[WAVE_PTS];

    if (!seeded) {
        const float GA = 3.14159265358979f * (3.0f - sqrtf(5.0f));
        for (int i = 0; i < WAVE_PTS; i++) {
            float t  = ((float)i + 0.5f) / (float)WAVE_PTS;
            uy[i]    = 1.0f - 2.0f * t;
            float rr = sqrtf(1.0f - uy[i] * uy[i]);
            float ph = GA * (float)i;
            ux[i]    = cosf(ph) * rr;
            uz[i]    = sinf(ph) * rr;
        }
        seeded = true;
    }

    if (R <= 0.0f) return;

    const float f   = (float)H * 1.5f;
    const float cyw = cosf(cam_yaw),   syw = sinf(cam_yaw);
    const float cpw = cosf(cam_pitch), spw = sinf(cam_pitch);

    for (int i = 0; i < WAVE_PTS; i++) {
        float wx = ux[i] * R, wy = uy[i] * R, wz = uz[i] * R;
        float x1 = cyw * wx - syw * wz;
        float z1 = syw * wx + cyw * wz;
        float y2 = cpw * wy - spw * z1;
        float z2 = spw * wy + cpw * z1;
        float dz = cam_dist + z2;
        if (dz < 0.5f) dz = 0.5f;
        g_wav_depth[i] = dz;
        float s = f / dz;
        sx[i] = (float)W * 0.5f + pan_x + x1 * s;
        sy[i] = (float)H * 0.5f - pan_y - y2 * s;
        order[i] = i;
    }
    qsort(order, WAVE_PTS, sizeof(int), cmp_wav_depth);

    const float RD = (float)RADIUS * 1.8f;

    for (int oi = 0; oi < WAVE_PTS; oi++) {
        int i = order[oi];
        float ft = (g_wav_depth[i] - (cam_dist - RD)) / (2.0f * RD);
        if (ft < 0.0f) ft = 0.0f; else if (ft > 1.0f) ft = 1.0f;
        Uint8 alpha = (Uint8)(205.0f - 135.0f * ft);
        SDL_SetRenderDrawColor(ren, 150, 235, 170, alpha);
        SDL_RenderPoint(ren, sx[i], sy[i]);
    }
}

/* ================================================================
 * Sampled active cells — sparse subset of the REAL grid.
 * ================================================================ */
#define CELL_STRIDE 16
#define C_N         (((L) + (CELL_STRIDE) - 1) / (CELL_STRIDE))
#define CELL_NSAMP  (C_N * C_N * C_N)

static int csamp_x[CELL_NSAMP], csamp_y[CELL_NSAMP], csamp_z[CELL_NSAMP];

static int draw_active_cells(SDL_Renderer *ren, int W, int H) {
    static bool seeded = false;
    if (!seeded) {
        int n = 0;
        for (int i = 0; i < C_N; i++)
        for (int j = 0; j < C_N; j++)
        for (int k = 0; k < C_N; k++) {
            csamp_x[n] = i * CELL_STRIDE;
            csamp_y[n] = j * CELL_STRIDE;
            csamp_z[n] = k * CELL_STRIDE;
            n++;
        }
        seeded = true;
    }

    int count = 0;
    for (int n = 0; n < CELL_NSAMP; n++) {
        int x = csamp_x[n], y = csamp_y[n], z = csamp_z[n];
        if (!grid[x][y][z].active) continue;
        count++;

        float px, py, pd;
        proj_point((float)(x - MID), (float)(y - MID),
                   (float)(z - MID), W, H, &px, &py, &pd);
        SDL_SetRenderDrawColor(ren, 255, 225, 110, 210);
        SDL_FRect rc = { px - 1.25f, py - 1.25f, 2.5f, 2.5f };
        SDL_RenderFillRect(ren, &rc);
    }
    return count;
}

/* ================================================================
 * ELECTION OF THE DIRECTION — ant-view tournament on w.
 *
 * Every cell carries a predefined value w in [0, 9L-1].  At the
 * expansion limit every ACTIVE shell cell is seeded with the
 * payload (score << 24) | code, score = hash((w + seed) mod 9L,
 * code); the maximum diffuses by "adopt the greatest neighbour".
 * code occupies the low bits, so the total order is strict: one
 * unique winner, independent of scan order.  The winner becomes
 * the momentum vector m — the spiral's rotation axis.
 * ================================================================ */
#define CODE_BITS 24
#define CODE_MASK ((1ull << CODE_BITS) - 1)
#define ELEC_PASSES_TOTAL    (2 * RADIUS + 8)
#define ELEC_PASSES_PER_TICK 8

static unsigned long long *bid = NULL;   /* per-cell payloads */
static unsigned int wseed  = 0;          /* global w dial     */
static int elec_left       = 0;          /* passes remaining  */
static int has_m           = 0;
static int elec_code       = -1;
static int mom_x = 0, mom_y = 0, mom_z = 0;

static unsigned int w_base_of(int x, int y, int z) {
    unsigned int h = (unsigned int)x * 73856093u
                   ^ (unsigned int)y * 19349663u
                   ^ (unsigned int)z * 83492791u;
    h ^= h >> 13; h *= 1274126177u; h ^= h >> 16;
    return h % (9u * (unsigned int)L);
}

static unsigned int elec_score(unsigned int weff, unsigned long long code) {
    unsigned int h = weff * 0x9E3779B1u;
    h ^= (unsigned int)(code >> 16);
    h *= 0x85EBCA6Bu;
    h ^= (unsigned int)(code & 0xFFFFu);
    h *= 0xC2B2AE35u;
    h ^= h >> 15;
    return h;
}

static void election_start(void) {
    long long i = 0;
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++, i++) {
        if (!grid[x][y][z].active) { bid[i] = 0; continue; }
        unsigned long long code = (unsigned long long)i;
        unsigned int weff =
            ((unsigned int)grid[x][y][z].w + wseed) % (9u * (unsigned int)L);
        bid[i] = ((unsigned long long)elec_score(weff, code) << CODE_BITS)
                 | code;
    }
    elec_left = ELEC_PASSES_TOTAL;
}

static void elec_diffuse_pass(int rev) {
    int x0 = rev ? L - 1 : 0, xs = rev ? -1 : 1;
    int y0 = rev ? L - 1 : 0, ys = rev ? -1 : 1;
    int z0 = rev ? L - 1 : 0, zs = rev ? -1 : 1;

    for (int xi = 0; xi < L; xi++) {
        int x = x0 + xi * xs;
        for (int yi = 0; yi < L; yi++) {
            int y = y0 + yi * ys;
            long long row = ((long long)x * L + y) * L;
            for (int zi = 0; zi < L; zi++) {
                int z = z0 + zi * zs;
                long long i = row + z;
                unsigned long long v = bid[i];
                unsigned long long m = v;
                if (x > 0)     { unsigned long long n = bid[i - (long long)L * L]; if (n > m) m = n; }
                if (x < L - 1) { unsigned long long n = bid[i + (long long)L * L]; if (n > m) m = n; }
                if (y > 0)     { unsigned long long n = bid[i - L]; if (n > m) m = n; }
                if (y < L - 1) { unsigned long long n = bid[i + L]; if (n > m) m = n; }
                if (z > 0)     { unsigned long long n = bid[i - 1]; if (n > m) m = n; }
                if (z < L - 1) { unsigned long long n = bid[i + 1]; if (n > m) m = n; }
                if (m != v) bid[i] = m;
            }
        }
    }
}

static int pending_apply = 0;   /* freshly elected axis waiting */

static void election_advance(void) {
    if (elec_left <= 0) return;

    int pass_idx = ELEC_PASSES_TOTAL - elec_left;
    int k = elec_left < ELEC_PASSES_PER_TICK ? elec_left
                                             : ELEC_PASSES_PER_TICK;
    for (int j = 0; j < k; j++)
        elec_diffuse_pass((pass_idx + j) & 1);
    elec_left -= k;

    if (elec_left == 0) {
        unsigned long long best = 0;
        long long bi = -1;
        const long long total = (long long)L * L * L;
        for (long long i = 0; i < total; i++) {
            if (bid[i] > best) { best = bid[i]; bi = i; }
        }
        if (bi >= 0) {
            int z = (int)(bi % L);
            int y = (int)((bi / L) % L);
            int x = (int)(bi / ((long long)L * L));
            mom_x = x - MID;
            mom_y = y - MID;
            mom_z = z - MID;
            elec_code = (int)(best & CODE_MASK);
            has_m = 1;
            pending_apply = 1;   /* hand the new direction over */
        }
    }
}

/* ================================================================
 * Orchestration: hand the elected direction to the walker.
 * The axis is frozen while a walk is alive; a pending election is
 * applied as soon as the current walk finishes.
 * ================================================================ */
static int walker_live = 0;

static void orchestrate(void) {
    if (!pending_apply) return;
    if (walker_live && !spiral_done) return;   /* axis frozen mid-walk */

    spiral_set_axis(mom_x, mom_y, mom_z);      /* normalise to |A|=L/2 */
    spiral_init();                             /* seed walker on m     */
    walker_live = 1;
    pending_apply = 0;
}

/* ================================================================
 * BROADCAST — ant-view diffusion covering the whole grid.
 *
 * Each newly computed spiral point stamps the current tick into its
 * own cell of bgrid (0 = never reached).  A copy-the-greatest-
 * neighbour pass then carries the monotonic values outward until
 * every cell of the lattice has been reached — global coverage by
 * purely local moves.  Visualised as faint 1-px dust over a sparse
 * sample: hue encodes WHICH value arrived last (golden-ratio walk,
 * different colour per wave), brightness how recently.
 * ================================================================ */
#define BCAST_EVERY   2      /* diffusion pass cadence (ticks)       */
#define BCAST_SSTRIDE 12     /* visual sampling stride               */
#define BG_N          (((L) + (BCAST_SSTRIDE) - 1) / (BCAST_SSTRIDE))
#define BCAST_NSAMP   (BG_N * BG_N * BG_N)
#define BCAST_YOUNG   36     /* ticks a fresh arrival stays glowing  */

static unsigned int *bgrid = NULL;   /* L^3 broadcast values, 0 = none */
static float bsamp_x[BCAST_NSAMP], bsamp_y[BCAST_NSAMP], bsamp_z[BCAST_NSAMP];

static void hsv2rgb(float h, float s, float v,
                    float *R, float *G, float *B) {
    float i = floorf(h * 6.0f);
    int   k = (int)i;
    if (k > 5) k = 5;
    if (k < 0) k = 0;
    float f = h * 6.0f - i;
    float p = v * (1.0f - s);
    float q = v * (1.0f - s * f);
    float t = v * (1.0f - s * (1.0f - f));
    switch (k) {
    case 0: *R = v; *G = t; *B = p; break;
    case 1: *R = q; *G = v; *B = p; break;
    case 2: *R = p; *G = v; *B = t; break;
    case 3: *R = p; *G = q; *B = v; break;
    case 4: *R = t; *G = p; *B = v; break;
    default:*R = v; *G = p; *B = q; break;
    }
}

static void bcast_inject(int x, int y, int z) {
    bgrid[((long long)x * L + y) * L + z] = (unsigned int)tick + 1u;
}

static void bcast_diffuse(void) {
    int p  = ((unsigned)tick / BCAST_EVERY) & 1;
    int x0 = p ? L - 1 : 0, xs = p ? -1 : 1;
    int y0 = p ? L - 1 : 0, ys = p ? -1 : 1;
    int z0 = p ? L - 1 : 0, zs = p ? -1 : 1;

    for (int xi = 0; xi < L; xi++) {
        int x = x0 + xi * xs;
        for (int yi = 0; yi < L; yi++) {
            int y = y0 + yi * ys;
            long long row = ((long long)x * L + y) * L;
            for (int zi = 0; zi < L; zi++) {
                int z = z0 + zi * zs;
                long long i = row + z;
                unsigned int v = bgrid[i];
                unsigned int m = v;
                if (x > 0)     { unsigned int n = bgrid[i - (long long)L * L]; if (n > m) m = n; }
                if (x < L - 1) { unsigned int n = bgrid[i + (long long)L * L]; if (n > m) m = n; }
                if (y > 0)     { unsigned int n = bgrid[i - L]; if (n > m) m = n; }
                if (y < L - 1) { unsigned int n = bgrid[i + L]; if (n > m) m = n; }
                if (z > 0)     { unsigned int n = bgrid[i - 1]; if (n > m) m = n; }
                if (z < L - 1) { unsigned int n = bgrid[i + 1]; if (n > m) m = n; }
                if (m != v) bgrid[i] = m;
            }
        }
    }
}

static void draw_broadcast(SDL_Renderer *ren, int W, int H) {
    static bool seeded = false;
    if (!seeded) {
        int n = 0;
        for (int i = 0; i < BG_N; i++)
        for (int j = 0; j < BG_N; j++)
        for (int k = 0; k < BG_N; k++) {
            bsamp_x[n] = (float)(i * BCAST_SSTRIDE) - MID;
            bsamp_y[n] = (float)(j * BCAST_SSTRIDE) - MID;
            bsamp_z[n] = (float)(k * BCAST_SSTRIDE) - MID;
            n++;
        }
        seeded = true;
    }

    const float RD = (float)RADIUS * 1.8f;

    for (int i = 0; i < BCAST_NSAMP; i++) {
        long long gx = (long long)(bsamp_x[i] + MID);
        long long gy = (long long)(bsamp_y[i] + MID);
        long long gz = (long long)(bsamp_z[i] + MID);
        unsigned int v = bgrid[(gx * L + gy) * L + gz];
        if (v == 0) continue;

        float px, py, pd;
        proj_point(bsamp_x[i], bsamp_y[i], bsamp_z[i], W, H, &px, &py, &pd);

        float h = fmodf((float)v * 0.61803398875f, 1.0f);
        float r, g, b;
        hsv2rgb(h, 0.55f, 1.0f, &r, &g, &b);

        int age = (int)tick - (int)v;
        float fr = 1.0f - (float)age / (float)BCAST_YOUNG;
        if (fr < 0.0f) fr = 0.0f;
        if (fr > 1.0f) fr = 1.0f;
        float base = 42.0f + 100.0f * fr;

        float ft = (pd - (cam_dist - RD)) / (2.0f * RD);
        if (ft < 0.0f) ft = 0.0f; else if (ft > 1.0f) ft = 1.0f;
        Uint8 a = (Uint8)(base * (0.65f + 0.35f * ft));

        SDL_SetRenderDrawColor(ren,
            (Uint8)(r * 255.0f), (Uint8)(g * 255.0f), (Uint8)(b * 255.0f), a);
        SDL_RenderPoint(ren, px, py);
    }
}

/* ================================================================
 * Spiral arm voxels (painter's algorithm)
 * ================================================================ */
static float g_vox_depth[MAX_SPIRAL_PTS];

static int cmp_vox_depth(const void *a, const void *b) {
    int ia = *(const int *)a, ib = *(const int *)b;
    return (g_vox_depth[ia] > g_vox_depth[ib]) -
           (g_vox_depth[ia] < g_vox_depth[ib]);
}

static void draw_arm(SDL_Renderer *ren, int W, int H) {
    int np = spiral_n;
    if (np <= 1) return;

    static int         vox_order [MAX_SPIRAL_PTS];
    static SDL_Vertex  vox_verts [MAX_SPIRAL_PTS * 12];
    static int         vox_idx   [MAX_SPIRAL_PTS * 18];

    float cy = cosf(cam_yaw),  sy = sinf(cam_yaw);
    float cp = cosf(cam_pitch), sp = sinf(cam_pitch);

    float nz[6] = {
         cp * sy,               /* +X */
        -cp * sy,               /* -X */
         sp,                    /* +Y */
        -sp,                    /* -Y */
         cp * cy,               /* +Z */
        -cp * cy,               /* -Z */
    };
    int vis[6];
    for (int f = 0; f < 6; f++) vis[f] = (nz[f] > 0.02f);

    const float hs = 0.175f;    /* half-size: visible gaps */
    const float cc[8][3] = {
        {-hs,-hs,-hs},{ hs,-hs,-hs},{-hs, hs,-hs},{ hs, hs,-hs},
        {-hs,-hs, hs},{ hs,-hs, hs},{-hs, hs, hs},{ hs, hs, hs},
    };
    static const int fq[6][4] = {
        {1,3,7,5},  /* +X */
        {0,4,6,2},  /* -X */
        {2,6,7,3},  /* +Y */
        {0,1,5,4},  /* -Y */
        {4,5,7,6},  /* +Z */
        {0,2,3,1},  /* -Z */
    };

    for (int i = 0; i < np; i++) {
        float wx = (float)(spiral_pts[i].x - MID);
        float wy = (float)(spiral_pts[i].y - MID);
        float wz = (float)(spiral_pts[i].z - MID);
        float z1 = sy * wx + cy * wz;
        g_vox_depth[i] = sp * wy + cp * z1;
        vox_order[i] = i;
    }
    qsort(vox_order, (size_t)np, sizeof(int), cmp_vox_depth);

    const float RD = (float)RADIUS * 1.8f;
    int nv = 0, ni = 0;

    for (int oi = 0; oi < np; oi++) {
        int i = vox_order[oi];
        float wx = (float)(spiral_pts[i].x - MID);
        float wy = (float)(spiral_pts[i].y - MID);
        float wz = (float)(spiral_pts[i].z - MID);

        float br, bg, bb;
        {
            float t = (float)i / (float)(np - 1);
            float h = t * 6.0f, x = 1.0f - fabsf(fmodf(h, 2.0f) - 1.0f);
            if (h < 1.0f)      { br = 1.0f; bg = x; bb = 0.0f; }
            else if (h < 2.0f) { br = x; bg = 1.0f; bb = 0.0f; }
            else if (h < 3.0f) { br = 0.0f; bg = 1.0f; bb = x; }
            else if (h < 4.0f) { br = 0.0f; bg = x; bb = 1.0f; }
            else if (h < 5.0f) { br = x; bg = 0.0f; bb = 1.0f; }
            else               { br = 1.0f; bg = 0.0f; bb = x; }
        }
        bool actv = grid[spiral_pts[i].x][spiral_pts[i].y][spiral_pts[i].z].active != 0;
        if (actv) { br = 1.0f; bg = 1.0f; bb = 1.0f; }

        float z2 = g_vox_depth[i];
        float ft = (z2 + RD) / (2.0f * RD);
        if (ft < 0.0f) ft = 0.0f; else if (ft > 1.0f) ft = 1.0f;
        float fade = 0.55f + 0.45f * ft;

        for (int f = 0; f < 6; f++) {
            if (!vis[f]) continue;
            float shade = 0.50f + 0.50f * nz[f];
            SDL_FColor col = { br * shade * fade,
                               bg * shade * fade,
                               bb * shade * fade,
                               1.0f };
            int base = nv;
            for (int k = 0; k < 4; k++) {
                int ci = fq[f][k];
                float pwx = wx + cc[ci][0];
                float pwy = wy + cc[ci][1];
                float pwz = wz + cc[ci][2];
                float z1v = sy * pwx + cy * pwz;
                float x1v = cy * pwx - sy * pwz;
                float y2v = cp * pwy - sp * z1v;
                float z2v = sp * pwy + cp * z1v;
                float dz  = cam_dist + z2v;
                if (dz < 0.5f) dz = 0.5f;
                float sc  = (float)H * 1.5f / dz;
                vox_verts[nv].position.x = (float)W * 0.5f + pan_x + x1v * sc;
                vox_verts[nv].position.y = (float)H * 0.5f - pan_y - y2v * sc;
                vox_verts[nv].color       = col;
                vox_verts[nv].tex_coord.x = 0.0f;
                vox_verts[nv].tex_coord.y = 0.0f;
                nv++;
            }
            vox_idx[ni++] = base + 0;
            vox_idx[ni++] = base + 1;
            vox_idx[ni++] = base + 2;
            vox_idx[ni++] = base + 0;
            vox_idx[ni++] = base + 2;
            vox_idx[ni++] = base + 3;
        }
    }

    if (nv > 0 && ni > 0)
        SDL_RenderGeometry(ren, NULL, vox_verts, nv, vox_idx, ni);
}

/* ================================================================
 * Main render
 * ================================================================ */
static void render(SDL_Renderer *ren, int W, int H) {
    SDL_SetRenderDrawColor(ren, 8, 8, 16, 255);
    SDL_RenderClear(ren);

    /* ---- lowest-priority layer first: the broadcast dust ---- */
    if (show_cast) draw_broadcast(ren, W, H);

    unsigned int pulse_r2 = pulse_from_time((unsigned int)tick);
    int cur_r = isqrt((int)pulse_r2);

    if (show_axes) {
        draw_sphere_frames(ren, W, H, (float)RADIUS);
        draw_world_axes(ren, W, H, (float)RADIUS);
    }
    if (show_sphere && cur_r > 0)
        draw_wave_sphere(ren, W, H, (float)cur_r);

    int act = 0;
    if (show_cells) act = draw_active_cells(ren, W, H);

    /* ---- the elected momentum vector m: centre -> winning cell --- */
    if (has_m) {
        float px, py, pd, qx, qy, qd;
        proj_point(0.0f, 0.0f, 0.0f, W, H, &px, &py, &pd);
        proj_point((float)mom_x, (float)mom_y, (float)mom_z,
                   W, H, &qx, &qy, &qd);
        SDL_SetRenderDrawColor(ren, 80, 255, 255, 230);
        SDL_RenderLine(ren, px, py, qx, qy);
        SDL_SetRenderDrawColor(ren, 255, 90, 255, 255);
        SDL_FRect rc = { qx - 3.0f, qy - 3.0f, 6.0f, 6.0f };
        SDL_RenderFillRect(ren, &rc);
        SDL_SetRenderDrawColor(ren, 255, 255, 255, 255);
        SDL_FRect rc2 = { qx - 1.5f, qy - 1.5f, 3.0f, 3.0f };
        SDL_RenderFillRect(ren, &rc2);
    }

    if (show_arm) draw_arm(ren, W, H);

    /* ------------------------------------------------------------
     * HUD
     * ------------------------------------------------------------ */
    {
        char l1[96], l2[96], l3[96], l4[96];

        SDL_SetRenderDrawColor(ren, 200, 200, 200, 255);
        draw_text(ren, "LMB ROTATE  RMB PAN  WHEEL ZOOM  SPACE PAUSE  1-5 SPEED  ESC QUIT",
                  10, 10);
        SDL_SetRenderDrawColor(ren, 170, 170, 170, 255);
        draw_text(ren, "A AUTO  W WSEED  X CAST  S SPHERE  C CELLS  B ARM  G AXES  R RESET",
                  10, 26);

        snprintf(l1, sizeof(l1),
                 "AXIS %d %d %d  M %d %d %d  NORM %.3f",
                 axis_x, axis_y, axis_z,
                 mom_x, mom_y, mom_z,
                 sqrt((double)mom_x * mom_x +
                      (double)mom_y * mom_y +
                      (double)mom_z * mom_z));
        snprintf(l2, sizeof(l2),
                 "WSEED %u  ACTIVE SAMPLES %d/%d",
                 wseed, act, CELL_NSAMP);
        if (elec_left > 0)
            snprintf(l3, sizeof(l3),
                     "ELECTING pass %d/%d",
                     ELEC_PASSES_TOTAL - elec_left, ELEC_PASSES_TOTAL);
        else if (has_m)
            snprintf(l3, sizeof(l3),
                     "ELECT code %d", elec_code);
        else
            snprintf(l3, sizeof(l3),
                     "(election at expansion peak)");
        snprintf(l4, sizeof(l4),
                 "SPIRAL %d/%d  R %d  TICK %d",
                 spiral_n, MAX_SPIRAL_PTS, cur_r, tick);

        SDL_SetRenderDrawColor(ren, 120, 220, 220, 255);
        draw_text(ren, l1, 10, 42);
        SDL_SetRenderDrawColor(ren, 180, 180, 180, 255);
        draw_text(ren, l2, 10, 58);
        SDL_SetRenderDrawColor(ren, 255, 140, 255, 255);
        draw_text(ren, l3, 10, 74);
        SDL_SetRenderDrawColor(ren, 120, 255, 160, 255);
        draw_text(ren, l4, 10, 90);

        if (paused) {
            SDL_SetRenderDrawColor(ren, 255, 220, 80, 255);
            draw_text(ren, "PAUSED", 10, 106);
        }
    }

    SDL_RenderPresent(ren);
}

/* ================================================================
 * Event handling
 * ================================================================ */
static void handle(SDL_Event *e) {
    static bool dragging_orbit = false;
    static bool dragging_pan   = false;
    static float last_x = 0, last_y = 0;

    switch (e->type) {
    case SDL_EVENT_QUIT:
        exit(0);
        break;

    case SDL_EVENT_MOUSE_BUTTON_DOWN:
        if (e->button.button == SDL_BUTTON_LEFT) {
            dragging_orbit = true;
            last_x = e->button.x;
            last_y = e->button.y;
        } else if (e->button.button == SDL_BUTTON_RIGHT) {
            dragging_pan = true;
            last_x = e->button.x;
            last_y = e->button.y;
        }
        break;

    case SDL_EVENT_MOUSE_BUTTON_UP:
        if (e->button.button == SDL_BUTTON_LEFT)  dragging_orbit = false;
        if (e->button.button == SDL_BUTTON_RIGHT) dragging_pan   = false;
        break;

    case SDL_EVENT_MOUSE_MOTION:
        if (dragging_orbit) {
            float dx = e->motion.x - last_x;
            float dy = e->motion.y - last_y;
            cam_yaw   -= dx * 0.008f;
            cam_pitch += dy * 0.008f;
            if (cam_pitch >  1.5f) cam_pitch =  1.5f;
            if (cam_pitch < -1.5f) cam_pitch = -1.5f;
            last_x = e->motion.x;
            last_y = e->motion.y;
        } else if (dragging_pan) {
            float scale = cam_dist / 600.0f;
            pan_x += (e->motion.x - last_x) * scale;
            pan_y += (e->motion.y - last_y) * scale;
            last_x = e->motion.x;
            last_y = e->motion.y;
        }
        break;

    case SDL_EVENT_MOUSE_WHEEL:
        cam_dist *= expf(-0.15f * e->wheel.y);
        if (cam_dist < RADIUS * 1.2f)    cam_dist = RADIUS * 1.2f;
        if (cam_dist > RADIUS * 12.0f)   cam_dist = RADIUS * 12.0f;
        break;

    case SDL_EVENT_KEY_DOWN:
        switch (e->key.key) {
        case SDLK_ESCAPE: exit(0); break;
        case SDLK_SPACE:  paused = !paused; break;
        case SDLK_A:      auto_rotate = !auto_rotate; break;
        case SDLK_W:      wseed++; election_start(); break;
        case SDLK_X:      show_cast = !show_cast; break;
        case SDLK_S:      show_sphere = !show_sphere; break;
        case SDLK_C:      show_cells = !show_cells; break;
        case SDLK_B:      show_arm = !show_arm; break;
        case SDLK_G:      show_axes = !show_axes; break;
        case SDLK_R:
            cam_yaw   = -0.9f;
            cam_pitch = -0.5f;
            cam_dist  = 4.4f * (float)RADIUS;
            pan_x = pan_y = 0.0f;
            break;
        case SDLK_1: steps_pp = 1; break;
        case SDLK_2: steps_pp = 2; break;
        case SDLK_3: steps_pp = 4; break;
        case SDLK_4: steps_pp = 8; break;
        case SDLK_5: steps_pp = 16; break;
        default: break;
        }
        break;

    default:
        break;
    }
}

/* ================================================================
 * main
 * ================================================================ */
int main(int argc, char **argv) {
    (void)argc; (void)argv;

    grid      = (Cell (*)[L][L])malloc(sizeof(Cell) * L * L * L);
    grid_next = (Cell (*)[L][L])malloc(sizeof(Cell) * L * L * L);
    if (!grid || !grid_next) {
        puts("Out of memory (need about 260 MB for the grids)");
        return 1;
    }
    memset(grid,      0, sizeof(Cell) * L * L * L);
    memset(grid_next, 0, sizeof(Cell) * L * L * L);

    /* election payload grid (ant-view tournament state) */
    bid = (unsigned long long *)calloc((size_t)L * L * L,
                                       sizeof(unsigned long long));
    if (!bid) {
        puts("Out of memory (election grid)");
        return 1;
    }

    /* broadcast value grid (ant-view state), 0 = never reached */
    bgrid = (unsigned int *)calloc((size_t)L * L * L, sizeof(unsigned int));
    if (!bgrid) {
        puts("Out of memory (broadcast grid)");
        return 1;
    }

    init();

    /* predefine the per-cell w values in [0, 9L-1] */
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        grid[x][y][z].w = (unsigned short)w_base_of(x, y, z);
    }

    /* hold the walker until the FIRST election delivers a direction:
     * pretend the initial (default-axis) walk already finished */
    spiral_done = 1;

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        puts("SDL_Init error");
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow(
        "Spiral CA — self-directed (axis elected from w at the expansion limit)",
        1024, 768, SDL_WINDOW_RESIZABLE);
    if (!win) {
        puts("Window error");
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *ren = SDL_CreateRenderer(win, NULL);
    if (!ren) {
        puts("Renderer error");
        SDL_DestroyWindow(win);
        SDL_Quit();
        return 1;
    }

    int running = 1;
    while (running) {
        SDL_Event e;
        while (SDL_PollEvent(&e)) {
            handle(&e);
            if (e.type == SDL_EVENT_QUIT) running = 0;
        }

        if (auto_rotate) cam_yaw += 0.02f;

        if (!paused) {
            Uint64 t0 = SDL_GetTicks();
            for (int s = 0; s < steps_pp; s++) {
                pulse_step();       /* isotropic breathing          */
                tick++;
                election_advance(); /* ant-view tournament progress */
                orchestrate();      /* hand elected m to the walker */

                /* expansion-limit detector: ascending ->
                 * descending phase crossing of the pulse */
                {
                    const unsigned int span =
                        (unsigned int)R_MAX * R_MAX * 92u / 100u;
                    const unsigned int per = span * 2u;
                    unsigned int pc = ((unsigned)tick * (unsigned)PULSE_STEP) % per;
                    unsigned int pp =
                        (((unsigned)tick - 1u) * (unsigned)PULSE_STEP) % per;
                    if ((pp < span) && (pc >= span))
                        election_start();
                }

                if (walker_live && !spiral_done) {
                    int before = spiral_n;
                    spiral_step();  /* walk on the elected axis */
                    /* stamp each fresh spiral point into bgrid */
                    while (before < spiral_n) {
                        bcast_inject(spiral_pts[before].x,
                                     spiral_pts[before].y,
                                     spiral_pts[before].z);
                        before++;
                    }
                }
                if ((tick % BCAST_EVERY) == 0)
                    bcast_diffuse();   /* ant-view broadcast pass */
            }
            Uint64 frame_ms = SDL_GetTicks() - t0;
            if (frame_ms < 13) SDL_Delay(13 - (Uint32)frame_ms);
        } else {
            SDL_Delay(13);
        }

        int W, H;
        SDL_GetWindowSize(win, &W, &H);
        if (W < 100) W = 100;
        if (H < 100) H = 100;
        render(ren, W, H);
    }

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();

    free(grid);
    free(grid_next);
    free(bid);
    free(bgrid);
    return 0;
}
