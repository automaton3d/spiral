/*
 *  momentum.c — standalone cellular automaton: the pulsating sphere,
 *               seeded with a NULL momentum vector m.
 *
 *  Same CA structure as the spiral family (spiral.h): an L^3 lattice
 *  whose cells carry the squared distance r^2 to the center, diffused
 *  by BFS with the classic "sum of odds" rule
 *
 *      (n+1)^2 - n^2 = 2n + 1
 *
 *  and driven by a triangular pulse p(t) in r^2 space, so a spherical
 *  shell of cells (|r^2 - p(t)| <= tolerance) breathes in and out.
 *
 *  MOMENTUM VECTOR m:
 *    The spiral programs use an arbitrary direction (the rotation
 *    axis) chosen up front.  Here that role belongs to the momentum
 *    vector m, which starts NULL — m = (0,0,0).
 *
 *    ELECTION AT THE EXPANSION LIMIT: every cell carries a
 *    predefined value w in [0, 9L-1].  When the pulse hits its
 *    expansion limit (turn-around of the triangular wave), exactly
 *    ONE point of the outermost active shell is elected and m takes
 *    its position.  The election is a neighbour-to-neighbour
 *    tournament (ant view): each cell adopts the greatest payload of
 *    its six neighbours, payload = (score << 24) | code, where score
 *    mixes the cell's w with its canonical index.  The strict total
 *    order makes the winner unique and scan-order independent, and
 *    changing any w (or the global seed, key W) elects a different
 *    surface point.
 *
 *  Visualisation: interactive 3D viewer showing only the pulsating
 *  sphere — a Fibonacci-sphere point cloud at the nominal wavefront
 *  radius, plus a sparse sample of the REAL active cells of the grid.
 *
 *  Controls:
 *    LMB drag   — orbit camera (yaw / pitch)
 *    RMB drag   — pan
 *    WHEEL      — zoom
 *    1 .. 5     — simulation steps per frame
 *    SPACE      — pause / resume
 *    A          — toggle auto-rotation
 *    W          — bump w seed and re-elect m immediately
 *    S          — toggle analytic wavefront sphere
 *    C          — toggle sampled active cells
 *    G          — toggle axes & bounding sphere
 *    R          — reset camera
 *    ESC        — quit
 */

#include "spiral.h"          /* same CA structure: L, Cell, pulse, isqrt */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define TAU (6.283185307179586476925286766559f)

/* ============================================================
 * Momentum vector m — born NULL, intentionally dormant.
 * ============================================================ */
static int mom_x = 0, mom_y = 0, mom_z = 0;

static double mom_norm(void) {
    return sqrt((double)mom_x * mom_x +
                (double)mom_y * mom_y +
                (double)mom_z * mom_z);
}

/* ---- predefined per-cell value w in [0, 9L-1] -------------------
 * Host-side hash of the position (init-time arithmetic may multiply;
 * the CA runtime only reads w). --------------------------------- */
static unsigned int w_base_of(int x, int y, int z) {
    unsigned int h = (unsigned int)x * 73856093u
                   ^ (unsigned int)y * 19349663u
                   ^ (unsigned int)z * 83492791u;
    h ^= h >> 13; h *= 1274126177u; h ^= h >> 16;
    return h % (9u * (unsigned int)L);
}

/* ---- election scoring -------------------------------------------
 * Mixes the effective w of a cell with its canonical index into a
 * 32-bit score; the full payload orders candidates lexicographically
 * by (score, code), which is a strict total order. --------------- */
#define CODE_BITS 24
#define CODE_MASK ((1ull << CODE_BITS) - 1)

static unsigned int elec_score(unsigned int weff, unsigned long long code) {
    unsigned int h = weff * 0x9E3779B1u;
    h ^= (unsigned int)(code >> 16);
    h *= 0x85EBCA6Bu;
    h ^= (unsigned int)(code & 0xFFFFu);
    h *= 0xC2B2AE35u;
    h ^= h >> 15;
    return h;
}

/* ============================================================
 * CA state (local names; the shared spiral core is not linked)
 * ============================================================ */
static Cell (*mgrid)[L][L]  = NULL;
static Cell (*mnext)[L][L]  = NULL;
static unsigned int mtick   = 0;

static void m_init(void) {
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        Cell *c  = &mgrid[x][y][z];
        c->r = 0; c->r2 = INF_R2; c->active = 0; c->spin = 0;
        c->w = (unsigned short)w_base_of(x, y, z);
        Cell *cn = &mnext[x][y][z];
        cn->r = 0; cn->r2 = INF_R2; cn->active = 0; cn->spin = 0;
    }
    /* seed the centre */
    mgrid[MID][MID][MID].r2 = 0;
    mgrid[MID][MID][MID].r  = 0;
}

/* ---- BFS wavefront: sum of odds, into the next buffer ---- */
static const int ddx[6] = {1, -1, 0,  0, 0,  0};
static const int ddy[6] = {0,  0, 1, -1, 0,  0};
static const int ddz[6] = {0,  0, 0,  0, 1, -1};

static void m_wavefront(void) {
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        mnext[x][y][z].r2 = mgrid[x][y][z].r2;
    }

    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        if (mgrid[x][y][z].r2 == INF_R2) continue;

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

            unsigned int new_r2 = mgrid[x][y][z].r2 + diff;
            if (new_r2 < mnext[nx][ny][nz].r2) {
                mnext[nx][ny][nz].r2 = new_r2;
            }
        }
    }
}

/* ---- one CA tick: BFS + copy-back + active shell flags ---- */
static void m_pulse_step(void) {
    m_wavefront();
    mnext[MID][MID][MID].r2 = 0;

    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        unsigned int old_r2 = mgrid[x][y][z].r2;
        unsigned int new_r2 = mnext[x][y][z].r2;
        mgrid[x][y][z].r2 = new_r2;
        if (new_r2 != INF_R2 && old_r2 == INF_R2) {
            mgrid[x][y][z].r = isqrt((int)new_r2);
        }
    }

    unsigned int pulse_r2 = pulse_from_time(mtick);
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++) {
        unsigned int r2 = mgrid[x][y][z].r2;
        if (r2 == INF_R2) {
            mgrid[x][y][z].active = 0;
        } else {
            unsigned int delta = (r2 > pulse_r2)
                               ? (r2 - pulse_r2)
                               : (pulse_r2 - r2);
            mgrid[x][y][z].active = (delta <= PULSE_TOLERANCE) ? 1 : 0;
        }
    }
}

/* ================================================================
 * Election of m at the expansion limit — ant-view tournament.
 *
 * State: one 64-bit payload per cell (bid); 0 = empty.  At the
 * expansion limit every ACTIVE shell cell is seeded with
 * (score << 24) | code; the maximum then diffuses through the
 * copy-the-greatest-neighbour rule until it has crossed the whole
 * shell (2*RADIUS+8 passes guarantee coverage).  The final maximum
 * is unique because code breaks every score tie.
 * ================================================================ */
#define ELEC_PASSES_TOTAL    (2 * RADIUS + 8)
#define ELEC_PASSES_PER_TICK 8

static unsigned long long *bid = NULL;   /* per-cell payloads */
static unsigned int wseed  = 0;          /* global w dial     */
static int elec_left       = 0;          /* passes remaining  */
static int has_m           = 0;
static int elec_code       = -1;

static void election_start(void) {
    long long i = 0;
    for (int x = 0; x < L; x++)
    for (int y = 0; y < L; y++)
    for (int z = 0; z < L; z++, i++) {
        if (!mgrid[x][y][z].active) { bid[i] = 0; continue; }
        unsigned long long code = (unsigned long long)i;
        unsigned int weff =
            ((unsigned int)mgrid[x][y][z].w + wseed) % (9u * (unsigned int)L);
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

static void election_advance(void) {
    if (elec_left <= 0) return;

    int pass_idx = ELEC_PASSES_TOTAL - elec_left;
    int k = elec_left < ELEC_PASSES_PER_TICK ? elec_left
                                             : ELEC_PASSES_PER_TICK;
    for (int j = 0; j < k; j++)
        elec_diffuse_pass((pass_idx + j) & 1);
    elec_left -= k;

    if (elec_left == 0) {
        /* read out the tournament winner (unique by construction) */
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
        }
    }
}

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
 * Sampled active cells — a sparse, fixed subset of the REAL grid
 * (every CELL_STRIDE-th site), proving the CA state visually.
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
        if (!mgrid[x][y][z].active) continue;
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
 * Main render
 * ================================================================ */
static void render(SDL_Renderer *ren, int W, int H) {
    SDL_SetRenderDrawColor(ren, 8, 8, 16, 255);
    SDL_RenderClear(ren);

    unsigned int pulse_r2 = pulse_from_time(mtick);
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

    /* ------------------------------------------------------------
     * HUD
     * ------------------------------------------------------------ */
    {
        char l1[96], l2[96], l3[96];

        SDL_SetRenderDrawColor(ren, 200, 200, 200, 255);
        draw_text(ren, "LMB ROTATE  RMB PAN  WHEEL ZOOM  SPACE PAUSE  1-5 SPEED  ESC QUIT",
                  10, 10);
        SDL_SetRenderDrawColor(ren, 170, 170, 170, 255);
        draw_text(ren, "A AUTO  W WSEED  S SPHERE  C CELLS  G AXES  R RESET",
                  10, 26);

        snprintf(l1, sizeof(l1),
                 "M %d %d %d  NORM %.3f  R %d  TICK %d",
                 mom_x, mom_y, mom_z, mom_norm(), cur_r, mtick);
        snprintf(l2, sizeof(l2),
                 "ACTIVE SAMPLES %d/%d  STEPS/FRAME %d",
                 act, CELL_NSAMP, steps_pp);
        if (elec_left > 0)
            snprintf(l3, sizeof(l3),
                     "ELECTING pass %d/%d  WSEED %u",
                     ELEC_PASSES_TOTAL - elec_left, ELEC_PASSES_TOTAL, wseed);
        else if (has_m)
            snprintf(l3, sizeof(l3),
                     "WSEED %u  ELECT code %d", wseed, elec_code);
        else
            snprintf(l3, sizeof(l3),
                     "WSEED %u  (election at expansion peak)", wseed);
        SDL_SetRenderDrawColor(ren, 120, 220, 220, 255);
        draw_text(ren, l1, 10, 42);
        SDL_SetRenderDrawColor(ren, 180, 180, 180, 255);
        draw_text(ren, l2, 10, 58);
        SDL_SetRenderDrawColor(ren, 255, 140, 255, 255);
        draw_text(ren, l3, 10, 74);

        if (paused) {
            SDL_SetRenderDrawColor(ren, 255, 220, 80, 255);
            draw_text(ren, "PAUSED", 10, 90);
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
        case SDLK_S:      show_sphere = !show_sphere; break;
        case SDLK_C:      show_cells = !show_cells; break;
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

    mgrid = (Cell (*)[L][L])malloc(sizeof(Cell) * L * L * L);
    mnext = (Cell (*)[L][L])malloc(sizeof(Cell) * L * L * L);
    if (!mgrid || !mnext) {
        puts("Out of memory (need about 260 MB for the grids)");
        return 1;
    }
    memset(mgrid, 0, sizeof(Cell) * L * L * L);
    memset(mnext, 0, sizeof(Cell) * L * L * L);

    /* election payload grid (ant-view tournament state) */
    bid = (unsigned long long *)calloc((size_t)L * L * L,
                                       sizeof(unsigned long long));
    if (!bid) {
        puts("Out of memory (election grid)");
        return 1;
    }

    m_init();

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        puts("SDL_Init error");
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow(
        "Momentum CA — pulsating sphere, m = 0 (LMB-orbit, RMB-pan, Wheel-zoom, ESC-quit)",
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
                m_pulse_step();
                mtick++;
                election_advance();   /* ant-view tournament progress */

                /* expansion-limit detector: the phase of the
                 * triangular pulse crosses from ascending to
                 * descending -> the outermost shell is present */
                {
                    const unsigned int span =
                        (unsigned int)R_MAX * R_MAX * 92u / 100u;
                    const unsigned int per = span * 2u;
                    unsigned int pc = (mtick * (unsigned)PULSE_STEP) % per;
                    unsigned int pp =
                        ((mtick - 1u) * (unsigned)PULSE_STEP) % per;
                    if ((pp < span) && (pc >= span))
                        election_start();
                }
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

    free(mgrid);
    free(mnext);
    free(bid);
    return 0;
}
