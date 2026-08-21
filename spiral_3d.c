/*
 *  spiral_3d.c — interactive 3D viewer for the spiral CA.
 *
 *  Uses the SAME cellular automaton as spiral.c (BFS sum-of-odds wavefront +
 *  Bresenham helix walker around an arbitrary axis).  The CA core is compiled
 *  separately with -DNO_SDL (see Makefile / Make.nmake / build_3d.bat);
 *  this file only adds the 3D rendering.
 *
 *  Rendering: the spiral arm is drawn as true 3D voxels (small cubes,
 *  smaller than one lattice cell so neighbours leave visible gaps), with
 *  every camera-facing face shaded by its normal, painter-sorted far→near.
 *
 *  Controls:
 *    LMB drag   — orbit camera (yaw / pitch)
 *    RMB drag   — pan
 *    WHEEL      — zoom
 *    1 .. 5     — simulation steps per frame
 *    SPACE      — pause / resume
 *    A          — toggle auto-rotation
 *    W          — toggle pulse-wave rings
 *    G          — toggle axes & bounding sphere
 *    R          — reset camera
 *    ESC        — quit
 */

#include "spiral.h"          /* includes SDL3/SDL.h for us (NO_SDL not set) */
#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <stdbool.h>

#define TAU (6.283185307179586476925286766559f)

/* ============================================================
 * Camera state
 * ============================================================ */
static float cam_yaw   = -0.9f;     /* rotation around world Y        */
static float cam_pitch = -0.5f;     /* rotation around world X        */
static float cam_dist  = 4.4f * (float)RADIUS;
static float pan_x = 0.0f, pan_y = 0.0f;

static bool auto_rotate  = false;
static bool show_wave    = true;
static bool show_axes    = true;
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
 * Projection helper (world coords already centred on lattice origin)
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
 * 3D wireframe helpers (drawn in world space, radius in grid units)
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

static void draw_wave_rings(SDL_Renderer *ren, int W, int H, int cur_r) {
    /* wavefront rings in the three coordinate planes */
    float r = (float)cur_r;
    draw_circle(ren, W, H, 0, 0, 0, r, 120, 255, 120);   /* XY  (green)   */
    draw_circle(ren, W, H, 0, 0, 0, r, 120, 180, 255);   /* XZ  (blue)    */
    draw_circle(ren, W, H, 0, 0, 0, r, 200, 140, 255);   /* YZ  (violet)  */
}

static void draw_sphere_frames(SDL_Renderer *ren, int W, int H, float R) {
    /* three great circles on the bounding sphere */
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

        /* tick marks */
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

static void draw_rotation_axis(SDL_Renderer *ren, int W, int H) {
    /* axis vector already has |AXIS| = R_MAX (L/2); scale to RADIUS */
    float an = sqrtf((float)axis_x * axis_x +
                     (float)axis_y * axis_y +
                     (float)axis_z * axis_z);
    if (an < 1e-6f) return;
    float s = (float)RADIUS / an;
    float ax = (float)axis_x * s;
    float ay = (float)axis_y * s;
    float az = (float)axis_z * s;

    float px0, py0, pd0, px1, py1, pd1;
    proj_point(-ax, -ay, -az, W, H, &px0, &py0, &pd0);
    proj_point( ax,  ay,  az, W, H, &px1, &py1, &pd1);

    SDL_SetRenderDrawColor(ren, 60, 220, 220, 255);
    SDL_RenderLine(ren, px0, py0, px1, py1);

    /* end caps */
    SDL_SetRenderDrawColor(ren, 100, 255, 255, 255);
    SDL_FRect rc0 = { px0 - 3.0f, py0 - 3.0f, 7.0f, 7.0f };
    SDL_FRect rc1 = { px1 - 3.0f, py1 - 3.0f, 7.0f, 7.0f };
    SDL_RenderFillRect(ren, &rc0);
    SDL_RenderFillRect(ren, &rc1);

    /* small centre marker */
    float cpx, cpy, cpd;
    proj_point(0, 0, 0, W, H, &cpx, &cpy, &cpd);
    SDL_SetRenderDrawColor(ren, 255, 255, 140, 255);
    SDL_FRect rc = { cpx - 1.5f, cpy - 1.5f, 3.0f, 3.0f };
    SDL_RenderFillRect(ren, &rc);
}

/* ================================================================
 * Voxel pipeline (painter's algorithm)
 * ================================================================ */
static float g_vox_depth[MAX_SPIRAL_PTS];

static int cmp_vox_depth(const void *a, const void *b) {
    int ia = *(const int *)a, ib = *(const int *)b;
    return (g_vox_depth[ia] > g_vox_depth[ib]) -
           (g_vox_depth[ia] < g_vox_depth[ib]);
}

/* ================================================================
 * Main render
 * ================================================================ */
static void render(SDL_Renderer *ren, int W, int H) {
    SDL_SetRenderDrawColor(ren, 8, 8, 16, 255);
    SDL_RenderClear(ren);

    unsigned int pulse_r2 = pulse_from_time((unsigned int)tick);
    int cur_r = isqrt((int)pulse_r2);

    if (show_axes) {
        draw_sphere_frames(ren, W, H, (float)RADIUS);
        draw_world_axes(ren, W, H, (float)RADIUS);
    }
    if (show_wave) draw_wave_rings(ren, W, H, cur_r);
    draw_rotation_axis(ren, W, H);

    /* ------------------------------------------------------------
     * The spiral itself: true 3D voxels (small cubes with every
     * camera-facing face shaded by its normal), painter-sorted
     * far→near, coloured by arc position; white when the cell is on
     * the currently active wavefront.
     * ------------------------------------------------------------ */
    {
        int np = spiral_n;
        if (np > 1) {
            static int         vox_order [MAX_SPIRAL_PTS];
            static SDL_Vertex  vox_verts [MAX_SPIRAL_PTS * 12];
            static int         vox_idx   [MAX_SPIRAL_PTS * 18];

            float cy = cosf(cam_yaw),  sy = sinf(cam_yaw);
            float cp = cosf(cam_pitch), sp = sinf(cam_pitch);

            /* view-space z-component of the 6 face normals;
             * positive = faces the camera (visible).  View matrix's
             * last row (viewer direction) is (cp*sy, sp, cp*cy). */
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

            const float hs = 0.35f;     /* <1 → visible gap between voxels */
            const float cc[8][3] = {
                {-hs,-hs,-hs},{ hs,-hs,-hs},{-hs, hs,-hs},{ hs, hs,-hs},
                {-hs,-hs, hs},{ hs,-hs, hs},{-hs, hs, hs},{ hs, hs, hs},
            };
            /* CCW-from-outside quads (indices into cc) */
            static const int fq[6][4] = {
                {1,3,7,5},  /* +X */
                {0,4,6,2},  /* -X */
                {2,6,7,3},  /* +Y */
                {0,1,5,4},  /* -Y */
                {4,5,7,6},  /* +Z */
                {0,2,3,1},  /* -Z */
            };

            /* view-space depth per voxel, far → near */
            for (int i = 0; i < np; i++) {
                float wx = (float)(spiral_pts[i].x - MID);
                float wy = (float)(spiral_pts[i].y - MID);
                float wz = (float)(spiral_pts[i].z - MID);
                float z1 = sy * wx + cy * wz;
                g_vox_depth[i] = sp * wy + cp * z1;
                vox_order[i] = i;
            }
            qsort(vox_order, (size_t)np, sizeof(int), cmp_vox_depth);

            const float RD = (float)RADIUS * 1.8f;  /* depth-fade range */
            int nv = 0, ni = 0;

            for (int oi = 0; oi < np; oi++) {
                int i = vox_order[oi];
                float wx = (float)(spiral_pts[i].x - MID);
                float wy = (float)(spiral_pts[i].y - MID);
                float wz = (float)(spiral_pts[i].z - MID);

                /* arc-position base colour (rainbow) */
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

                /* depth fade from the voxel centre */
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
    }

    /* ------------------------------------------------------------
     * HUD
     * ------------------------------------------------------------ */
    {
        char axis_txt[96];
        snprintf(axis_txt, sizeof(axis_txt),
                 "AXIS %d %d %d  R=%d  N=%d/%d  TICK=%d",
                 axis_x, axis_y, axis_z, cur_r, spiral_n, MAX_SPIRAL_PTS, tick);

        SDL_SetRenderDrawColor(ren, 200, 200, 200, 255);
        draw_text(ren, "LMB ROTATE  RMB PAN  WHEEL ZOOM  SPACE PAUSE  A AUTO  W WAVE  G AXES  R RESET  ESC QUIT",
                  10, 10);
        SDL_SetRenderDrawColor(ren, 120, 220, 220, 255);
        draw_text(ren, axis_txt, 10, 26);

        if (paused) {
            SDL_SetRenderDrawColor(ren, 255, 220, 80, 255);
            draw_text(ren, "PAUSED", 10, 42);
        }
        if (spiral_done) {
            SDL_SetRenderDrawColor(ren, 120, 255, 120, 255);
            draw_text(ren, "SPIRAL COMPLETE", 10, 42);
        }
        {
            char sp[24];
            snprintf(sp, sizeof(sp), "STEPS/FRAME %d", steps_pp);
            SDL_SetRenderDrawColor(ren, 180, 180, 180, 255);
            draw_text(ren, sp, W - 180, 10);
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
        case SDLK_W:      show_wave = !show_wave; break;
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
    if (argc >= 4) {
        spiral_set_axis(atoi(argv[1]), atoi(argv[2]), atoi(argv[3]));
    } else {
        spiral_set_axis(AXIS_X, AXIS_Y, AXIS_Z);
    }

    /* allocate grids (exactly as the headless test does) */
    grid      = (Cell (*)[L][L])malloc(sizeof(Cell) * L * L * L);
    grid_next = (Cell (*)[L][L])malloc(sizeof(Cell) * L * L * L);
    if (!grid || !grid_next) {
        printf("Out of memory (need ~%lu MB)\n",
               (unsigned long)(2 * sizeof(Cell) * L * L * L / (1024 * 1024)));
        return 1;
    }
    memset(grid,      0, sizeof(Cell) * L * L * L);
    memset(grid_next, 0, sizeof(Cell) * L * L * L);

    init();

    if (!SDL_Init(SDL_INIT_VIDEO)) {
        printf("SDL_Init error: %s\n", SDL_GetError());
        return 1;
    }

    SDL_Window *win = SDL_CreateWindow(
        "Spiral CA — interactive 3D viewer (LMB-orbit, RMB-pan, Wheel-zoom, ESC-quit)",
        1024, 768, SDL_WINDOW_RESIZABLE);
    if (!win) {
        printf("Window error: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *ren = SDL_CreateRenderer(win, NULL);
    if (!ren) {
        printf("Renderer error: %s\n", SDL_GetError());
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
                if (spiral_done) break;
                step_all();
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
    return 0;
}