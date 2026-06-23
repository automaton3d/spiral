/*
 * integrated.h — sinc wave CA + cylindrical helix + pulsating wavefront
 *
 * Three mechanisms share one grid:
 *   1) Sinc wave: 3D integer-only wave equation → emergent sin(r)/r
 *   2) Cylindrical helix: walker paints spiral geometry during init
 *   3) Pulsating wavefront: BFS distance propagation with pulsing shell
 *
 * AND triple detection: trig ∧ active ∧ spin
 *   trig   = Bresenham trigger from sinc wave (oscillating)
 *   active = on pulsating wavefront shell (sweeping)
 *   spin   = on spiral curve (static after init)
 *
 * Runtime CA operations: addition, subtraction, shift, comparison.
 * No multiplication, no division, no loops in cell rules.
 * (Multiplications in compile-time constants, init, and rendering only.)
 */

#ifndef INTEGRATED_H_
#define INTEGRATED_H_

/* --- Compiler portability --- */
#ifdef __CUDACC__
#define HD __host__ __device__
#else
#define HD
#endif

#if defined(_MSC_VER) && !defined(__cplusplus)
#define SINLINE static __inline
#else
#define SINLINE static inline
#endif

/* --- SDL (only for host rendering) --- */
#ifndef NO_SDL
#include <SDL3/SDL.h>
#endif

#include <stdint.h>

/* =================================================================
 * Grid dimensions — change L here to scale the entire simulation
 * ================================================================= */
#ifndef L
#define L 81
#endif

#define INF_R2       0xFFFFFFFFu
#define GRID_SPACING 8
#define MID          (L / 2)
#define R_MAX        (L / 2)

/* =================================================================
 * Geometry constants
 * ================================================================= */
#define RADIUS       (L / 2 - 2)
#define DIFF_SHIFT   4
#define SHELL_R      (L * 12 / 100)
#define SHELL_W      (L / 10)
#define CORE_R       3
#define SHELL_TARGET 16384

/* =================================================================
 * Scaled constants (L-invariant behaviour)
 * ================================================================= */
#define PULSE_TOLERANCE  1
#define SPIRAL_TOLERANCE 1
#define PULSE_STEP       (((L) + 15) / 30)
#define ABSORB_W         ((RADIUS / 27) > 2 ? (RADIUS / 27) : 2)
#define DIFF_DIV_SHIFT   ((RADIUS >= 384) ? 6 : (RADIUS >= 192) ? 5 : \
                          (RADIUS >=  96) ? 4 : (RADIUS >=  40) ? 3 : 2)
#define SMOOTH_W         ((RADIUS / 20) > 1 ? (RADIUS / 20) : 1)
#define YELLOW_VIS_TOL   ((RADIUS / 35) > 1 ? (RADIUS / 35) : 1)
#define TTL_DECAY_MASK   ((RADIUS >= 384) ? 127 : (RADIUS >= 192) ? 63 : \
                          (RADIUS >=  96) ?  31 : (RADIUS >=  40) ? 15 : 7)
#define VEL_DAMP_SHIFT   ((RADIUS >= 384) ? 9 : (RADIUS >= 192) ? 8 : \
                          (RADIUS >=  96) ? 7 : (RADIUS >=  40) ? 6 : 5)

/* =================================================================
 * Cylindrical helix parameters (compile-time constants)
 *
 * R_CYL/r_max ≈ 0.18103 (numerical fit, L-independent).
 * Integer approximation: 93/512 ≈ 0.18164 (<0.4% error).
 * Cylinder axis at (MID, MID + R_CYL), parallel to z.
 * ================================================================= */
#define R_CYL       ((R_MAX * 93 + 256) >> 9)
#define CYL_X0      MID
#define CYL_Y0      (MID + R_CYL)
#define R_CYL_SQ    (R_CYL * R_CYL)
#define CYL_CIRC    (R_CYL << 3)
#define CYL_TOTAL   (CYL_CIRC + R_MAX)

/* =================================================================
 * Display
 * ================================================================= */
#define WINDOW_W     1850
#define WINDOW_H     750
#define GRAPH_HEIGHT 400
#define GRAPH_SCALE_X (700 / (RADIUS > 1 ? RADIUS : 1))

/* =================================================================
 * Stability detection (for sinc convergence)
 * ================================================================= */
#define STABILITY_THRESHOLD 35
#define STABILITY_FRAMES    180

/* =================================================================
 * Graph normalization
 * ================================================================= */
#define PROFILE_PEAK_REF  (SHELL_TARGET * 3)

/* =================================================================
 * Flat indexing
 * ================================================================= */
#define IDX(x,y,z) ((x)*L*L + (y)*L + (z))

/* =================================================================
 * Unified Cell struct — sinc wave + spiral + wavefront
 * ================================================================= */
typedef struct {
    /* sinc wave CA */
    int u;                /* displacement */
    int v;                /* velocity */
    int acc;              /* Bresenham accumulator */
    int sinc_p;           /* emergent sinc numerator */
    int sinc_q;           /* emergent sinc denominator */
    /* shared geometry (filled by wavefront BFS) */
    int r;                /* integer radius from center */
    unsigned int r2;      /* Euclidean distance-squared (INF_R2 = unvisited) */
    unsigned int active;  /* 1 if on pulsating shell this tick */
    /* sinc trigger + persistence */
    unsigned char ttl;
    unsigned char ttl_triple; /* AND triple TTL (trig ∧ active ∧ spin) */
    unsigned char trig;   /* 1 if Bresenham triggered this tick */
    /* spiral geometry (static after init) */
    unsigned char spin;   /* 1 if on the cylindrical helix */
} Cell;

/* =================================================================
 * Spiral point storage
 * ================================================================= */
typedef struct { int x, y, z; } SpiralPt;

#define MAX_SPIRAL_PTS (R_MAX * 3)
extern SpiralPt spiral_pts[];
extern int      spiral_n;

/* =================================================================
 * Grid pointers
 * ================================================================= */
extern Cell (*grid)[L][L];
extern Cell (*grid_next)[L][L];
extern int tick;

/* AND hit counts */
extern int and_count[L];
extern int and_double_count;
extern int and_double_total;
extern int and_double_last;
extern int and_triple_count;
extern int and_triple_total;
extern int and_triple_last;

/* =================================================================
 * Integer square root (bit-by-bit, no multiplication)
 * ================================================================= */
SINLINE HD int isqrt(int n) {
    if (n <= 0) return 0;
    int result = 0;
    int bit = 1 << 30;
    while (bit > n) bit >>= 2;
    while (bit != 0) {
        if (n >= result + bit) {
            n -= result + bit;
            result = (result >> 1) + bit;
        } else {
            result >>= 1;
        }
        bit >>= 2;
    }
    return result;
}

/* =================================================================
 * Pulsating sweep: triangular wave in r2 space
 * ================================================================= */
SINLINE HD unsigned int pulse_from_time(unsigned int t) {
    const unsigned int min_r2 = 0;
    const unsigned int max_r2 =
        (unsigned int)((unsigned int)R_MAX * R_MAX * 92 / 100);
    const unsigned int step = PULSE_STEP;
    unsigned int span = max_r2 - min_r2;
    if (span == 0) return min_r2;
    unsigned int period = span + span;
    unsigned int phase  = (t * step) % period;
    if (phase < span)
        return min_r2 + phase;
    else
        return max_r2 - (phase - span);
}

/* =================================================================
 * Public API
 * ================================================================= */
void init(void);
void step_all(void);

#ifndef USE_CUDA
void sinc_step(void);
void pulse_step(void);
#endif

#ifndef NO_SDL
void render_frame(SDL_Renderer *ren);
#endif

/* spiral geometry generator (called by init) */
void generate_spiral(void);

#endif /* INTEGRATED_H_ */
