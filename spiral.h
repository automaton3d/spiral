/*
 * spiral.h — header for emergent spiral CA on pulsating wavefront
 *            (arbitrary rotation axis version)
 *
 * The spiral is a cylindrical helix traced by a walker:
 *   1) BFS wavefront propagation (sum-of-odds for r2)
 *   2) Walker advances 1 cell/tick around a cylinder of radius R_CYL
 *      whose axis is parallel to the integer vector
 *      AXIS = (axis_x, axis_y, axis_z) with |AXIS| = L/2, and whose
 *      centre line passes at distance R_CYL from the lattice centre,
 *      so the walk starts at r = 0 and ends at r = L/2 (theta: 0..pi).
 *   3) A Bresenham accumulator interleaves orbital steps (rotation
 *      around AXIS) with climb steps (6-neighbour DDA along AXIS).
 *
 * Runtime CA operations: addition, subtraction, shift, comparison.
 * NO multiplication, NO division, NO floats, NO lookup of trig.
 * All axis-dependent constants (cross-product increments D[m], the
 * incremental table K[k][m] = 2 D_k . D_m, the offset vector P and
 * the target |r x AXIS|^2) are computed ONCE at init time, exactly
 * like the old 2*dx+1 constants were.
 *
 * Radial measure: for a point u = r - P (P = offset of the cylinder
 * axis line from the centre), the squared distance to the axis is
 *      dist^2 = |u x AXIS|^2 / |AXIS|^2
 * so the walker compares |u x AXIS|^2 against R_CYL^2 * |AXIS|^2.
 * The cross product c = u x AXIS is itself maintained incrementally:
 * a unit move e_i changes c by the constant vector e_i x AXIS, and
 * |c|^2 is updated by a running increment G[m] whose own update is a
 * constant from the K table — the exact generalisation of sum-of-odds.
 */

#ifndef SPIRAL_H_
#define SPIRAL_H_

/* --- MSVC portability --- */
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
#define L 221

#define INF_R2    0xFFFFFFFFu
#define MID       (L / 2)
#define R_MAX     (L / 2)
#define RADIUS    (L / 2 - 2)
#define RADIUS_SQ ((long long)RADIUS * RADIUS)

/* --- Pulsating sweep parameters (L-invariant) --- */
#define PULSE_TOLERANCE  1
#define PULSE_STEP       (((L) + 15) / 30)

/* --- Cylindrical helix parameters ---
 * Compile-time integer approximation: 93/512 ~= 0.18164 of r_max. */
#define R_CYL       ((R_MAX * 93 + 256) >> 9)
#define R_CYL_SQ    (R_CYL * R_CYL)

/* Default rotation axis: (0,0,+L/2), matching the original +z helix.
 * Override at runtime: spiral ax ay az  (any direction; the vector is
 * rescaled so that |AXIS| = L/2). */
#define AXIS_X      0
#define AXIS_Y      0
#define AXIS_Z      R_MAX

/* Discrete (Manhattan) circumference of a radius-R_CYL circle:
 * axial steps cover 2*pi*R Euclidean distance in about 8R steps. */
#define CYL_CIRC    (R_CYL << 3)

/* --- Display --- */
#define WINDOW_W  (L + 500 + 80)
#define WINDOW_H  (L + 280)

/* Yellow ring visual tolerance (thin) */
#define YELLOW_VIS_TOL  ((RADIUS / 35) > 1 ? (RADIUS / 35) : 1)

/* =================================================================
 * Cell struct — minimal for spiral CA
 * ================================================================= */
typedef struct {
    int r;                /* integer radius from center */
    unsigned int r2;      /* Euclidean distance-squared (INF_R2 = unvisited) */
    unsigned char active; /* 1 if on pulsating shell */
    unsigned char spin;   /* 1 if on the spiral arm */
} Cell;

/* --- Spiral point storage (for walker + rendering) --- */
typedef struct { int x, y, z; } SpiralPt;

#define MAX_SPIRAL_PTS (R_MAX * 12)
extern SpiralPt spiral_pts[];
extern int spiral_n;
extern int spiral_done;

/* --- Grid pointers (heap-allocated) --- */
extern Cell (*grid)[L][L];
extern Cell (*grid_next)[L][L];
extern int tick;
extern int axis_x, axis_y, axis_z;   /* rotation axis, |AXIS| = L/2 */

/* =================================================================
 * Integer square root (bit-by-bit, no multiplication)
 * ================================================================= */
SINLINE int isqrt(int n) {
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
 * (uses multiplication — host-only, not part of the CA FSM)
 * ================================================================= */
SINLINE unsigned int pulse_from_time(unsigned int t) {
    const unsigned int min_r2 = 0;
    const unsigned int max_r2 = (unsigned int)((unsigned int)R_MAX * R_MAX * 92 / 100);
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

/* --- Host-side functions (spiral.c) --- */
void init(void);
void step_all(void);
void pulse_step(void);
void spiral_set_axis(int ax, int ay, int az);
void spiral_init(void);
void spiral_step(void);

#ifndef NO_SDL
void render_frame(SDL_Renderer *ren);
#endif

#endif /* SPIRAL_H_ */
