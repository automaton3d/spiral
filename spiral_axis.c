/*
 * spiral_axis.c — eccentric helical path from centre to an arbitrary point A,
 *                 using only integer arithmetic, cardinal (6-neighbour) steps,
 *                 and no floats.
 *
 * Geometry:
 *   C = (MID,MID,MID) is the grid centre.
 *   A = (AX,AY,AZ) is the endpoint, with |A-C| ~ L/2.
 *   D = A-C is the axis direction (and total displacement).
 *   The helix is an eccentric cylinder around an axis parallel to D, offset from C
 *   by a perpendicular vector OFF0 of length R.  One full rotation brings the
 *   walker from C back to the same side of the cylinder, ending at A.
 *
 * The offset vector lives in the plane perpendicular to D.  We build an integer
 * orthogonal basis (U,V) of that plane and walk a digital circle in the
 * anisotropic metric  a*u^2 + b*v^2 = R^2, where a=|U|^2, b=|V|^2.  Consecutive
 * ideal points P_i = (C-OFF0) + axis_DDA(i) + u_i*U + v_i*V are connected by a
 * 3D Bresenham line, producing only cardinal 6-neighbour steps.
 */

#include <stdio.h>
#include <stdlib.h>
#include <stdint.h>
#include <string.h>

#define L 221
#define MID (L/2)
#define RMAX (L/2 - 2)
#define MAX_PATH (L*L*4)

typedef struct { int x, y, z; } Vec;
typedef struct { int u, v; } UV;

static inline int absi(int a) { return a < 0 ? -a : a; }

static int gcd(int a, int b) {
    a = absi(a); b = absi(b);
    while (b) { int t = a % b; a = b; b = t; }
    return a;
}

static int64_t dot(Vec a, Vec b) {
    return (int64_t)a.x*b.x + (int64_t)a.y*b.y + (int64_t)a.z*b.z;
}
static Vec cross(Vec a, Vec b) {
    Vec r;
    r.x = a.y*b.z - a.z*b.y;
    r.y = a.z*b.x - a.x*b.z;
    r.z = a.x*b.y - a.y*b.x;
    return r;
}
static int64_t norm2(Vec v) {
    return dot(v,v);
}
static Vec add(Vec a, Vec b) {
    Vec r = {a.x+b.x, a.y+b.y, a.z+b.z};
    return r;
}
static Vec sub(Vec a, Vec b) {
    Vec r = {a.x-b.x, a.y-b.y, a.z-b.z};
    return r;
}
static Vec mulk(Vec a, int k) {
    Vec r = {a.x*k, a.y*k, a.z*k};
    return r;
}

static Vec reduce(Vec v) {
    int g = gcd(gcd(v.x, v.y), v.z);
    if (g == 0) return v;
    if (g < 0) g = -g;
    Vec r = {v.x/g, v.y/g, v.z/g};
    return r;
}

static int isqrt64(int64_t n) {
    if (n <= 0) return 0;
    int64_t result = 0;
    int64_t bit = 1LL << 62;
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
    return (int)result;
}

/* Find an orthogonal integer basis (U,V) of the plane perpendicular to D.
 * U is the shortest of the three obvious primitive perpendicular vectors;
 * V = reduce(U x D) is automatically orthogonal to both U and D. */
static void perp_basis(Vec D, Vec *U, Vec *V) {
    Vec c[3] = {
        {0, D.z, -D.y},
        {D.z, 0, -D.x},
        {-D.y, D.x, 0}
    };
    Vec cands[3];
    int nc = 0;
    for (int i = 0; i < 3; i++) {
        if (norm2(c[i]) > 0) cands[nc++] = reduce(c[i]);
    }

    /* choose the shortest nonzero candidate as U */
    Vec best = cands[0];
    for (int i = 1; i < nc; i++)
        if (norm2(cands[i]) > 0 && norm2(cands[i]) < norm2(best))
            best = cands[i];
    if (norm2(best) == 0) best = (Vec){0,1,0}; /* should not happen */
    *U = best;

    *V = reduce(cross(*U, D));
    if (norm2(*V) == 0) {
        *V = reduce(cross(D, *U));
    }
    if (norm2(*V) == 0) {
        /* fallback: pick any independent candidate */
        for (int i = 0; i < nc; i++)
            if (norm2(cross(*U, cands[i])) > 0) { *V = cands[i]; break; }
    }
}

/* Bresenham 3D line, 6-neighbour, from a to b.  Appends each new cell to path. */
typedef struct { Vec *pts; int n; int cap; } Path;

static void push(Path *p, Vec v) {
    if (p->n >= p->cap) return;
    p->pts[p->n++] = v;
}

static int64_t dist2(Vec a, Vec b) {
    int64_t dx = a.x - b.x, dy = a.y - b.y, dz = a.z - b.z;
    return dx*dx + dy*dy + dz*dz;
}

static void line3d(Path *p, Vec a, Vec b) {
    Vec v = a;
    push(p, v);
    if (v.x == b.x && v.y == b.y && v.z == b.z) return;
    const int d[6][3] = {{1,0,0},{-1,0,0},{0,1,0},{0,-1,0},{0,0,1},{0,0,-1}};
    while (1) {
        int best = -1;
        int64_t best_d = (int64_t)1 << 62;
        for (int k = 0; k < 6; k++) {
            Vec n = {v.x + d[k][0], v.y + d[k][1], v.z + d[k][2]};
            int64_t dd = dist2(n, b);
            if (dd < best_d) { best_d = dd; best = k; }
        }
        if (best < 0) break;
        v.x += d[best][0]; v.y += d[best][1]; v.z += d[best][2];
        push(p, v);
        if (v.x == b.x && v.y == b.y && v.z == b.z) break;
    }
}

/* Digital circle in the (u,v) anisotropic metric: a*u^2 + b*v^2 = R^2.
 * Returns the number of offset points and fills uv[0..n-1] in CCW order.
 * It walks until it returns to the start (or max_steps reached). */
static int circle_uv(int64_t a, int64_t b, int64_t R2, UV *uv, int max_steps) {
    /* start at a point on the ellipse with u>0, v=0 (if possible) */
    int u0 = 0;
    while ((int64_t)(u0+1)*(u0+1)*a <= R2) u0++;
    int v0 = 0;
    int u = u0, v = v0;
    int n = 0;
    uv[n++] = (UV){u, v};
    for (int step = 0; step < max_steps && n < max_steps; step++) {
        int64_t best_err = (int64_t)1 << 62;
        int best_du = 0, best_dv = 0;
        int64_t best_tan = -(1LL << 62);
        /* tangent vector for the ellipse in (u,v) space is (-b*v, a*u) */
        int64_t tu = -b * v;
        int64_t tv =  a * u;
        int cand[4][2] = {{1,0},{-1,0},{0,1},{0,-1}};
        for (int k=0;k<4;k++) {
            int du = cand[k][0], dv = cand[k][1];
            int nu = u + du, nv = v + dv;
            int64_t err = a*(int64_t)nu*nu + b*(int64_t)nv*nv - R2;
            if (err < 0) err = -err;
            int64_t tang = (int64_t)tu*du + (int64_t)tv*dv;
            if (tang <= 0) continue; /* CCW only */
            if (err < best_err || (err == best_err && tang > best_tan)) {
                best_err = err; best_tan = tang; best_du = du; best_dv = dv;
            }
        }
        if (best_du == 0 && best_dv == 0) {
            /* no CCW candidate; try any direction that reduces error */
            for (int k=0;k<4;k++) {
                int du = cand[k][0], dv = cand[k][1];
                int nu = u+du, nv = v+dv;
                int64_t err = a*(int64_t)nu*nu + b*(int64_t)nv*nv - R2;
                if (err < 0) err = -err;
                if (err < best_err) {
                    best_err = err; best_du = du; best_dv = dv;
                }
            }
        }
        if (best_du == 0 && best_dv == 0) break;
        u += best_du; v += best_dv;
        uv[n++] = (UV){u, v};
        if (u == u0 && v == v0) { n--; break; } /* closed */
    }
    return n;
}

int main(int argc, char **argv) {
    Vec A;
    if (argc >= 4) {
        A.x = atoi(argv[1]);
        A.y = atoi(argv[2]);
        A.z = atoi(argv[3]);
    } else {
        A = (Vec){MID, MID, MID + RMAX}; /* default: +z, like old program */
    }
    Vec C = {MID, MID, MID};
    Vec D = sub(A, C);
    int64_t D2 = norm2(D);
    if (D2 == 0) {
        printf("A cannot be the centre.\n");
        return 1;
    }

    Vec U, V;
    perp_basis(D, &U, &V);
    printf("D=(%d,%d,%d) |D|^2=%lld\n", D.x, D.y, D.z, (long long)D2);
    printf("U=(%d,%d,%d) |U|^2=%lld\n", U.x, U.y, U.z, (long long)norm2(U));
    printf("V=(%d,%d,%d) |V|^2=%lld\n", V.x, V.y, V.z, (long long)norm2(V));
    printf("U.V=%lld\n", (long long)dot(U,V));

    int64_t a = norm2(U);
    int64_t b = norm2(V);
    /* choose helix radius.  It must be larger than the longest basis vector
     * of the perpendicular plane in order to resolve a digital circle.
     * We use twice the longer of |U|,|V|, clamped to a sensible window. */
    int64_t max_basis2 = a > b ? a : b;
    int max_basis = isqrt64(max_basis2);
    if (max_basis < 1) max_basis = 1;
    int R = max_basis * 2;
    if (R < 6) R = 6;
    if (R > MID/2) R = MID/2;               /* keep the half-turn inside the grid */
    int64_t R2 = (int64_t)R * R;
    printf("Helix radius R=%d, a=%lld b=%lld R2=%lld\n", R, (long long)a, (long long)b, (long long)R2);

    UV uv[4096];
    int nuv = circle_uv(a, b, R2, uv, 4096);
    if (nuv <= 0) {
        printf("Could not generate circle.\n");
        return 1;
    }
    printf("Circle has %d offset steps.\n", nuv);

    /* start offset: choose the first circle point */
    Vec off0 = add(mulk(U, uv[0].u), mulk(V, uv[0].v));
    Vec axis0 = sub(C, off0);   /* so that axis0 + off0 = C */

    Path p;
    p.pts = (Vec*)malloc(MAX_PATH * sizeof(Vec));
    if (!p.pts) { printf("oom\n"); return 1; }
    p.n = 0; p.cap = MAX_PATH;

    /* DDA state for the axis advance over nuv steps */
    Vec axis = axis0;
    Vec axis_err = {0,0,0};
    Vec prevP = add(axis, off0);

    for (int i = 0; i <= nuv; i++) {
        /* offset for this angle */
        UV cur = uv[i % nuv];
        Vec off = add(mulk(U, cur.u), mulk(V, cur.v));
        /* axis position: axis0 + (i * D) / nuv via DDA */
        if (i > 0) {
            axis_err = add(axis_err, D);
            /* move one cell in each coordinate when accumulator crosses nuv */
            while (axis_err.x >= nuv) { axis.x++; axis_err.x -= nuv; }
            while (axis_err.x <= -nuv) { axis.x--; axis_err.x += nuv; }
            while (axis_err.y >= nuv) { axis.y++; axis_err.y -= nuv; }
            while (axis_err.y <= -nuv) { axis.y--; axis_err.y += nuv; }
            while (axis_err.z >= nuv) { axis.z++; axis_err.z -= nuv; }
            while (axis_err.z <= -nuv) { axis.z--; axis_err.z += nuv; }
        }
        Vec P = add(axis, off);
        line3d(&p, prevP, P);
        prevP = P;
    }

    printf("Path points: %d\n", p.n);
    printf("Start: (%d,%d,%d)  End: (%d,%d,%d)\n",
           p.pts[0].x, p.pts[0].y, p.pts[0].z,
           p.pts[p.n-1].x, p.pts[p.n-1].y, p.pts[p.n-1].z);
    printf("Expected end: (%d,%d,%d)\n", A.x, A.y, A.z);

    /* basic sanity: distance to target */
    Vec end = p.pts[p.n-1];
    int64_t d2 = (int64_t)(end.x-A.x)*(end.x-A.x) +
                 (int64_t)(end.y-A.y)*(end.y-A.y) +
                 (int64_t)(end.z-A.z)*(end.z-A.z);
    printf("Distance^2 to A: %lld\n", (long long)d2);

    /* render a small ASCII slice at z=MID to see the spiral projection */
    char grid[70][70];
    memset(grid, ' ', sizeof(grid));
    int ox = 35, oy = 35;
    for (int i = 0; i < p.n; i++) {
        int x = p.pts[i].x - MID + ox;
        int y = p.pts[i].y - MID + oy;
        if (x >= 0 && x < 70 && y >= 0 && y < 70)
            grid[y][x] = '.';
    }
    grid[oy][ox] = 'C';
    grid[A.y - MID + oy][A.x - MID + ox] = 'A';
    for (int y = 0; y < 35; y++) {
        for (int x = 0; x < 70; x++) putchar(grid[y][x]);
        putchar('\n');
    }

    free(p.pts);
    return 0;
}
