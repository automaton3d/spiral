/*
 * ac_rotation.c — Integer-only 3D rotation around the z-axis
 *
 * Method: Paeth 3-shear decomposition with Bresenham accumulators.
 *
 *   R(θ) = Shear_x(-tan θ/2) · Shear_y(sin θ) · Shear_x(-tan θ/2)
 *
 * The angle θ is defined by a pair (sin_p, sin_q) where sin(θ) = sin_p/sin_q.
 * All arithmetic is integer: addition, subtraction, comparison only.
 * Bresenham propagation replaces multiplication.
 *
 * Bit count is preserved exactly (bijective shears).
 *
 * Build:
 *   gcc -O2 -o ac_rotation ac_rotation.c -lm
 *   cl /Ox /Fe:ac_rotation.exe ac_rotation.c
 */

#include <stdio.h>
#include <stdlib.h>
#include <string.h>

#ifndef L
#define L 128
#endif

#define MID  (L / 2)
#define TOTAL ((size_t)L * L * L)

typedef unsigned char Cell;

#define IDX(x, y, z) ((size_t)(z) * L * L + (size_t)(y) * L + (size_t)(x))

/* ------------------------------------------------------------------ */
/*  Integer square root (Bresenham-style, no multiply in the loop)    */
/* ------------------------------------------------------------------ */
static int isqrt_int(int n)
{
    if (n <= 0) return 0;
    int r = 0, bit = 1 << 30;
    while (bit > n) bit >>= 2;
    while (bit) {
        if (n >= r + bit) { n -= r + bit; r = (r >> 1) + bit; }
        else              { r >>= 1; }
        bit >>= 2;
    }
    return r;
}

/* ------------------------------------------------------------------ */
/*  Compute cumulative shift per row/col using Bresenham              */
/*  shifts[center] = 0                                                */
/*  shifts[center ± d] = ±floor(d * p / q)                           */
/*  Only uses: add, sub, compare — no multiply.                       */
/* ------------------------------------------------------------------ */
static void compute_shifts(int *shifts, int p, int q, int center, int sign)
{
    int acc = 0, total = 0;
    shifts[center] = 0;
    for (int d = 1; d < L; d++) {
        acc += p;
        while (acc >= q) { acc -= q; total++; }
        if (center + d < L) shifts[center + d] = sign * total;
        if (center - d >= 0) shifts[center - d] = -sign * total;
    }
}

/* ------------------------------------------------------------------ */
/*  Shear content in x-direction: row y shifts by shifts[y]           */
/* ------------------------------------------------------------------ */
static void shear_x(Cell *grid, const int *shifts)
{
    Cell *tmp = (Cell *)calloc(TOTAL, sizeof(Cell));
    for (int z = 0; z < L; z++) {
        for (int y = 0; y < L; y++) {
            int s = shifts[y];
            for (int x = 0; x < L; x++) {
                int sx = x - s;
                if (sx >= 0 && sx < L)
                    tmp[IDX(x, y, z)] = grid[IDX(sx, y, z)];
            }
        }
    }
    memcpy(grid, tmp, TOTAL * sizeof(Cell));
    free(tmp);
}

/* ------------------------------------------------------------------ */
/*  Shear content in y-direction: column x shifts by shifts[x]        */
/* ------------------------------------------------------------------ */
static void shear_y(Cell *grid, const int *shifts)
{
    Cell *tmp = (Cell *)calloc(TOTAL, sizeof(Cell));
    for (int z = 0; z < L; z++) {
        for (int x = 0; x < L; x++) {
            int s = shifts[x];
            for (int y = 0; y < L; y++) {
                int sy = y - s;
                if (sy >= 0 && sy < L)
                    tmp[IDX(x, y, z)] = grid[IDX(x, sy, z)];
            }
        }
    }
    memcpy(grid, tmp, TOTAL * sizeof(Cell));
    free(tmp);
}

/* ------------------------------------------------------------------ */
/*  Paeth 3-shear rotation around z                                   */
/*                                                                    */
/*  Input: (sin_p, sin_q) where sin(θ) = sin_p / sin_q               */
/*                                                                    */
/*  Derived:                                                          */
/*    cos(θ) ≈ isqrt(sin_q² - sin_p²) / sin_q                       */
/*    tan(θ/2) = sin_p / (sin_q + cos_num)                           */
/*                                                                    */
/*  Shear 1,3: x by -tan(θ/2) · (y - cy)                            */
/*  Shear 2  : y by +sin(θ)   · (x - cx)                            */
/* ------------------------------------------------------------------ */
void rotate_z(Cell *grid, int sin_p, int sin_q)
{
    int center = MID;

    /* tan(θ/2) = sin_p / (sin_q + cos_num) */
    int cos_num = isqrt_int(sin_q * sin_q - sin_p * sin_p);
    int tan_p = sin_p;
    int tan_q = sin_q + cos_num;

    int *shifts = (int *)malloc(L * sizeof(int));

    printf("  Rotation: sin=%d/%d  tan_half=%d/%d\n", sin_p, sin_q, tan_p, tan_q);

    /* Shear 1: x by -tan(θ/2) · dy */
    compute_shifts(shifts, tan_p, tan_q, center, -1);
    shear_x(grid, shifts);

    /* Shear 2: y by +sin(θ) · dx */
    compute_shifts(shifts, sin_p, sin_q, center, +1);
    shear_y(grid, shifts);

    /* Shear 3: x by -tan(θ/2) · dy  (same as shear 1) */
    compute_shifts(shifts, tan_p, tan_q, center, -1);
    shear_x(grid, shifts);

    free(shifts);
}

/* ------------------------------------------------------------------ */
/*  Init: centered cube L/2 × L/2 × L/2 with asymmetric flag         */
/* ------------------------------------------------------------------ */
static void init_cube(Cell *grid)
{
    int half = L / 4;
    int start = MID - half;
    int end = start + L / 2;

    for (int z = start; z < end; z++)
        for (int y = start; y < end; y++)
            for (int x = start; x < end; x++)
                grid[IDX(x, y, z)] = 1;

    /* asymmetric flag so rotation is visually obvious */
    for (int z = start; z < start + 2 && z < L; z++)
        for (int y = start; y < start + 2 && y < L; y++)
            for (int x = end; x < end + 3 && x < L; x++)
                grid[IDX(x, y, z)] = 1;
}

/* ------------------------------------------------------------------ */
/*  Count active cells                                                */
/* ------------------------------------------------------------------ */
static int count_bits(const Cell *grid)
{
    int n = 0;
    for (size_t i = 0; i < TOTAL; i++) n += (grid[i] & 1);
    return n;
}

/* ------------------------------------------------------------------ */
/*  Save to .dat for visualization                                    */
/* ------------------------------------------------------------------ */
static void save_dat(const Cell *grid, const char *filename)
{
    FILE *f = fopen(filename, "w");
    fprintf(f, "# L = %d\n# Bits: %d\nx y z\n", L, count_bits(grid));
    for (int z = 0; z < L; z++)
        for (int y = 0; y < L; y++)
            for (int x = 0; x < L; x++)
                if (grid[IDX(x, y, z)] & 1)
                    fprintf(f, "%d %d %d\n", x, y, z);
    fclose(f);
}

/* ------------------------------------------------------------------ */
/*  Main                                                              */
/* ------------------------------------------------------------------ */
int main(int argc, char **argv)
{
    int sin_p = 1, sin_q = 2;  /* default: sin(θ) = 1/2 → θ = 30° */

    if (argc >= 3) {
        sin_p = atoi(argv[1]);
        sin_q = atoi(argv[2]);
    }

    if (sin_p <= 0 || sin_q <= 0 || sin_p > sin_q) {
        printf("Usage: %s [sin_p sin_q]  where 0 < sin_p <= sin_q\n", argv[0]);
        return 1;
    }

    Cell *grid = (Cell *)calloc(TOTAL, sizeof(Cell));

    printf("L=%d  Grid=%zu cells\n", L, TOTAL);
    init_cube(grid);
    save_dat(grid, "before_rotation.dat");

    int bits_before = count_bits(grid);
    printf("Bits before: %d\n", bits_before);

    rotate_z(grid, sin_p, sin_q);

    int bits_after = count_bits(grid);
    printf("Bits after:  %d\n", bits_after);
    printf("Preserved:   %s (delta=%d)\n",
           bits_before == bits_after ? "YES" : "NO",
           bits_after - bits_before);

    save_dat(grid, "after_rotation.dat");
    printf("Output: before_rotation.dat, after_rotation.dat\n");

    free(grid);

#ifdef _WIN32
    printf("\nPressione Enter para sair...\n");
    getchar();
#endif
    return 0;
}
