/*
 * cordic.h — header para emergent spiral CA on pulsating wavefront
 */

#ifndef CORDIC_H_
#define CORDIC_H_

/* --- MSVC portability --- */
#if defined(_MSC_VER) && !defined(__cplusplus)
#define SINLINE static __inline
#else
#define SINLINE static inline
#endif

/* --- SDL (apenas para renderização gráfica) --- */
#ifndef NO_SDL
#define SDL_MAIN_HANDLED  /* Impede que a SDL altere o ponto de entrada main() */
#include <SDL3/SDL.h>
#endif

#include <stdint.h>

/* --- Parâmetros da grade --- */
#define L 7
#define MID (L/2)
#define CELL_SIZE 20

/* --- Estrutura da célula --- */
typedef struct {
    int x, y, z;   /* coordenadas absolutas */
    int cx, cy;    /* controle interno CORDIC */
    int r;         /* contador radial */
    int active;    /* bit de ativação */
} Cell;

/* --- Prototipagem das funções --- */
void inicializar(Cell grid[L][L][L]);
void evoluir(Cell grid[L][L][L]);
void project_iso(int x, int y, int z, int *px, int *py);

#endif /* CORDIC_H_ */
