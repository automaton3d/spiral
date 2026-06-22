/*
 * cordic.c — implementação principal da emergent spiral CA
 */

#include "cordic.h"
#include <stdio.h>

/* --- Inicialização da grade --- */
void inicializar(Cell grid[L][L][L]) {
    for (int i=0;i<L;i++)
    for (int j=0;j<L;j++)
    for (int k=0;k<L;k++) {
        grid[i][j][k].x = i;
        grid[i][j][k].y = j;
        grid[i][j][k].z = k;
        grid[i][j][k].cx = 0;
        grid[i][j][k].cy = 0;
        grid[i][j][k].r = 0;
        grid[i][j][k].active = 0;
    }
}

/* --- Evolução da grade --- */
void evoluir(Cell grid[L][L][L]) {
    for (int i=0;i<L;i++)
    for (int j=0;j<L;j++)
    for (int k=0;k<L;k++) {
        int dx = grid[i][j][k].x - MID;
        int dy = grid[i][j][k].y - MID;

        for (grid[i][j][k].r = 0; grid[i][j][k].r < L/2; grid[i][j][k].r++) {
            int cond1 = (grid[i][j][k].cx * grid[i][j][k].cx == dx * dx);
            int cond2 = (grid[i][j][k].cy * grid[i][j][k].cy == dy * dy);
            int cond3 = (grid[i][j][k].r * grid[i][j][k].r == dx * dx + dy * dy);

            if (cond1 && cond2 && cond3) {
                grid[i][j][k].active = 1;
            }

            grid[i][j][k].cx++;
            grid[i][j][k].cy++;
        }
    }
}

/* --- Projeção isométrica simples --- */
void project_iso(int x, int y, int z, int *px, int *py) {
    *px = (x - y) * CELL_SIZE + 300;
    *py = (x + y) / 2 * CELL_SIZE - z * CELL_SIZE + 300;
}

/* --- Função principal --- */
int main(void) {
    Cell current[L][L][L];
    inicializar(current);
    evoluir(current);

#ifndef NO_SDL

    if (!SDL_Init(SDL_INIT_VIDEO))
    {
        printf("Erro SDL_Init: %s\n", SDL_GetError());
        return 1;
    }
    
    SDL_Window *win = SDL_CreateWindow(
        "Helicoidal 3D",
        600, 600,
        SDL_WINDOW_RESIZABLE);
    if (!win) {
        printf("Erro SDL_CreateWindow: %s\n", SDL_GetError());
        SDL_Quit();
        return 1;
    }

    SDL_Renderer *ren = SDL_CreateRenderer(win, NULL);
    if (!ren) {
        printf("Erro SDL_CreateRenderer: %s\n", SDL_GetError());
        SDL_DestroyWindow(win);
        SDL_Quit();
        return 1;
    }

int running = 1;
while (running) {
    SDL_Event e;
    while (SDL_PollEvent(&e)) {
        if (e.type == SDL_EVENT_QUIT) running = 0;
    }

    evoluir(current);  // recalcula a cada frame

    SDL_SetRenderDrawColor(ren, 0, 0, 0, 255);
    SDL_RenderClear(ren);

    for (int i=0;i<L;i++)
    for (int j=0;j<L;j++)
    for (int k=0;k<L;k++) {
        if (current[i][j][k].active) {
            int px, py;
            project_iso(i, j, k, &px, &py);
            SDL_SetRenderDrawColor(ren, 255, 0, 0, 255);
            SDL_FRect rect = { (float)px, (float)py, 6.0f, 6.0f };
            SDL_RenderFillRect(ren, &rect);
        }
    }

    SDL_RenderPresent(ren);
    SDL_Delay(16);
}

    SDL_DestroyRenderer(ren);
    SDL_DestroyWindow(win);
    SDL_Quit();
#else
    printf("Renderização SDL desativada (NO_SDL definido).\n");
#endif

    return 0;
}
