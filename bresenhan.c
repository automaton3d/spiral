/* bresenhan.c - Espiral Esférica começando na origem */
#include <SDL3/SDL.h>
#include <stdio.h>
#include <stdlib.h>
#include <math.h>

#define R 1000
#define MAXPTS 12000
#define WINDOW_WIDTH  1400
#define WINDOW_HEIGHT 920
#define PI 3.141592653589793

typedef struct { int x, y, z; } Point;

Point points[MAXPTS];
int point_count = 0;

void generate_spherical_spiral(void) {
    point_count = 0;
    int steps = 3000;

    for (int i = 0; i <= steps; i++) {
        double t = (double)i / steps;
        
        // phi de  π/2 (equador, origem) até 0 (polo norte)
        double phi = (PI / 2.0) * (1.0 - t);     
        double theta = 2.0 * PI * t;              // exatamente 1 volta

        double radius_xy = R * sin(phi);

        int x = (int)(radius_xy * cos(theta));
        int y = (int)(radius_xy * sin(theta));
        int z = (int)(R * cos(phi));              // de 0 até +R

        points[point_count++] = (Point){x, y, z};
    }
}

void project_iso(float wx, float wy, float wz, float scale, float ox, float oy,
                 float *sx, float *sy) {
    *sx = (wx - wy) * 0.7071f * scale + ox;
    *sy = (wx + wy) * 0.4082f * scale - wz * 0.8165f * scale + oy;
}

int main(void) {
    SDL_Init(SDL_INIT_VIDEO);
    SDL_Window *window = SDL_CreateWindow("Espiral Esférica - Começa na Origem (1 volta)",
                                        WINDOW_WIDTH, WINDOW_HEIGHT, 0);
    SDL_Renderer *renderer = SDL_CreateRenderer(window, NULL);

    generate_spherical_spiral();
    printf("Gerados %d pontos - Começa em (0,0,0)\n", point_count);

    float scale = 0.29f;
    float offset_x = WINDOW_WIDTH * 0.5f;
    float offset_y = WINDOW_HEIGHT * 0.52f;

    int running = 1;
    SDL_Event event;

    while (running) {
        while (SDL_PollEvent(&event)) {
            if (event.type == SDL_EVENT_QUIT) running = 0;
            if (event.type == SDL_EVENT_KEY_DOWN) {
                if (event.key.key == SDLK_ESCAPE) running = 0;
                if (event.key.key == SDLK_LEFT)  scale *= 1.06f;
                if (event.key.key == SDLK_RIGHT) scale *= 0.94f;
                if (event.key.key == SDLK_UP)    offset_y -= 12;
                if (event.key.key == SDLK_DOWN)  offset_y += 12;
            }
        }

        SDL_SetRenderDrawColor(renderer, 5, 5, 25, 255);
        SDL_RenderClear(renderer);

        float z_scale = 1.0f;

        for (int i = 0; i < point_count - 1; i++) {
            float t = (float)i / point_count;
            Uint8 r = (Uint8)(40 + t * 200);
            Uint8 g = (Uint8)(100 + t * 150);
            Uint8 b = (Uint8)(255 - t * 100);

            float iso_x1, iso_y1, iso_x2, iso_y2;
            project_iso((float)points[i].x, (float)points[i].y, (float)points[i].z * z_scale,
                        scale, offset_x, offset_y, &iso_x1, &iso_y1);
            project_iso((float)points[i+1].x, (float)points[i+1].y, (float)points[i+1].z * z_scale,
                        scale, offset_x, offset_y, &iso_x2, &iso_y2);

            SDL_SetRenderDrawColor(renderer, r, g, b, 255);
            SDL_RenderLine(renderer, iso_x1, iso_y1, iso_x2, iso_y2);
        }

        // Eixos
        float axis_len = R * 1.45f;
        float ox, oy;
        project_iso(0,0,0, scale, offset_x, offset_y, &ox, &oy);

        float ex, ey;
        project_iso(axis_len, 0, 0, scale, offset_x, offset_y, &ex, &ey);
        SDL_SetRenderDrawColor(renderer, 255, 70, 70, 255);
        SDL_RenderLine(renderer, ox, oy, ex, ey);

        project_iso(0, axis_len, 0, scale, offset_x, offset_y, &ex, &ey);
        SDL_SetRenderDrawColor(renderer, 70, 255, 70, 255);
        SDL_RenderLine(renderer, ox, oy, ex, ey);

        project_iso(0, 0, axis_len, scale, offset_x, offset_y, &ex, &ey);
        SDL_SetRenderDrawColor(renderer, 100, 200, 255, 255);
        SDL_RenderLine(renderer, ox, oy, ex, ey);

        SDL_RenderPresent(renderer);
        SDL_Delay(16);
    }

    SDL_Quit();
    return 0;
}