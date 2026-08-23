# Espiral emergente em torno de um eixo arbitrário

Implementação consolidada do automato celular que gera uma hélice cilíndrica
em torno de um eixo inteiro arbitrário. A direção de rotação não é imposta:
é **eleita pelo próprio CA** a partir do campo `w` no instante em que a onda
pulsante atinge seu limite de expansão.

A descrição matemática detalhada está em `ca_essence.tex` / `ca_essence.pdf`.

## Arquivos principais

| arquivo | função |
|---|---|
| `spiral.h` | constantes, tipos (`Cell`, `SpiralPt`) e protótipos do core |
| `spiral.c` | core do CA: onda pulsante (`pulse_step`) e caminhante helicoidal (`spiral_step`) |
| `spiral_auto.c` | visualizador 3D interativo (SDL3), orquestração e eleição do eixo |
| `build_3d.bat` | atalho para compilar `spiral_auto.exe` com MSVC/nmake |
| `Makefile.nmake` | build Windows (MSVC) para `spiral_auto.exe` |
| `SDL3.dll` | biblioteca SDL3 pré-compilada para Windows |
| `ca_essence.tex` / `ca_essence.pdf` | manuscrito com a matemática do automato |

## Ideia central

1. **Respiração isotrópica:** uma onda triangular em `r²` aciona uma casca
   esférica (`active = 1`) que expande e contrai.
2. **Eleição da direção:** no pico de expansão, cada célula ativa propaga o
   maior *payload* `(score << 24) | code` entre seus 6 vizinhos faciais. A
   célula vencedora define o vetor `m` (eixo de rotação).
3. **Hélice:** o caminhante de Bresenham 3D traça a espiral em torno de `m`,
   reescalado para `|m| = L/2`, usando apenas `+`, `-`, `<<` e comparações
   (sem multiplicação, divisão, float ou trigonométricas no laço do CA).
4. **Broadcast (ant view):** cada ponto da espiral carimba seu `tick` numa
   segunda grade; uma difusão local leva esse valor a toda a rede, gerando
   o "pó" colorido por ciclo de respiração.

## Compilar

### Windows — MSVC

```cmd
nmake /f Makefile.nmake
rem ou
build_3d.bat
```

Requer SDL3 instalado em `E:\vcpkg\installed\x64-windows` (veja `SDL_INC` e
`SDL_LIB` no `Makefile.nmake` e ajuste se necessário).

### Linux / macOS — GCC

O visualizador precisa do SDL3. Com o SDL3 e `pkg-config` disponíveis:

```sh
gcc -O2 spiral_auto.c spiral.c -o spiral_auto $(pkg-config --cflags --libs sdl3) -lm
./spiral_auto
```

Para compilar apenas o core sem SDL (objeto de biblioteca):

```sh
gcc -DNO_SDL -O2 -c spiral.c -o spiral_nosdl.o
```

## Executar

```sh
./spiral_auto
```

Janela padrão: 1024×768, redimensionável.

## Controles

| tecla / mouse | ação |
|---|---|
| **LMB** + arrastar | orbitar câmera (perspectiva) |
| **RMB** + arrastar | pan |
| **wheel** | zoom |
| **1 … 5** | passos de simulação por frame (1, 2, 4, 8, 16) |
| **Espaço** | pausar / continuar |
| **A** | ligar/desligar rotação automática |
| **W** | incrementar `wseed` e re-eleger o eixo imediatamente |
| **D** | ligar/desligar *broadcast dust* |
| **S** | ligar/desligar esfera de onda analítica |
| **C** | ligar/desligar amostra das células ativas |
| **B** | ligar/desligar braço da espiral (voxels) |
| **G** | ligar/desligar eixos e esfera delimitadora |
| **X / Y / Z** | vista ortográfica ao longo do eixo correspondente |
| **I** | vista isométrica ortográfica |
| **R** | resetar câmera |
| **ESC** | sair |

## Teste geométrico headless

O arquivo `test_axis.c` (presente na branch `devin/arbitrary-axis`) pode ser
usado para validar o core sem GUI:

```sh
git show origin/devin/arbitrary-axis:test_axis.c > test_axis.c
gcc -DNO_SDL -O2 spiral.c test_axis.c -o test_axis -lm
./test_axis
```

A saída reporta, para vários eixos, o raio do cilindro aproximado, o raio
final e o ângulo total percorrido.

## Estado atual

- Core de eixo arbitrário funcional (validado com `test_axis.c`).
- Visualizador 3D com eleição automática do eixo, broadcast dust e vistas
  ortográficas X/Y/Z/I implementados.
- Build Linux/SDL3 depende de SDL3 instalado; não há `Makefile` GCC no
  momento.
