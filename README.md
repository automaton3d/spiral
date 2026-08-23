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
| `test_axis.c` | teste geométrico headless do helicoide para vários eixos |
| `Makefile` | build GCC/Linux do teste headless (e visualizador quando SDL3 estiver disponível) |
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

### Linux / macOS — GCC

Teste headless (não requer SDL):

```sh
make test
```

O target padrão também compila e roda o teste:

```sh
make
```

Visualizador 3D (requer SDL3 e `pkg-config`):

```sh
make spiral_auto
./spiral_auto
```

Caso `pkg-config` não encontre o SDL3, ajuste `SDL3_CFLAGS` e `SDL3_LIBS`
no `Makefile`.

### Windows — MSVC

```cmd
nmake /f Makefile.nmake
rem ou
build_3d.bat
```

Requer SDL3 instalado em `E:\vcpkg\installed\x64-windows` (veja `SDL_INC` e
`SDL_LIB` no `Makefile.nmake` e ajuste se necessário).

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

O `test_axis.c` verifica o helicoide sem GUI para vários eixos, reportando o
raio do cilindro aproximado, o raio final e o ângulo total percorrido.

```sh
make test
```

Exemplo de saída esperada (com `L = 201`, `R_CYL = 18`, `RADIUS = 98`):

```
axis(   0,   0,   1)->(   0,   0, 100) |A|=100.0  pts= 241  cyl r: 16.20..19.62 ...
axis(   1,   1,   1)->(  58,  58,  58) |A|=100.5  pts= 307  cyl r: 13.02..22.72 ...
...
```

Os valores devem mostrar raio final próximo de `RADIUS` e curva total
superior a `2π` rad (~6.28) para os eixos testados.

## Estado atual

- Core de eixo arbitrário funcional (validado por `test_axis.c`).
- Visualizador 3D com eleição automática do eixo, broadcast dust e vistas
  ortográficas X/Y/Z/I implementados.
- Teste headless e `Makefile` GCC/Linux restaurados.
- Build do visualizador em Linux depende do SDL3 instalado.
