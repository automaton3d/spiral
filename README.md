\# Espiral emergente em torno de um eixo arbitrário



Versão generalizada de `spiral.c` / `spiral.h`: o helicoide cilíndrico emerge

em torno de qualquer eixo `AXIS = (ax, ay, az)` com `|AXIS| = L/2`, e não mais

apenas em torno de `+z`.



\## Ideia central (sem multiplicação em tempo de execução)



Para um ponto `u = v - P` (com `P` perpendicular a `AXIS`, `|P| = R\_CYL`,

deslocamento da linha do cilindro em relação ao centro da rede):



```

dist²(u, eixo) = |u × AXIS|² / |AXIS|²

```



Logo o caminhante compara `E = |c|²` (com `c = u × AXIS`) contra

`TGT = R\_CYL² · |AXIS|²` — nenhuma divisão, nenhum float.



O truque que elimina a multiplicação é a atualização incremental em dois níveis

(generalização exata da "soma de ímpares" `2·dx + 1`):



| grandeza | atualização por passo `m` | custo |

|---|---|---|

| `c = u × AXIS` | `c += D\[m]`, com `D\[m] = e\_m × AXIS` constante | soma |

| `E = |c|²` | `E += G\[m]` | soma |

| `G\[k]` | `G\[k] += K\[m]\[k]`, `K\[m]\[k] = 2·D\_m·D\_k` constante | soma |

| `q2 = |v|²` | `q2 ± (2·v\_i ± 1)` | shift + soma |



`D\[m]`, `K\[m]\[k]`, `TGT`, `TOL` e `P` são calculados \*\*uma única vez\*\* em

`spiral\_set\_axis()` / `spiral\_init()` (tempo de inicialização, exatamente como

as constantes `2·dx+1` da versão antiga). O laço da CA usa só `+`, `-`, `<<`

e comparações, em vizinhança local de 6 células.



\## Escolha do passo



\- \*\*Bresenham de passo\*\*: `climb\_total = |ax|+|ay|+|az|` passos de subida

&#x20; distribuídos entre `orbit\_total` passos orbitais.

\- \*\*Passo orbital\*\*: entre os 6 vizinhos, mantém-se apenas os de sentido de

&#x20; rotação correto (projeção tangencial `t = AXIS × u = -c`, cujas componentes

&#x20; já são `-c\_x, -c\_y, -c\_z` — de graça); entre esses, escolhe-se o de maior

&#x20; avanço angular dentro da banda radial `|E' - TGT| ≤ TOL`; fora da banda,

&#x20; o de menor erro radial.

\- \*\*Passo de subida\*\*: DDA multidimensional ao longo de `AXIS`

&#x20; (`dda\[i] += |a\_i|`, maior vence, `dda\[i] -= |a|₁`), seguindo o eixo inclinado

&#x20; com moves cardinais.

\- \*\*Fim\*\*: quando `q2 ≥ RADIUS²`, isto é `r = 0 → L/2` (θ: 0 → π).



\## Compilar e rodar



```sh

\# com SDL3 (visualização)

gcc -O2 spiral.c -o spiral $(pkg-config --cflags --libs sdl3) -lm

./spiral 3 -2 5        # eixo arbitrário; o vetor é reescalado para |AXIS| = L/2



\# teste geométrico headless

gcc -DNO\_SDL -O2 spiral.c test\_axis.c -o test\_axis -lm \&\& ./test\_axis

```



\## Resultado do teste (L = 221, R\_CYL = 20, RADIUS = 108)



Para eixos `(0,0,1)`, `(1,1,1)`, `(3,-2,5)`, `(-7,4,1)`, `(1,2,0)`, `(2,3,6)`

a distância ao eixo do cilindro fica em `18.1 … 22.0` (alvo 20, erro ≤ \~2

células, limite de quantização da rede), o raio final chega a `r ≈ 108 = L/2`

e o braço completa ao menos uma volta (θ ≥ 2π de fase orbital) em todos os casos.



