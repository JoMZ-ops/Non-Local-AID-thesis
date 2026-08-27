# Notas de lectura — Polonyi, arXiv:1701.04068v4

*The Abraham-Lorentz force and electrodynamics at the classical electron radius*
(J. Polonyi, Strasbourg / CNRS-IPHC, versión v4 fechada 6 jun 2019, 22 páginas).

Estas notas fijan la transcripción de las ecuaciones que implementa `nlaid/`.
**Verificar contra el PDF original antes de citar cualquiera en la tesis.**

---

## 1. Tesis central

La fuerza de Abraham-Lorentz no es un defecto de la electrodinámica clásica que
haya que parchear: es un **efecto de cutoff UV mal entendido**.

1. Al eliminar el campo EM (resolver Maxwell con la línea de mundo como fuente y
   reinsertar), los pares de patas externas se acoplan y **se forman loops**.
   Son de orden `O(hbar^0)` — clásicos. De ahí la divergencia UV en una teoría
   clásica (introducción, punto A, p. 3).
2. El loop diverge, requiere regulador, y **el regulador debe ser
   Lorentz-invariante**: el campo cercano divergente amplifica cualquier ruptura
   de simetría (punto F, p. 4).
3. La fuerza de A-L es **cutoff-independiente pero generada por el cutoff**: la
   estructura de una **anomalía**. El integrando diverge como `O(ell^-1)` pero el
   resultado es finito, por convergencia no uniforme. Polonyi lo compara
   explícitamente con la anomalía quiral (p. 8).

## 2. El regulador — la no-localidad

Se reemplaza la delta en la función de Green retardada, `delta(x^2) -> delta_B(x^2)`,
sujeta a tres condiciones (p. 6):

| # | Condición | Razón física |
|---|---|---|
| (i) | `int dz delta_B(z) = 1` | preservar el flujo del campo radiado |
| (ii) | `delta_B(0) = 0` | separar puntos singulares en `Theta(x^0) delta_B(x^2)` |
| (iii) | `delta_B(z) = 0` para `z < 0` | **suprimir la interacción superlumínica** |

La condición (iii) es el mecanismo entero: una trayectoria *runaway* tendría que
superar `c` para escapar, y ahí el regulador apaga la autointeracción.

Dos realizaciones:

- **Ec. (4), desplazada:** `delta_B(x^2) = delta(x^2 - ell^2)`
- **Ec. (5), suavizada:** `delta_B(x^2) = Theta(x^2) x^2 exp(-sqrt(x^2)/ell) / (12 ell^4)`

La separación no local es el **intervalo invariante de Minkowski**:

> "the retarded source point, x', is found by the condition `ell^2 = (x - x')^2`"
> (p. 9, bajo la ec. 16)

> "smears the electromagnetic field over the invariant distance `ds^2 = ell^2`"
> (conclusiones, p. 19)

**No** es `x^2 - x'^2`. La distinción cambia el núcleo de la autointeracción.

## 3. Ecuaciones implementadas

### Ec. (6) — forma general

    xddot = 4 r_0B int_{-inf}^{s} ds' delta_B'((x-x')^2)
            { (x-x')(xdot.xdot') - [xdot.(x-x')] xdot' }

### Ec. (14) — susceptibilidad (bloque 1)

Con `F^r_omega = 1/[(omega + i eps)^2 chi^r_omega]`:

    chi^r_omega = 1 + r0 [ (2/3) i omega
                  - (2/omega^2) int_{-inf}^{0} du delta_B(u^2) N(omega u)/u^2 ]

    N(phi) = (1 + i phi - phi^2) e^{-i phi} - 1 + phi^2/2 - (2/3) i phi^3

Estable y causal si `chi^r` es analítica y sin ceros en `Im omega > 0` (p. 8).

### Ec. (16) — regulador desplazado, retardo finito (bloque 2)

    xddot = r0 (m/m_B) (1/A^2) [ (B+1)/A * V1 + V2 ]
    A  = xdot'.(x-x')
    B  = xddot'.(x'-x)
    V1 = (x-x')(xdot.xdot') - [xdot.(x-x')] xdot'
    V2 = (x-x')(xdot.xddot') + [xdot.(x'-x)] xddot'

con `x'` fijado por `ell^2 = (x-x')^2`.

### Ec. (17) — regulador suavizado, memoria infinita (bloque 3)

    xddot = (r0/3 ell^4)(m/m_B) int_{-inf}^{0} du
            (1 - sqrt(w)/(2 ell)) exp(-sqrt(w)/ell)
            { (x-x')(xdot.xdot') + [xdot.(x'-x)] xdot' }

con `w = (x-x')^2` y `x' = x(s+u)`.

### Ec. (18) — condición de renormalización

    chi^r_omega = 1 + (2/3) i r0 omega

## 4. Resultados derivados a mano (usados como tests)

Estas cantidades **no** están en el paper; se derivaron para validar el código
contra valores cerrados en lugar de contra otra corrida del propio código.

**Serie del numerador de la ec. (14).** Los cuatro primeros órdenes se cancelan
idénticamente con la sustracción, dejando

    N(phi) = sum_{n>=4} (-i phi)^n (n-1)^2 / n!

Se comprueba sumando `(n-1)^2 x^n/n! = (x^2 - x + 1) e^x` con `x = -i phi`, que
reproduce `(1 + i phi - phi^2) e^{-i phi}`. El orden dominante es
`N = (3/8) phi^4 - (2/15) i phi^5 + O(phi^6)`, **coincidiendo exactamente** con la
afirmación del paper (p. 8) de que la función racional es `(3/8) omega^4 u^2 (1 + O(omega u))`.
Es la verificación de que la transcripción de la ec. (14) es correcta.

**Momentos del regulador suavizado, ec. (5):**

| Cantidad | Valor |
|---|---|
| `int dz delta_B(z)` | `1` (confirma la condición (i)) |
| `int_0^inf dz z^{-1/2} delta_B(z)` | `1/(3 ell)` |
| `delta_m / m` vía ec. (10) | `r0 / (6 ell)` |
| `I2 = int_{-inf}^0 du delta_B(u^2) u^2` | `2 ell` |
| `chi^r_omega` a orden `omega^2` | `1 + r0[(2/3) i omega - (3/2) omega^2 ell]` |

Para el regulador desplazado, ec. (4): `I2 = ell/2`, `delta_m/m = r0/(2 ell)`, y
`chi^r_omega` tiene forma cerrada vía `int_{-inf}^0 du delta_B(u^2) f(u) = f(-ell)/(2 ell)`.

**Errata de extracción.** Extractores automáticos de texto leen la ec. (10) como
`int dz sqrt(z) delta_B(z)`, que es dimensionalmente imposible. La imagen de la
p. 7 muestra `int dz/sqrt(z) delta_B(z)`. Verificado contra el PDF.

**Conservación exacta de la normalización.** Contrayendo los lados derechos de
las ecs. (6), (16) y (17) con `xdot_mu` se obtiene **cero idénticamente**. Por
ejemplo para la (6):

    xdot.{(x-x')(xdot.xdot') - [xdot.(x-x')] xdot'}
        = [xdot.(x-x')](xdot.xdot') - [xdot.(x-x')](xdot'.xdot) = 0

Consecuencia: `xdot^2 = 1` no es una restricción a imponer sino una consecuencia
exacta, y su deriva numérica mide directamente el error de integración.
`WorldLine.norm_drift` la reporta.

**Semiplano de convergencia de la ec. (14).** Sustituyendo `u = -v`, el integrando
decae como `v^2 exp[-v(1/ell + Im omega)]`, de modo que la representación converge
para `Im omega > -1/ell` — todo el semiplano superior. La condición `Im omega > 0`
del paper (p. 8) es suficiente, no el límite real.

## 5. Naturaleza matemática de cada ecuación

| Ec. | Tipo | Consecuencia numérica |
|---|---|---|
| (14) | integral definida en `omega` complejo | cuadratura + búsqueda de ceros |
| (16) | DDE **neutra** con retardo dependiente del estado | root-finding + predictor-corrector por paso |
| (17) | Volterra no lineal de **2da especie** (reducida), núcleo dependiente de la solución | marcha explícita hacia adelante |

La (17) **no** es una ecuación integral de segunda especie en sentido estricto —
lo que aparece fuera de la integral es `xddot`, una derivada de la incógnita, así
que es integro-diferencial. Pero tomando `u = xdot` e integrando una vez, y
colapsando la integral doble por Fubini, `(x, u)` sí satisface un sistema no
lineal de Volterra de segunda especie. Eso es lo que autoriza la marcha explícita.

Tres salvedades:

1. **El núcleo depende de la solución** (`sqrt((x-x')^2)` en la exponencial). No
   hay núcleo resolvente ni serie de Neumann; los teoremas aplicables son de la
   familia no lineal, no del caso lineal clásico.
2. **El límite inferior es `-inf`.** La prehistoria es dato conocido y contribuye
   como forzante, lo que vuelve el problema bien puesto como **problema de valor
   inicial**. Es la razón de que Polonyi insista en el formalismo CTP (punto E, p. 4).
3. **El núcleo es regular**, no débilmente singular: cuando `s' -> s`, la
   exponencial tiende a 1 y el corchete se anula linealmente. No hay singularidad
   tipo Abel, así que la cuadratura compuesta converge a orden pleno.

## 6. Figuras del paper

- **Fig. 1** (p. 10): `|xddot| r0` vs `s/r0`, suavizado, `r0/ell = 3`.
  (a) `m/m_B = 1.95, 1.98, 2.0` — relajación exponencial monótona.
  (b) `m/m_B = -3.8, -3.91, -4.1` — oscilación con envolvente exponencial.
- **Fig. 2** (p. 11): diagrama de fases en `(r0/ell, m/m_B)`. Región sombreada
  estable; líneas sólidas = condición de renormalización (18); **línea punteada =
  teoría linealizada fijada por el contraterm (10)**.
- **Fig. 3** (p. 12): calidad del ajuste a la condición de renormalización.

## 7. Lo que queda fuera de `nlaid/`

La sección III (QED: acción efectiva CTP, regularización del propagador completo,
acción de línea de mundo, ecs. 19–40) es el andamiaje formal que justifica la
sección II. No se implementa: el objetivo del código es la dinámica clásica
efectiva de las ecs. (14), (16) y (17).
