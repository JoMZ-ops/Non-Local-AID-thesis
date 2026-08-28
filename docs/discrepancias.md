# Discrepancias y hallazgos sobre el paper

Registro de puntos donde el código **no** reproduce lo que afirma Polonyi
(arXiv:1701.04068v4). No están resueltos. Se documentan aquí en lugar de
ajustar el código hasta que "salga" la figura del paper.

---

## D1. Carácter de la solución para `m_B > 0` (Fig. 1a) — **ABIERTA**

**El paper dice** (p. 10): *"The acceleration changes in a monotonous,
exponential manner after some transient period, depending on the initial
conditions, if `m_B > 0` as shown in Fig. 1 (a)"*, con `r0/ell = 3` y
`m/m_B = 1.95, 1.98, 2.0`.

**El código da**: relajación exponencial, sí — con tasas `-0.741, -0.726, -0.717`
respectivamente, variando suavemente con `m/m_B` — pero **oscilatoria**, no
monótona: `|xddot|` presenta mínimos profundos y periódicos en escala log, es
decir la componente espacial con signo cambia de signo repetidamente.

**CORRECCIÓN (lectura de la Fig. 1 a resolución suficiente).** En una versión
previa de este documento se dio D1 por resuelta, argumentando que el
"monotonous" del paper describía la envolvente de una señal oscilatoria. **Eso
era incorrecto.** La Fig. 1(a) real, en escala lineal, muestra curvas suaves y
no oscilatorias sobre `s/r0` hasta **300**, con `a r0 ~ 0.1..0.5` todavía al
final del rango. Las tres curvas son: `1.95` decae, `1.98` sube a ~0.18 hacia
`s/r0 ~ 150` y decae lentamente, `2.0` explota fuera de escala. El borde de
estabilidad está entre 1.98 y 2.0.

El código da tasa de decaimiento `-0.72` y período de oscilación `2.44`. A
`s/r0 = 300` eso serían 123 oscilaciones y un decaimiento de `e^-216`. No es
el mismo régimen dinámico.

**Lo que sigue en pie.** El bloque 1 (espectral, ec. 14) y el bloque 3
(integración no lineal, ec. 17) coinciden entre sí al 0.17% en la tasa y a
cuatro cifras en la frecuencia, por caminos independientes:

| cantidad | bloque 1 | bloque 3 |
|---|---|---|
| tasa de decaimiento | `-0.71793` | `-0.7167` |
| espaciado entre mínimos | `1.2207` (medio período) | `1.2207` |

Que **ambos** discrepen del paper de la misma manera descarta un error de
integración temporal y apunta a un error común: una transcripción compartida
por las ecs. (14) y (17), o una interpretación distinta de los parámetros.

**Pista.** La Fig. 1(b) oscila con período `~1.1`, es decir `omega ~ 5.7..6.3`.
El bloque 1 encuentra ceros subdominantes en `omega = +-5.70 - 2.65i` y
`+-5.88 - 2.64i`. Esa frecuencia **sí está** en el espectro calculado, pero no
como modo dominante. Sugiere una diferencia de escala o de rama, no ecuaciones
distintas.

**Confirmado independientemente.** El contraterm de la ec. (10) da
`m/m_B = 1/(1 - r0/6 ell) = 2.0000` exactamente a `r0/ell = 3`, y el borde de
estabilidad de la Fig. 1(a) cae entre 1.98 y 2.0. La fórmula del contraterm
está bien; el desacuerdo está en la dinámica, no en la normalización de masa.

## D2. Estabilidad marginal en `m/m_B = -3.91` (Fig. 1b)

**El paper dice** (p. 10): con `m_B < 0` la aceleración es *"oscillatory with
exponentially exploding or decreasing envelope"*, y dibuja `-3.91` como la línea
sólida, es decir el valor que satisface la condición de renormalización (18) —
por tanto cercano al borde de estabilidad.

**El código da**: crecimiento violento. `|xddot|` pasa de `~0.5` a `1.1e6` en
`s = 12`, con tasa de envolvente `~ +1.45`. No es marginal en absoluto.

Ese tramo **sí está convergido**: `ds = 5e-3` y `ds = 2e-3` coinciden a cuatro
cifras (`1.0948e6` vs `1.0951e6`) con deriva `1.2e-4`. Más allá de `s ~ 13` la
deriva crece a `O(1)` y la integración deja de ser fiable — por eso las curvas
se truncan en la figura.

## Hipótesis a revisar (en orden de plausibilidad)

1. **La condición de renormalización no se está imponiendo.** El paper *no* elige
   `m/m_B` libremente: lo **fija** resolviendo la ec. (18),
   `chi^r_omega = 1 + (2/3) i r0 omega`, monitoreando la relajación a `s` grande
   (p. 10). Los valores de la Fig. 1 bracketean esa solución. Comparar corridas a
   `m/m_B` arbitrario contra la Fig. 1 puede ser sencillamente la comparación
   equivocada. **Esta es la vía más probable y la siguiente pieza a implementar.**

2. **La condición inicial difiere.** El paper prescribe una trayectoria y apaga la
   fuente; aquí se usa un pulso externo `C^infinito` (ver `worldline.py` para por
   qué: el corte abrupto degrada el orden del integrador). El propio paper avisa
   que *"The precise value of m/m_B at the stability edge is found to be slightly
   dependent on the initial, prescribed trajectory"* (p. 10). Eso explicaría un
   corrimiento del borde, pero probablemente no un cambio de carácter monótono
   a oscilatorio.

3. **Convención de signo o normalización en `m/m_B`.** Revisar si el `m/m_B` de la
   Fig. 1 es el mismo parámetro que multiplica `r0` en la ec. (17).

4. **Tensión interna en el propio paper.** El rango del eje vertical de la Fig. 2
   es aproximadamente `[-1.2, 0.7]`, que **no contiene** ninguno de los valores
   usados en la Fig. 1 (`1.95..2.0` y `-3.8..-4.1`). Puede ser lectura mía de una
   figura de baja resolución, o puede indicar que las dos figuras usan
   parametrizaciones distintas. **Verificar en el PDF original a alta resolución
   antes de sacar conclusiones.**

## Prueba cruzada pendiente

El contraste más fuerte disponible, aún no hecho: comparar la **tasa de
relajación medida** por el bloque 3 contra los **ceros de `chi^r_omega` en el
semiplano inferior** calculados por el bloque 1. Ambos bloques predicen la misma
cantidad por caminos completamente independientes (integración temporal no lineal
vs. análisis espectral linealizado). Si coinciden en el régimen de acoplamiento
débil, valida los dos a la vez; si no, localiza el error.


---

## D1 — hipótesis descartada: reescalado del cutoff

Se probó si la Fig. 1(a) correspondería a un `r0/ell` distinto del nuestro (por
ejemplo si el paper mezclara escalas entre reguladores, cf. F1 abajo). **No.**
El modo dominante de `chi^r` en el semiplano inferior, regulador suavizado:

| `r0/ell` | tasa `Im(omega)` | decaimiento a `s/r0 = 300` |
|---|---|---|
| 0.5 | `-0.3095` | `4.7e-41` |
| 1.0 | `-0.5338` | `2.8e-70` |
| 1.5 | `-0.7500` | `1.9e-98` |
| 2.0 | `-0.8346` | `1.8e-109` |
| 3.0 | `-0.7179` | `2.9e-94` |

Ningún cutoff da una relajación lo bastante lenta. La Fig. 1(a) muestra
`a r0 ~ 0.1..0.18` todavía en `s/r0 = 300`, lo que exige `|tasa| <~ 0.003`.
La discrepancia no es un reescalado de `ell`.

---

# Hallazgos sobre el paper

## F1. La Fig. 2(b) usa el contraterm del regulador DESPLAZADO

**Verificado por digitalización de la figura del PDF** (calibración: eje `y=0`
en la fila 446, 367 px por unidad; eje `x` con 59.36 px por unidad de `r0/ell`,
retícula de marcas menores unitarias confirmada en ambos paneles).

La línea punteada de **ambos** paneles sigue

    m/m_B = 1 / (1 - r0/(2 ell))

es decir el contraterm del regulador **desplazado**, ec. (4), con
`delta_m = e^2 / (2 c^2 ell)`.

Panel (a), regulador desplazado — como corresponde:

| `r0/ell` | medido | `1/(1-x/2)` | error |
|---|---|---|---|
| 7 | `-0.395` | `-0.400` | 0.005 |
| 8 | `-0.338` | `-0.333` | 0.004 |
| 10 | `-0.250` | `-0.250` | 0.000 |
| 12 | `-0.201` | `-0.200` | 0.001 |
| 15 | `-0.155` | `-0.154` | 0.001 |
| 18 | `-0.125` | `-0.125` | 0.000 |
| 20 | `-0.112` | `-0.111` | 0.001 |

Error medio `0.0017`. **La fórmula del contraterm es correcta.**

Panel (b), regulador **suavizado** — pero la punteada sigue la MISMA ley:

| `r0/ell` | medido | `1/(1-x/2)` (despl.) | `1/(1-x/6)` (suav.) |
|---|---|---|---|
| 11 | `-0.217` | `-0.222` | `-1.200` |
| 12 | `-0.201` | `-0.200` | `-1.000` |
| 13 | `-0.181` | `-0.182` | `-0.857` |
| 14 | `-0.166` | `-0.167` | `-0.750` |
| 15 | `-0.155` | `-0.154` | `-0.667` |
| 16 | `-0.142` | `-0.143` | `-0.600` |
| 17 | `-0.134` | `-0.133` | `-0.545` |

Error medio contra la ley del desplazado: `0.002`. Contra la del suavizado:
`0.66`. La superposición gráfica está en `figures/fig2b_verificacion_contraterm.png`
(naranja = ley del desplazado, verde = ley del suavizado; la naranja tapa
exactamente la punteada azul del paper).

**Aplicando la ec. (10) al regulador suavizado, ec. (5), el resultado correcto es**

    int_0^inf dz z^{-1/2} delta_B(z) = 1/(3 ell)   =>   delta_m = e^2/(6 c^2 ell)

verificado analíticamente (sustitución `z = t^2 ell^2`, integral gamma) y
numéricamente por cuadratura, junto con la normalización `int delta_B = 1` que
fija la definición del regulador. La curva correspondiente es la verde, que
cruza `m/m_B = -1` en `r0/ell = 12` en lugar de `r0/ell = 4`.

**Interpretación.** O bien es un desliz al graficar (la misma curva dibujada en
los dos paneles), o bien Polonyi define un `ell` efectivo distinto para el
regulador suavizado — por ejemplo `ell_ef = ell/3`, que igualaría ambos
`delta_m`. El paper no menciona tal reescalado. **Verificar antes de citar
cualquiera de las dos curvas en la tesis.**
