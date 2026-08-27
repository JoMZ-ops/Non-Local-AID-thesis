# Discrepancias abiertas con el paper

Registro de puntos donde el código **no** reproduce lo que afirma Polonyi
(arXiv:1701.04068v4). No están resueltos. Se documentan aquí en lugar de
ajustar el código hasta que "salga" la figura del paper.

---

## D1. Carácter de la solución para `m_B > 0` (Fig. 1a) — **RESUELTA**

**El paper dice** (p. 10): *"The acceleration changes in a monotonous,
exponential manner after some transient period, depending on the initial
conditions, if `m_B > 0` as shown in Fig. 1 (a)"*, con `r0/ell = 3` y
`m/m_B = 1.95, 1.98, 2.0`.

**El código da**: relajación exponencial, sí — con tasas `-0.741, -0.726, -0.717`
respectivamente, variando suavemente con `m/m_B` — pero **oscilatoria**, no
monótona: `|xddot|` presenta mínimos profundos y periódicos en escala log, es
decir la componente espacial con signo cambia de signo repetidamente.

**Resolución.** No es una discrepancia. El bloque 1 (análisis espectral
linealizado, ec. 14) localiza el modo dominante a `r0/ell = 3` en

    omega = +-2.5736 - 0.71793 i

es decir un **par complejo conjugado**: la relajación es necesariamente
oscilatoria, con tasa `Im omega = -0.71793` y período `2 pi / 2.5736 = 2.4414`.

Contraste con el bloque 3 (integración temporal no lineal, ec. 17) a
`m/m_B = 2.00`, por un camino completamente independiente:

| cantidad | bloque 1 (espectral) | bloque 3 (no lineal) |
|---|---|---|
| tasa de decaimiento | `-0.71793` | `-0.7167` |
| espaciado entre mínimos de `|xddot|` | `1.2207` (medio período) | `1.2207` |

Coinciden al 0.17% en la tasa y a cuatro cifras en la frecuencia. La palabra
"monotonous" del paper describe la **envolvente**, no la señal: su Fig. 1(a)
abarca `s/r0` de 50 a 300, donde caben ~100 oscilaciones de período 2.44 y solo
la envolvente es resoluble a esa escala de dibujo.

Vale la pena notar además que el contraterm de la ec. (10) da
`m/m_B = 1/(1 - r0/6ell) = 2.0000` exactamente a `r0/ell = 3` — el valor central
de la Fig. 1(a), que el paper bracketea con 1.95 y 1.98. Confirma
independientemente la fórmula del contraterm.

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
