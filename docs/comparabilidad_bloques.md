# Comparabilidad entre bloques: dónde vive `m/m_B` y dónde no

## El problema, en una línea

El plano de fases de la Fig. 2 tiene **dos** ejes, `(r₀/ℓ, m/m_B)`. Comparar un
resultado del bloque 1 contra uno de los bloques 2 o 3 fijándose solo en `r₀/ℓ`
es comparar dos puntos distintos del plano, y da conclusiones opuestas sobre
estabilidad con la misma ecuación.

Apareció al graficar la ec. (16) (`seccion_A/fig_ec16_retardo.py`): a
`r₀/ℓ = 0.8` la misma ecuación **relaja** con `m/m_B = 1` y **crece** con el
`m/m_B = 5/3` que fija el contratérmino de la ec. (10).

## La regla

| rutina | ec. | ¿`m/m_B` es libre? |
|---|---|---|
| `susceptibility` | (14) | **no lo lee**: la renormalización ya está dentro |
| `dispersion` | (\*) | lo lee, pero **debe quedarse en 1** para valer `ω²χ^r_ω` |
| `integrate_linear` | (\*) | no lo lee: divide por `1 − δm/m`, que es la misma renormalización |
| `integrate_delay` | (16) | **sí**, es el segundo eje del plano de fases |
| `integrate_memory` | (17) | **sí**, ídem |
| `block4_renorm` | (18) | **sí**: justamente lo busca, es la incógnita |

Dicho de otro modo: en el **bloque 1** la masa desnuda no es un parámetro, es
una consecuencia — `χ^r_ω` describe la teoría ya renormalizada. En los
**bloques 2 y 3** `m/m_B` es una coordenada, y el valor consistente con el
contratérmino es `m_over_mB_counterterm()`, es decir `1/(1 − δm/m)`, que es la
**línea punteada de la Fig. 2**.

La trampa está en que `dispersion` acepta el parámetro: pasarle el valor del
contratérmino *no* corrige nada, lo duplica, porque la resta `φ²/2` de
`_f_disp` ya quitó el contratérmino.

## Evidencia numérica

### 1. `dispersion` a `m/m_B = 1` es `χ`, y a cualquier otro valor no

Regulador suavizado, `r₀/ℓ = 3`, cero dominante:

| ruta | cero |
|---|---|
| `dispersion`, `m/m_B = 1` | `2.573620880600 − 0.717930440774 i` |
| `susceptibility` (ec. 14) | `2.573620879786 − 0.717930439793 i` |
| `dispersion`, `m/m_B = 2` (contratérmino) | `0 + 3.000000 i` |

Las dos primeras coinciden a **9 decimales**; la tercera es otra teoría.

### 2. La ec. (16) cambia de veredicto con el segundo eje

Regulador desplazado, mismo pulso, `s ≤ 25`:

| `r₀/ℓ` | `m/m_B` | resultado | tasa de los picos |
|---|---|---|---|
| 0.4 | 1.00 | relaja | — |
| 0.4 | 1.25 (contratérmino) | relaja | — |
| 0.8 | 1.00 | relaja | `−0.219` |
| 0.8 | 1.67 (contratérmino) | **crece** | `+0.264` |

Sobre la línea del contratérmino, el par `0.4` (estable) / `0.8` (inestable)
encierra el crítico `0.669` del bloque 1, como debe ser. Fuera de esa línea, no
hay razón para que lo haga.

Contraste con la vía espectral en el punto inestable (`r₀/ℓ = 0.8`, línea del
contratérmino):

| cantidad | ec. (16) integrada | ceros de `χ`, ec. (14) |
|---|---|---|
| tasa | `+0.264` | `Im ω = +0.230` |
| frecuencia | `π/0.4045 = 7.77` | `Re ω = 7.45` |

Mismo signo y mismo orden. La diferencia residual (15% en la tasa, 4% en la
frecuencia) es lo esperable: el pulso es finito y la (16) es no lineal, mientras
que la (14) es lineal e infinitesimal.

### 3. El paper mismo está sobre la línea del contratérmino

En la Fig. 1, Polonyi usa el regulador suavizado a `r₀/ℓ = 3` y barre
`m/m_B = 1.95, 1.98, 2.00`. El valor del contratérmino ahí es

    δm/m = r₀/(6ℓ) = 1/2   ⟹   m/m_B = 1/(1 − 1/2) = 2.000

exactamente. Las tres curvas del panel (a) rodean la línea punteada; el panel
(b), con `m/m_B = −3.80, −3.91, −4.10`, explora la rama `m_B < 0`. Es
confirmación independiente de que el eje vertical importa y de dónde cae la
línea.

## Qué NO queda afectado

La tabla de validación cruzada de `docs/ecuacion_linealizada.md` (tasa integrada
contra `Im ω` espectral, seis cutoffs) es consistente: `integrate_linear` divide
por `1 − δm/m` y `fig_trayectoria_lineal._cero_dominante` llama a `dispersion`
con `Params(ell=ell)`, o sea `m/m_B = 1`. Las dos rutas están en el mismo punto
—el renormalizado— y por eso coinciden a cuatro cifras. No hay nada que
corregir ahí.

Tampoco afecta a `seccion_A/fig_ec14_polos.py` ni a
`seccion_A/fig_ec15_vs_ec14.py`: los dos trabajan solo con `susceptibility`, que
no lee el parámetro.

## Reglas prácticas

1. **Al comparar bloque 1 contra bloque 2 o 3**, correr los bloques 2/3 en
   `m_over_mB = make_regulator(kind, ell).m_over_mB_counterterm()`. Con
   `m/m_B = 1` se está en otro punto del plano.
2. **Nunca pasar `m_over_mB ≠ 1` a `dispersion`** si lo que se quiere es
   `ω²χ^r_ω`. El parámetro existe para explorar la ecuación desnuda, no para
   corregir la renormalizada.
3. **Al reportar un resultado de estabilidad, dar siempre el par**
   `(r₀/ℓ, m/m_B)`. Solo el cutoff no identifica un punto.

<!-- ponytail: `dispersion` acepta un m_over_mB que casi siempre debe valer 1;
     invita al error que motivó este documento. La corrección mínima seria
     quitar el parametro y dejar la ec. (*) renormalizada como unica lectura,
     pero eso toca la API de nlaid/ y sus tests, asi que queda propuesto, no
     hecho. -->
