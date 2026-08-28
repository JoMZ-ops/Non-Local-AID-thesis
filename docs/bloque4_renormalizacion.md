# Bloque 4 — condición de renormalización y estructura de fases

## El problema con la ec. (18) tomada al pie de la letra

La condición de renormalización del paper es

    chi^r_omega = 1 + (2/3) i r0 omega                    (18)

Su **único cero** está en `omega = 3i/(2 r0)`, es decir en el semiplano
**superior**: es el polo runaway de Abraham-Lorentz. Una teoría que satisfaga
(18) exactamente no tiene ningún modo de relajación que monitorear.

Por tanto (18) no puede imponerse literalmente sobre una teoría estable. El
paper la aplica de forma operacional (p. 10):

> *"the renormalization condition, (18), can be satisfied by **monitoring the
> relaxation for large s**"*  [caso `m_B > 0`]
>
> *"The **relaxation of the envelope** is used to find the physical theory,
> obeying eq. (18) in this case"*  [caso `m_B < 0`]

Es decir: un procedimiento de *shooting*. Fijado `ell`, se integra la ecuación
no lineal, se mide la relajación a `s` grande, y se ajusta `m_B` hasta que
coincida con la predicción. Las líneas sólidas de la Fig. 2 son el resultado.

## Las dos condiciones implementadas

Como la lectura literal es imposible, se implementan dos condiciones que sí
están bien definidas, ambas por shooting sobre `m/m_B` a cutoff fijo:

**R1 — empalme con la teoría linealizada.** Tasa no lineal = `Im(omega)` del
cero dominante de `chi^r` (ec. 14). Es la lectura más fiel de "monitorear la
relajación": exige que la dinámica no lineal reproduzca la relajación que
predice la teoría linealizada ya renormalizada.

**R2 — borde de estabilidad.** Tasa no lineal = 0. Es la frontera de la región
sombreada de la Fig. 2, y el objeto más directamente comparable con el paper:
las tripletas de valores de la Fig. 1 bracketean justamente un cambio de signo
de la tasa.

## Por qué R2 biseca la clasificación y no la tasa

Cerca del borde de estabilidad la envolvente es plana **por definición**, de
modo que el ajuste exponencial pierde toda su capacidad de discriminación
justo en el punto que el shooting busca (`r2` cae a ~0.07 en el borde, frente
a ~0.99 a un 10% de distancia). Resolver `tasa = 0` con un ajuste que no es
informativo ahí produce raíces espurias.

Se biseca por tanto la **clasificación** estable/inestable, que sigue siendo
nítida en el borde. La clasificación usa tres señales robustas:

| señal | clasificación |
|---|---|
| deriva de `xdot^2` por encima del umbral | **inestable** (la integración se rompe porque la solución diverge) |
| medida fiable | el signo de la tasa decide |
| señal decaída al ruido de redondeo | **estable** (fuertemente amortiguada) |

Esta última distinción importa: una señal muerta y un borde de estabilidad
producen la misma firma superficial —una envolvente casi plana— y confundirlas
generaba raíces espurias en el propio punto inicial del bracket.

## Resultado: la frontera cae a `r_0B/ell` constante

El borde de estabilidad NO ocurre a `m/m_B` fijo sino a

    r_0B / ell = (r0/ell) * (m/m_B) = constante

donde `r_0B = e^2/(m_B c^2)` es el radio clásico **desnudo**. Es decir, el
criterio de estabilidad compara el acoplamiento desnudo con el cutoff, no con
`r0`. Para el regulador suavizado, ec. (17):

| `r0/ell` | 2 | 3 | 4 | 6 | 8 | 12 |
|---|---|---|---|---|---|---|
| `x (m/m_B)` rama + | 12.06 | 12.05 | 12.04 | 12.02 | 12.00 | 12.00 |
| `x (m/m_B)` rama − | −5.43 | −5.45 | −5.49 | −5.56 | — | — |

Constante al 0.5%. La Fig. 2 del paper tiene la **misma estructura**: sus
fronteras digitalizadas también dan `x·(m/m_B)` constante. Eso es un acuerdo
real sobre la *forma* de la frontera.

## Pero las constantes no coinciden

| | rama + | rama − |
|---|---|---|
| este código, suavizado ec. (17) | `+12.0` | `−5.5` |
| paper Fig. 2(b), suavizado | `+5.9` | `−10.7` |
| paper Fig. 2(a), desplazado | `+1.97` | `−3.52` |

Dos observaciones:

1. Los conjuntos de magnitudes `{12.0, 5.5}` y `{10.7, 5.9}` se parecen, pero
   con las ramas **intercambiadas**: en mi cálculo la rama de mayor magnitud
   es la positiva, en el paper es la negativa. Sugiere una diferencia de
   convención de signo en `m/m_B`, no un error de escala.
2. El cociente panel(b)/panel(a) del paper es `5.9/1.97 = 3.00` y
   `10.7/3.52 = 3.04`. **Exactamente 3**, en las dos ramas.

## El factor 3 de la Fig. 2 (ver también F1 en discrepancias.md)

Ajustando el factor de escala `c` que hace coincidir las dos fronteras
digitalizadas, `region_b(x) = region_a(x/c)`:

    c = 3.030,  RMS = 0.021      (con c = 1: RMS = 0.617)

Treinta veces mejor que sin reescalar. Y **3 es exactamente el factor que
iguala los contratérminos** de ambos reguladores: `delta_m` vale `r0/2ell`
(desplazado) y `r0/6ell` (suavizado), de modo que
`delta_m_desplazado(ell) = delta_m_suavizado(ell/3)`.

En mi cálculo independiente ese cociente **no** es 3: los cutoffs críticos de
la teoría linealizada (bloque 1) dan `4.001 / 0.669 = 5.98`, y los momentos de
los reguladores dan 3 para `delta_m` pero 4 para `I2`. Los dos reguladores no
están relacionados por un reescalado único de `ell`.

**Interpretación posible**: el panel (b) de la Fig. 2 podría ser el panel (a)
reescalado bajo el supuesto de que ambos reguladores son equivalentes al
igualar `delta_m`, en lugar de un cálculo suavizado independiente. Eso
encajaría con F1 (la línea punteada del panel (b) sigue la ley del
desplazado). **Es una hipótesis, no un hecho demostrado**: no se puede
descartar que la equivalencia sea real y que el desacuerdo esté en mi
implementación. Verificar antes de citar.
