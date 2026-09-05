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

## Borde de estabilidad R2: barrido completo en el cutoff

`scripts/scan_borde.py` resuelve, para cada cutoff, el `m/m_B` donde la tasa
de relajación se anula, en las dos ramas de signo de `m_B`. Datos en
`data/borde_smeared.json` y `data/borde_shifted.json`.

Suavizado, ec. (5) — `s_end = 14`, `ds = 1e-2`:

| r0/ell | m/m_B (+) | m/m_B (−) | r_0B/ell (+) | r_0B/ell (−) |
|---:|---:|---:|---:|---:|
|  2 |  6.0280 | −2.7150 | 12.056 | −5.430 |
|  3 |  4.0150 | −1.8180 | 12.045 | −5.454 |
|  4 |  3.0090 | −1.3720 | 12.036 | −5.488 |
|  6 |  2.0030 | −0.9260 | 12.018 | −5.556 |
|  8 |  1.5000 | −0.7115 | 12.000 | −5.692 |
| 10 |  1.2005 | −0.5719 | 12.005 | −5.719 |
| 12 |  1.0000 | −0.4804 | 12.000 | −5.765 |

Desplazado, ec. (4) — `s_end = 10`, `ds = 5e-3`:

| r0/ell | m/m_B (+) | m/m_B (−) | r_0B/ell (+) | r_0B/ell (−) |
|---:|---:|---:|---:|---:|
|  2 | 0.6267 | −0.4900 | 1.253 | −0.980 |
|  3 | 0.3412 | −0.3220 | 1.024 | −0.966 |
|  4 | 0.2530 | −0.2447 | 1.012 | −0.979 |
|  6 | 0.1704 | −0.1648 | 1.022 | −0.989 |
|  8 | 0.1259 | −0.1234 | 1.007 | −0.987 |
| 12 | 0.0878 | −0.0836 | 1.053 | −1.004 |

`r_0B = r0 (m/m_B)` es el radio clásico **desnudo**.

### Lo que se puede afirmar

1. **La frontera es una condición sobre `r_0B/ell`, no sobre `r0/ell`.** En las
   cuatro ramas el producto es estacionario mientras `m/m_B` recorre casi dos
   órdenes de magnitud (de 6.03 a 0.084). El criterio de estabilidad compara
   el acoplamiento desnudo con el cutoff.
2. **Suavizado, rama positiva: `r_0B/ell = 12.00`**, alcanzado
   monótonamente desde arriba. Es el resultado más limpio del barrido.
3. **Desplazado: `|r_0B/ell| ≈ 1` en ambas ramas**, con dispersión compatible
   con la tolerancia de bisección (`tol = 3e-3` en `m/m_B`, que a `r0/ell = 12`
   son 0.036 en el producto). El punto `r0/ell = 2` es el outlier de la rama
   positiva (1.253) y es exactamente donde `delta_m = m` para este regulador:
   el contratérmino diverge ahí, así que el punto está en el borde del dominio
   y no debe promediarse con los demás.

### Lo que queda abierto

- **La asimetría entre ramas del suavizado.** La rama negativa deriva de −5.43
  a −5.77 sin estacionarse. Extrapolando en `1/x` sobre los últimos 3–5 puntos
  el límite sale entre −5.89 y −5.97, **compatible con −6** (y por tanto con
  una razón exacta −2 frente a la rama positiva), pero la dispersión de los
  ajustes es del mismo orden que la separación respecto de −6. **No es
  concluyente a la tolerancia actual.** Lo decidiría bajar `tol` a 1e-4 y
  añadir `r0/ell = 20, 30`.
- **Por qué el desplazado es casi simétrico (+1.02 / −0.99) y el suavizado no
  (+12.0 / −5.8).** Sin explicación por ahora.

### Trampa operativa encontrada durante este barrido

`pip install -e .` había quedado apuntando a un clon desechable en `/tmp`. Como
`python3 scripts/foo.py` pone `scripts/` en `sys.path[0]` y **no** el directorio
de trabajo, `import nlaid` resolvía al paquete instalado y `RAIZ` apuntaba al
clon: los resultados se escribían allí en silencio. Los `.py` eran idénticos, de
modo que los números eran válidos, pero el `.json` del repo no se actualizaba.
`scan_borde.py` ahora imprime la ruta de salida al arrancar.
