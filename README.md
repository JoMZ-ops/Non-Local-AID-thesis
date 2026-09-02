# Non-Local AID — dinámica no local con regularización *point-splitting*

Implementación numérica de la dinámica efectiva de una carga puntual bajo
regularización por *point-splitting*, siguiendo J. Polonyi,
[arXiv:1701.04068v4](https://arxiv.org/abs/1701.04068),
*The Abraham-Lorentz force and electrodynamics at the classical electron radius*.

Punto de partida para la aplicación del método al modelo extendido de Starobinsky.

## Estructura

El código está organizado en **tres bloques que comparten un solo núcleo**, no
como tres programas independientes: un objeto `Params`, una jerarquía de
reguladores y una sola clase de historia de línea de mundo sirven a los tres.

| Módulo | Contenido |
|---|---|
| `nlaid/core.py` | métrica, `Params`, reguladores `delta_B` (ecs. 4 y 5) |
| `nlaid/worldline.py` | historia de línea de mundo, interpolación, punto retardado, fuente externa |
| `nlaid/block1_linear.py` | **Bloque 1** — susceptibilidad `chi^r_omega`, ec. (14), y búsqueda de ceros |
| `nlaid/block2_delay.py` | **Bloque 2** — regulador desplazado, ec. (16), retardo finito |
| `nlaid/block3_memory.py` | **Bloque 3** — regulador suavizado, ec. (17), memoria infinita |

Los tres bloques comparten `Params(ell, r0, m_over_mB, dim)`, de modo que cambiar
el cutoff o el regulador se propaga a todo el análisis sin duplicar parámetros.

## Convenciones

- Signatura `(+,-,-,-)`; líneas de mundo temporales tienen `x^2 > 0`.
- Unidades `c = 1`, `r0 = 1`: el cutoff `ell` se mide en unidades del radio
  clásico del electrón, y los ejes de la Fig. 2 del paper son
  `r0/ell` y `m/m_B` — exactamente los dos parámetros libres del modelo.
- Parametrización por tiempo propio `s`, con `xdot^2 = 1`.
- Por defecto `dim = 2` (1+1 dimensiones), que es donde el paper produce sus figuras.

## Instalación

```bash
pip install -e ".[dev]"
```

El `-e` importa: instala el paquete apuntando al repo, de modo que `import nlaid`
funciona desde cualquier directorio y los cambios al código se ven sin reinstalar.

### En Google Colab o Jupyter

```python
!git clone -b claude/nonlocal-memory-system-analysis-t9dnw3 https://github.com/JoMZ-ops/Non-Local-AID-thesis
!pip install -q -e Non-Local-AID-thesis
```

Los scripts localizan `figures/` y `data/` a partir de `nlaid.RAIZ`, no del
directorio de trabajo, así que corren igual desde donde sea:

```python
!python Non-Local-AID-thesis/scripts/fig_barrido_ell.py
```

o desde dentro del notebook. **`scripts/` no es un paquete instalable** — `pip`
solo instala `nlaid` —, así que un `import scripts...` falla con
`ModuleNotFoundError`. Lo correcto es la magia nativa de IPython/Colab, que
ejecuta el script en el espacio de nombres del notebook:

```python
%run Non-Local-AID-thesis/scripts/fig_barrido_ell.py
```

Después de eso, `compute`, `main` y `ELLS` quedan disponibles, y el barrido se
puede repetir con otros cutoffs sin editar el archivo:

```python
import numpy as np
main(ells=np.geomspace(4.0, 0.15, 7), force=True)
```

`notebooks/recorrido_completo.ipynb` ya trae la celda de arranque que hace todo
esto solo.

## Uso

```python
from nlaid.core import Params, make_regulator
from nlaid.block1_linear import dominant_pole, is_stable
from nlaid.block3_memory import integrate_memory

p = Params(ell=1/3, m_over_mB=-3.91)          # r0/ell = 3

# Bloque 1: ¿es estable la teoría linealizada?
reg = make_regulator("smeared", p.ell)
print(is_stable(reg, p), dominant_pole(reg, p))

# Bloque 3: integrar la dinámica completa
wl = integrate_memory(p, s_end=25.0, ds=5e-3)
print("deriva de xdot^2:", wl.norm_drift.max())
```

Figuras:

```bash
python scripts/fig_block1.py      # estabilidad linealizada + contraterm
python scripts/fig_blocks23.py    # reproducción de la Fig. 1 del paper
```

## Validación

`pytest` (41 tests). La estrategia es contrastar contra **valores cerrados
derivados a mano**, no contra otra corrida del propio código:

- Momentos analíticos de los reguladores (`int delta_B = 1`, `delta_m/m = r0/6ell`,
  `I2 = 2 ell`).
- Forma cerrada de `chi^r_omega` para el regulador desplazado.
- El polo *runaway* de Abraham-Lorentz en `omega = 3i/2` al remover el cutoff (ec. 15).
- Conteo de ceros por principio del argumento contrastado con la búsqueda por Muller.
- **Ortogonalidad exacta** `xddot . xdot = 0` de los lados derechos de las ecs.
  (16) y (17): identidad analítica que el código cumple a `1e-15`, y el detector
  más sensible de un error de transcripción.
- Orden de convergencia 2 en `ds` para ambos integradores; convergencia
  independiente en la ventana de memoria y en la cuadratura del bloque 3.

Los tests marcados `slow` son los de convergencia: `pytest -m "not slow"` los omite.

## Notas de lectura

`docs/notas_paper.md` fija la transcripción de las ecuaciones, los resultados
analíticos derivados para los tests, y la clasificación matemática de cada
ecuación. **Verificar contra el PDF original antes de citar en la tesis.**
