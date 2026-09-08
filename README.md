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
| `scripts/`, `seccion_A/` | figuras. **No son paquetes instalables**: se corren, no se importan |

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
!git clone -b Electrodynamics---Abraham-\&-Lorentz-Force https://github.com/JoMZ-ops/Non-Local-AID-thesis
!pip install -q -e Non-Local-AID-thesis
```

(El `\&` es por el shell: sin la barra, `&` manda el comando al fondo y el clon
se corta a la mitad. La rama de desarrollo es
`Electrodynamics---Abraham-&-Lorentz-Force`.)

Los scripts localizan `figures/` y `data/` a partir de `nlaid.RAIZ`, no del
directorio de trabajo, y **no lo modifican**, así que corren igual desde donde
sea sin dejar rastro en tu sesión:

```python
!python Non-Local-AID-thesis/scripts/fig_barrido_ell.py
```

o desde dentro del notebook. **`scripts/` no es un paquete instalable** — `pip`
solo instala `nlaid` —, así que un `import scripts...` falla con
`ModuleNotFoundError`. Lo correcto es la magia nativa de IPython/Colab, que
ejecuta el script en el espacio de nombres del notebook:

```python
%run Non-Local-AID-thesis/scripts/fig_barrido_ell.py
%run Non-Local-AID-thesis/seccion_A/fig_ec14_polos.py
```

**Instalas una vez y a partir de ahí la ruta relativa vale siempre**: ningún
script cambia el directorio de trabajo. Las rutas de salida se construyen desde
`nlaid.RAIZ` (`fig.savefig(RAIZ / "figures/...")`), así que el archivo aterriza
en el repo aunque invoques el script desde donde sea, y tu `cwd` queda intacto.

> Hasta el commit `fc017d9` cada script hacía `os.chdir(RAIZ)` al arrancar. El
> efecto secundario era que, tras el primer `%run`, el `cwd` del notebook pasaba
> a ser la raíz del repo y el segundo `%run` con ruta relativa fallaba con
> `Exception: File ... not found`. Si ves ese error, tu clon es anterior al
> arreglo: vuelve a clonar, o usa `f"{nlaid.RAIZ}/scripts/..."`, que es absoluta
> y funciona en cualquier versión.

Después de eso, `compute`, `main` y `ELLS` quedan disponibles, y el barrido se
puede repetir con otros cutoffs sin editar el archivo:

```python
import numpy as np
main(ells=np.geomspace(4.0, 0.15, 7), force=True)
```

`notebooks/recorrido_completo.ipynb` ya trae la celda de arranque que hace todo
esto solo.

### Por qué no ves las gráficas aparecer

**No es un error tuyo: es por diseño.** Todos los scripts empiezan con

```python
import matplotlib
matplotlib.use("Agg")
```

`Agg` es el backend *sin pantalla*: dibuja a un archivo, nunca a una ventana ni
a la salida del notebook. Ninguno de los scripts llama a `plt.show()`; todos
terminan en `fig.savefig(RAIZ / "figures/...png")` y anuncian la ruta. Es lo correcto
para un script que también debe correr por línea de comandos o en CI, donde no
hay pantalla que valga.

Para **ver** la figura en el notebook, muéstrala después de generarla:

```python
%run Non-Local-AID-thesis/seccion_A/fig_ec7_ec10.py
from IPython.display import Image
Image("Non-Local-AID-thesis/figures/seccionA_ec7_ec10.png")
```

Dos avisos prácticos:

1. Después de un `%run` de cualquiera de estos scripts, el backend del notebook
   queda en `Agg` y **tus propias gráficas inline dejan de aparecer**. Se
   arregla volviendo a ejecutar `%matplotlib inline`.
2. Para trabajar de verdad — probar, cambiar un parámetro, volver a mirar — no
   uses los scripts: usa el paquete. Ahí sí hay figura inline, porque la armas
   tú:

```python
%matplotlib inline
import matplotlib.pyplot as plt, numpy as np
from nlaid.core import Params, make_regulator
from nlaid.block1_linear import susceptibility

w = np.linspace(0.02, 9, 400)
for x in (2, 3, 4):                       # r0/ell
    ell = 1/x
    chi = susceptibility(w, make_regulator("smeared", ell), Params(ell=ell))
    plt.plot(w, chi.imag, label=f"$r_0/\\ell$ = {x}")
plt.xlabel("$\\omega$"); plt.ylabel("Im $\\chi^r_\\omega$"); plt.legend()
plt.show()
```

Los scripts de `figures/` son el **producto final**; `nlaid/` es la
**herramienta**. Para explorar, importa la herramienta.

`notebooks/seccion_A.ipynb` recorre así las ecs. (7), (10), (14), (15) y (16):
arma cada gráfica en el notebook, con celdas cortas para cambiar un parámetro y
volver a mirar. Es el punto de entrada recomendado si estás empezando.

### Dónde editar el código

| quiero… | dónde | por qué |
|---|---|---|
| correr y experimentar en celdas | Colab / Jupyter | rápido, sin instalar nada |
| **cambiar el código** | clon local + `pip install -e` | el editable install hace que el cambio se vea sin reinstalar |
| corregir una línea de documentación | GitHub web | no vale la pena clonar |

Lo que **no** conviene es copiar y pegar el contenido de los archivos en celdas:
se pierde el historial, y dos copias del mismo código divergen en una tarde.

En Colab los archivos del clon son **efímeros**: al cerrar la sesión se borran.
Sirve para correr, no para editar. Si igual quieres editar ahí, abre el archivo
con el explorador de archivos (panel izquierdo, doble clic) y guarda con
`Ctrl+S`: el cambio corre de inmediato, pero muere con la sesión si no lo
llevas al repositorio — ver **«Editaste en Colab y funciona: ¿cómo lo
guardas?»** más abajo.

Con la instalación editable (`pip install -e`), tras cambiar un archivo de
`nlaid/` basta con reiniciar el kernel, o mejor:

```python
%load_ext autoreload
%autoreload 2
```

y a partir de ahí cada celda vuelve a leer el archivo modificado sola.

### Editaste en Colab y funciona: ¿cómo lo guardas?

El clon de Colab es una copia efímera. Ya verificaste el cambio ahí; ahora hay
que llevarlo al repositorio. Tres rutas.

**1. Descargar y subir a GitHub.** La más simple, sin credenciales.

Panel de archivos de Colab → clic derecho sobre el archivo → *Descargar*. En
GitHub, en la rama de desarrollo: *Add file → Upload files*, arrastras el
archivo **dentro de la carpeta que le corresponde** (`scripts/`, `nlaid/`…),
mensaje de commit, *Commit changes*.

Sube el archivo entero, así que no hay riesgo de pegar de más o de menos, y
**funciona con los PNG de `figures/`**, que son binarios.

**2. Copiar y pegar en GitHub web.** Solo para archivos de texto: `Ctrl+A`,
`Ctrl+C` en el editor de Colab; en GitHub, ✏️ sobre el archivo, `Ctrl+A`,
`Ctrl+V`, commit. Más frágil que la ruta 1 y no sirve para figuras.

**3. `git push` desde el propio Colab.** Requiere un *personal access token* con
permiso de escritura (*Contents*) sobre el repo. Primero, actualizar el clon y
commitear:

```python
!git -C Non-Local-AID-thesis config user.name  "Tu Nombre"
!git -C Non-Local-AID-thesis config user.email "tu@correo"
!git -C Non-Local-AID-thesis add -A
!git -C Non-Local-AID-thesis commit -m "fig_masa: marca el cutoff critico"
!git -C Non-Local-AID-thesis pull --rebase
```

El `pull --rebase` no es opcional: si la rama avanzó desde que clonaste, el push
se rechaza. Y luego:

```python
import getpass, subprocess
tok = getpass.getpass("Token de GitHub: ")      # no se ve ni queda en el .ipynb
subprocess.run(["git", "-C", "Non-Local-AID-thesis", "push",
                f"https://{tok}@github.com/JoMZ-ops/Non-Local-AID-thesis.git",
                "HEAD:Electrodynamics---Abraham-&-Lorentz-Force"], check=True)
```

Tres precauciones, porque un token es una credencial:

1. **Nunca lo escribas literal en una celda**: el notebook se guarda con todo lo
   que teclaste dentro. De ahí `getpass`.
2. `subprocess.run` con lista, **no** `!git push https://TOKEN@...`: con `!` el
   comando pasa por el shell y el token puede quedar en el historial y en la
   salida de la celda.
3. Si el push falla, el error **puede incluir la URL con el token**. Borra la
   salida de esa celda si ocurre.

| situación | ruta |
|---|---|
| un cambio suelto, de vez en cuando | **1** |
| vas a editar seguido desde Colab | **3** (se configura una vez) |
| trabajo sostenido | clon local: `git push` sin ceremonia |

**Antes de subir, mira el resultado, no solo que no haya fallado.** «Corrió sin
error» y «quedó bien» son cosas distintas — una anotación puede salirse del
panel sin que Python se queje:

```python
from IPython.display import Image
Image("Non-Local-AID-thesis/figures/masa_vs_cutoff.png")
```

Si tocaste `nlaid/`, además `!python3 -m pytest -q Non-Local-AID-thesis`.

Y commitea el código junto con la figura que produce, en el mismo commit: así el
repo nunca tiene un PNG que no corresponde al código que lo generó.

(Montar Google Drive y clonar ahí evita reclonar en cada sesión, pero **no**
guarda nada en GitHub: para eso sigue haciendo falta una de las tres rutas.)

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
python scripts/fig_block1.py            # estabilidad linealizada + contraterm
python scripts/fig_blocks23.py          # reproducción de la Fig. 1 del paper
python seccion_A/fig_ec7_ec10.py        # integrandos de las ecs. (7) y (10)
python seccion_A/fig_ec14_polos.py      # ceros de chi en el plano complejo
python seccion_A/fig_ec15_vs_ec14.py    # el límite local contra la (14)
python seccion_A/fig_ec16_retardo.py    # geometría del retardo, ec. (16)
```

Cada uno imprime la ruta del PNG que escribió en `figures/`.

## Validación

`pytest` (66 tests). La estrategia es contrastar contra **valores cerrados
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
