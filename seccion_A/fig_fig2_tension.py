"""La Fig. 2 del paper, digitalizada por vectores, contra nuestro calculo.

Responde la hipotesis 4 de `docs/discrepancias.md`, que estaba anotada como
"lectura mia de una figura de baja resolucion; verificar antes de sacar
conclusiones". Verificada: los datos salen del content stream del PDF, no de
pixeles (ver `scripts/digitaliza_fig2.py`).

Panel A  las dos lineas punteadas superpuestas. Son la MISMA ruta vectorial en
         los dos paneles -- mismo numero de puntos y max|dif| = 0 exacto -- y
         siguen 1/(1 - x/2), el contratermino del regulador DESPLAZADO. En el
         panel (b), que es el suavizado, deberia ser 1/(1 - x/6). Es F1,
         ahora sin margen de lectura.
Panel B  el borde de estabilidad del panel (a), desplazado, contra el nuestro:
         misma forma hiperbolica y mismo orden, razon 0.6 a 1.1.
Panel C  lo mismo para el panel (b), suavizado: razon 3.4 a 300, y ni siquiera
         la misma forma. Ademas, el valor que la Fig. 1 usa a r0/ell = 3
         (m/m_B ~ 2) cae FUERA del eje vertical de la Fig. 2, que llega a
         0.853.

Necesita `data/fig2_vectorial.npz` (de `scripts/digitaliza_fig2.py`) y
`data/borde_*.json` (de `scripts/scan_borde.py`).

Uso:  python3 seccion_A/fig_fig2_tension.py
"""

import json
import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid import RAIZ

VECT = RAIZ / "data/fig2_vectorial.npz"
INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"
AZUL, NARANJA, ROJO, VERDE = "#2a78d6", "#eb6834", "#c1442e", "#2e7d5b"

# Los tres valores de la Fig. 1(a), a r0/ell = 3.
FIG1 = [1.95, 1.98, 2.00]


def borde(kind):
    d = json.load(open(RAIZ / f"data/borde_{kind}.json"))
    x = np.array([float(k) for k in d])
    return x, np.array([d[k]["positiva"] for k in d], dtype=float)


def panel_punteada(ax, d):
    for pan, c, lw, ls in (("a", AZUL, 5.0, "-"), ("b", NARANJA, 1.8, "-")):
        for k in (0, 1):
            p = d[f"{pan}_pun{k}"]
            ax.plot(p[:, 0], p[:, 1], lw=lw, ls=ls, color=c, alpha=.75,
                    label=f"punteada del panel ({pan})" if k == 1 else None)

    # La ley suavizada tiene su polo en x = 6: se corta, o matplotlib une las
    # dos ramas con una vertical que parece una curva y no lo es.
    x = np.linspace(2.05, 20, 400)
    y = 1 / (1 - x / 6)
    y[np.abs(1 - x / 6) < 0.08] = np.nan
    ax.plot(x, y, lw=1.6, ls=(0, (5, 2)), color=VERDE,
            label="$1/(1-x/6)$: contratérmino suavizado")
    ax.annotate("las dos punteadas son la MISMA ruta\n"
                "vectorial: $\\max|\\Delta| = 0$ exacto",
                xy=(8.5, -0.92), fontsize=8.5, color=INK)
    ax.annotate("lo que el panel (b)\ndebería mostrar", xy=(13.5, -0.52),
                fontsize=8.5, color=VERDE)

    ax.axhline(0, color=GRID, lw=1.2, zorder=0)
    ax.set_xlabel("$r_0/\\ell$")
    ax.set_ylabel("$m/m_B$")
    ax.set_title("A.  La línea punteada, ec. (10)\n"
                 "ambos paneles dibujan la ley del regulador desplazado",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, 20.5)
    ax.set_ylim(-1.6, 1.2)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")


def panel_borde(ax, d, pan, kind, titulo, nota):
    p = d[f"{pan}_sol0"]                       # rama m/m_B > 0
    o = np.argsort(p[:, 0])
    ax.plot(p[o, 0], p[o, 1], lw=2.2, color=AZUL, label=f"Fig. 2({pan}) del paper")

    xb, yb = borde(kind)
    ax.plot(xb, yb, "s--", lw=1.8, ms=5, color=NARANJA,
            label="nuestro borde (tasa = 0)")

    a = d[f"{pan}_area"]
    ax.axhline(a[3], color=MUTED, lw=1, ls=(0, (2, 3)))
    ax.annotate(f"tope del eje de la Fig. 2: {a[3]:.3f}", xy=(0.4, a[3] + 0.06),
                fontsize=8, color=MUTED)

    if kind == "smeared":
        ax.plot([3.0] * 3, FIG1, "o", ms=6, color=ROJO, zorder=6)
        ax.annotate("los tres valores de la Fig. 1(a)\n"
                    "($m/m_B$ = 1.95, 1.98, 2.00)\nno caben en la Fig. 2",
                    xy=(3.0, 1.98), xytext=(5.2, 2.6), fontsize=8.5, color=ROJO,
                    arrowprops=dict(arrowstyle="->", color=ROJO, lw=.9))

    ax.annotate(nota, xy=(0.4, 0.06 if kind == "shifted" else 3.9),
                fontsize=8.5, color=INK)
    ax.set_xlabel("$r_0/\\ell$")
    ax.set_ylabel("$m/m_B$")
    ax.set_title(titulo, loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, 20.5)
    ax.set_ylim(0, 0.9 if kind == "shifted" else 4.6)
    ax.legend(frameon=False, fontsize=8.5,
              loc="upper right" if kind == "shifted" else "center right")


def main():
    if not os.path.exists(VECT):
        sys.exit(f"falta {VECT}: corra antes scripts/digitaliza_fig2.py <pdf>")
    d = dict(np.load(VECT))

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8))
    panel_punteada(axes[0], d)
    panel_borde(axes[1], d, "a", "shifted",
                "B.  Panel (a), desplazado: coincide\n"
                "misma forma $\\sim C/x$, razón 0.6 a 1.1",
                "razón paper/nuestro entre\n0.58 y 1.07")
    panel_borde(axes[2], d, "b", "smeared",
                "C.  Panel (b), suavizado: no coincide\n"
                "razón 3.4 a 300, y ni siquiera la misma forma",
                "el nuestro es $12/x$ exacto;\nel del paper sube y baja")

    for ax in axes:
        ax.grid(True, color=GRID, lw=.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK)
        ax.yaxis.label.set_color(INK)

    fig.suptitle("La Fig. 2 del paper leída del PDF como vectores: qué coincide "
                 "y qué no", x=0.008, ha="left", fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salida = RAIZ / "figures/seccionA_fig2_tension.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)

    # Las razones que se citan en los titulos, impresas para poder auditarlas.
    for pan, kind in (("a", "shifted"), ("b", "smeared")):
        p = d[f"{pan}_sol0"]
        o = np.argsort(p[:, 0])
        xb, yb = borde(kind)
        m = (xb >= p[o, 0].min()) & (xb <= p[o, 0].max())
        r = yb[m] / np.interp(xb[m], p[o, 0], p[o, 1])
        print(f"  panel ({pan}) {kind:8s}: razón nuestro/paper de "
              f"{r.min():.2f} a {r.max():.2f}")


if __name__ == "__main__":
    main()
