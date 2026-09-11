"""Figura del bloque 4: borde de estabilidad calculado vs. Fig. 2 del paper.

Panel A  regulador desplazado, ec. (16), contra el panel (a) de la Fig. 2.
Panel B  regulador suavizado,  ec. (17), contra el panel (b) de la Fig. 2.

La frontera del paper sale de data/fig2_vectorial.npz, leida del content
stream del PDF por scripts/digitaliza_fig2.py: son las coordenadas que
Mathematica escribio, sin error de lectura. (La digitalizacion por pixeles,
data/fig2_digitalizada.npz, queda como registro historico; la usa
figures/fig2b_verificacion_contraterm.png.) Los bordes propios salen de
scripts/scan_borde.py.

Ambas fronteras caen a r_0B/ell constante, con r_0B = r0 (m/m_B) el radio
clasico DESNUDO: el criterio de estabilidad compara el acoplamiento desnudo
con el cutoff. Por eso se anota esa constante en cada curva -- es la forma
invariante de comparar, independiente del cutoff.

Uso:  python3 scripts/fig_block4.py
"""

import json
import os

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid import RAIZ

MIO = "#2a78d6"        # calculo propio
PAPER = "#eb6834"      # digitalizacion del paper
INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"


def carga_mio(reg):
    ruta = RAIZ / f"data/borde_{reg}.json"
    if not os.path.exists(ruta):
        return None
    d = json.load(open(ruta))
    x, p, n = [], [], []
    for k, v in sorted(d.items(), key=lambda kv: float(kv[0])):
        x.append(float(k))
        p.append(v.get("positiva"))
        n.append(v.get("negativa"))
    return np.array(x), np.array([np.nan if q is None else q for q in p]), \
        np.array([np.nan if q is None else q for q in n])


def cte(x, y):
    """r_0B/ell = (r0/ell)(m/m_B): media y dispersion.

    Es constante SOLO si la frontera escala como 1/x. Se devuelve tambien la
    desviacion porque el panel (b) del paper no escala asi, y dar solo la media
    ahi sugiere una regularidad que no existe.
    """
    m = np.isfinite(y)
    if not m.any():
        return np.nan, np.nan
    p = x[m] * y[m]
    return float(np.mean(p)), float(np.std(p))


def carga_paper(pan, n=600):
    """Las dos ramas de la Fig. 2 del paper, sobre una malla comun.

    `sol0` es la rama m/m_B > 0 y `sol1` la negativa; cada una trae su propia
    malla en x, asi que se interpolan sobre una comun para poder sombrear
    entre ellas.
    """
    d = np.load(RAIZ / "data/fig2_vectorial.npz")
    s0, s1 = d[f"{pan}_sol0"], d[f"{pan}_sol1"]
    o0, o1 = np.argsort(s0[:, 0]), np.argsort(s1[:, 0])
    lo = max(s0[:, 0].min(), s1[:, 0].min(), 0.5)
    x = np.linspace(lo, min(s0[:, 0].max(), s1[:, 0].max()), n)
    return x, np.interp(x, s0[o0, 0], s0[o0, 1]), np.interp(x, s1[o1, 0], s1[o1, 1])


def main():
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))

    for ax, pan, reg, titulo in (
        (axes[0], "a", "shifted", "A.  desplazado, ec. (16)  vs  Fig. 2(a)"),
        (axes[1], "b", "smeared", "B.  suavizado, ec. (17)  vs  Fig. 2(b)"),
    ):
        xd, sup, inf = carga_paper(pan)
        ax.plot(xd, sup, lw=2, color=PAPER, label="paper (vectorial)")
        ax.plot(xd, inf, lw=2, color=PAPER)
        ax.fill_between(xd, inf, sup, color=PAPER, alpha=.10)

        mio = carga_mio(reg)
        if mio is not None:
            x, p, n = mio
            ax.plot(x, p, "o-", lw=2, ms=5, color=MIO, label="este código")
            ax.plot(x, n, "o-", lw=2, ms=5, color=MIO)
            for k, (c, sd) in enumerate((cte(x, p), cte(x, n))):
                ax.annotate(f"$r_{{0B}}/\\ell = {c:+.2f} \\pm {sd:.2f}$",
                            xy=(0.44, 0.93 - 0.055 * k), xycoords="axes fraction",
                            color=MIO, fontsize=9)
        for k, (c, sd) in enumerate((cte(xd, sup), cte(xd, inf))):
            ax.annotate(f"$r_{{0B}}/\\ell = {c:+.2f} \\pm {sd:.2f}$",
                        xy=(0.44, 0.82 - 0.055 * k), xycoords="axes fraction",
                        color=PAPER, fontsize=9)

        ax.axhline(0, color=GRID, lw=1.2, zorder=1)
        ax.set_xlabel("$r_0/\\ell$")
        ax.set_ylabel("$m/m_B$")
        ax.set_title(titulo, loc="left", fontsize=11, color=INK, pad=10)
        ax.set_xlim(0, 21)
        ax.set_ylim(-3.2, 3.2)
        ax.legend(frameon=False, fontsize=9, loc="lower right")
        ax.grid(True, color=GRID, lw=.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)

    fig.suptitle("Bloque 4 — borde de estabilidad: cálculo propio vs. Fig. 2 de "
                 "Polonyi (2019), leída del PDF como vectores", x=0.011, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salida = RAIZ / "figures/block4_borde.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)


if __name__ == "__main__":
    main()
