"""Figura del bloque 4: borde de estabilidad calculado vs. Fig. 2 del paper.

Panel A  regulador desplazado, ec. (16), contra el panel (a) de la Fig. 2.
Panel B  regulador suavizado,  ec. (17), contra el panel (b) de la Fig. 2.

La frontera digitalizada del paper esta en data/fig2_digitalizada.npz
(calibracion documentada en data/fig2_digitalizada.json). Los bordes propios
salen de scripts/scan_borde.py.

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

MIO = "#2a78d6"        # calculo propio
PAPER = "#eb6834"      # digitalizacion del paper
INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"


def carga_mio(reg):
    ruta = f"data/borde_{reg}.json"
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
    """r_0B/ell = (r0/ell)(m/m_B), constante si la frontera escala como 1/x."""
    m = np.isfinite(y)
    return float(np.mean(x[m] * y[m])) if m.any() else np.nan


def main():
    dig = np.load("data/fig2_digitalizada.npz")
    fig, axes = plt.subplots(1, 2, figsize=(11.5, 4.8))

    for ax, pan, reg, titulo in (
        (axes[0], "a", "shifted", "A.  desplazado, ec. (16)  vs  Fig. 2(a)"),
        (axes[1], "b", "smeared", "B.  suavizado, ec. (17)  vs  Fig. 2(b)"),
    ):
        xd, sup, inf = dig[f"{pan}_x"], dig[f"{pan}_sup"], dig[f"{pan}_inf"]
        libre = (sup < 0.80) & (inf > -1.22)        # no recortado por el marco
        ax.plot(xd[libre], sup[libre], lw=2, color=PAPER, label="paper (digitalizado)")
        ax.plot(xd[libre], inf[libre], lw=2, color=PAPER)
        ax.fill_between(xd[libre], inf[libre], sup[libre], color=PAPER, alpha=.10)

        mio = carga_mio(reg)
        if mio is not None:
            x, p, n = mio
            ax.plot(x, p, "o-", lw=2, ms=5, color=MIO, label="este código")
            ax.plot(x, n, "o-", lw=2, ms=5, color=MIO)
            cp, cn = cte(x, p), cte(x, n)
            ax.annotate(f"$r_{{0B}}/\\ell = {cp:+.2f}$", xy=(0.52, 0.90),
                        xycoords="axes fraction", color=MIO, fontsize=9)
            ax.annotate(f"$r_{{0B}}/\\ell = {cn:+.2f}$", xy=(0.52, 0.83),
                        xycoords="axes fraction", color=MIO, fontsize=9)
        cps, cns = cte(xd[libre], sup[libre]), cte(xd[libre], inf[libre])
        ax.annotate(f"$r_{{0B}}/\\ell = {cps:+.2f}$", xy=(0.52, 0.74),
                    xycoords="axes fraction", color=PAPER, fontsize=9)
        ax.annotate(f"$r_{{0B}}/\\ell = {cns:+.2f}$", xy=(0.52, 0.67),
                    xycoords="axes fraction", color=PAPER, fontsize=9)

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

    fig.suptitle("Bloque 4 — borde de estabilidad: cálculo propio vs. Fig. 2 "
                 "de Polonyi (2019)", x=0.011, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig("figures/block4_borde.png", dpi=170, facecolor="#fcfcfb")
    print("figures/block4_borde.png")


if __name__ == "__main__":
    main()
