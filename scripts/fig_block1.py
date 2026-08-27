"""Figura del Bloque 1: estabilidad de la teoria linealizada, ec. (14).

Panel A  Im(omega) del cero dominante de chi^r contra r0/ell. La asintota
         3/2 es el polo runaway de Abraham-Lorentz de la ec. (15).
Panel B  Linea de contraterm m/m_B de la ec. (10): la linea punteada de la
         Fig. 2 del paper.

Uso:  python3 scripts/fig_block1.py
"""

import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid.core import Params, make_regulator
from nlaid.block1_linear import critical_cutoff, dominant_pole

warnings.simplefilter("ignore")

# Paleta categorica en orden fijo (validada para daltonismo: dE 24.7 protan).
COLOR = {"shifted": "#2a78d6", "smeared": "#eb6834"}
LABEL = {"shifted": "desplazado  ec. (4)", "smeared": "suavizado  ec. (5)"}
INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"

BOX = dict(re_max=12.0, im_hi=12.0)


CACHE = "figures/block1_datos.npz"


def compute(force=False):
    """Barrido de polos. Se cachea: el calculo domina el costo del replot."""
    import os
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        return d["xs"], {k: d[k] for k in ("shifted", "smeared")}, \
               {k: float(d["crit_" + k]) for k in ("shifted", "smeared")}

    xs = np.concatenate([np.linspace(0.3, 6, 40), np.linspace(6.5, 30, 40)])
    poles, crit = {}, {}
    for kind in ("shifted", "smeared"):
        poles[kind] = np.array([
            (lambda z: z.imag if z is not None else np.nan)(
                dominant_pole(make_regulator(kind, 1.0 / x), Params(ell=1.0 / x),
                              grid=180, **BOX))
            for x in xs])
        crit[kind] = critical_cutoff(kind, lo=0.3, hi=30.0, tol=2e-3,
                                     n=1200, **BOX)
        print(f"{kind:9s}  r0/ell critico = {crit[kind]:.3f}"
              f"   (ell_c = {1/crit[kind]:.3f} r0)")
    np.savez(CACHE, xs=xs, **poles,
             **{"crit_" + k: v for k, v in crit.items()})
    return xs, poles, crit


def main():
    xs, poles, crit = compute()

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    # --- Panel A ---------------------------------------------------------
    axA.axhline(1.5, color=MUTED, lw=1, ls=(0, (5, 4)), zorder=1)
    axA.annotate("polo de Abraham-Lorentz,  Im $\\omega = 3/2$   (ec. 15)",
                 xy=(9.5, 1.30), fontsize=8.5, color=MUTED)
    for kind in ("shifted", "smeared"):
        m = np.isfinite(poles[kind])
        axA.plot(xs[m], poles[kind][m], lw=2, color=COLOR[kind],
                 label=LABEL[kind], zorder=3)
        axA.axvline(crit[kind], color=COLOR[kind], lw=1, ls=":", alpha=.7, zorder=2)
        j = np.where(m)[0][-1]
        axA.annotate(LABEL[kind].split()[0], xy=(xs[j], poles[kind][j]),
                     xytext=(6, 3 if kind == "smeared" else -9),
                     textcoords="offset points",
                     color=COLOR[kind], fontsize=9, va="center")
        axA.annotate(f"$r_0/\\ell$ critico = {crit[kind]:.2f}",
                     xy=(crit[kind], 3.95 if kind == "shifted" else 3.55),
                     xytext=(5, 0), textcoords="offset points",
                     color=COLOR[kind], fontsize=8.5)
    axA.set_xlabel("$r_0/\\ell$   (cutoff $\\Lambda = 1/\\ell$ creciente $\\rightarrow$)")
    axA.set_ylabel("Im $\\omega$  del cero dominante de $\\chi^r_\\omega$")
    axA.set_title("A.  Tasa de crecimiento de la inestabilidad", loc="left",
                  fontsize=11, color=INK, pad=10)
    axA.set_xlim(0, 34); axA.set_ylim(0, 4.4)
    axA.legend(frameon=False, fontsize=9, loc="upper right",
               bbox_to_anchor=(1.0, 0.78))

    # --- Panel B ---------------------------------------------------------
    xb = np.linspace(0.3, 30, 900)
    for kind in ("shifted", "smeared"):
        y = np.array([make_regulator(kind, 1.0 / x).m_over_mB_counterterm()
                      for x in xb])
        y[np.abs(y) > 6] = np.nan          # corta la divergencia delta_m = m
        axB.plot(xb, y, lw=2, color=COLOR[kind], label=LABEL[kind])
    axB.axhline(0, color=GRID, lw=1.2, zorder=1)
    axB.annotate("$m_B > 0$", xy=(24, 0.55), fontsize=9, color=MUTED)
    axB.annotate("$m_B < 0$", xy=(24, -0.75), fontsize=9, color=MUTED)
    axB.set_xlabel("$r_0/\\ell$")
    axB.set_ylabel("$m/m_B$")
    axB.set_title("B.  Contraterm de masa, ec. (10)\n"
                  "(linea punteada de la Fig. 2 del paper)", loc="left",
                  fontsize=11, color=INK, pad=10)
    axB.set_xlim(0, 34); axB.set_ylim(-1.6, 1.6)
    for kind, yy in (("shifted", -0.18), ("smeared", -0.60)):
        axB.annotate(LABEL[kind].split()[0], xy=(30.2, yy),
                     color=COLOR[kind], fontsize=9, va="center")
    axB.legend(frameon=False, fontsize=9, loc="lower left")

    for ax in (axA, axB):
        ax.grid(True, color=GRID, lw=.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)

    fig.suptitle("Bloque 1 — estabilidad linealizada del regulador point-splitting "
                 "(Polonyi 2019, ec. 14)", x=0.011, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.94))
    fig.savefig("figures/block1_estabilidad.png", dpi=170,
                facecolor="#fcfcfb")
    print("\nfigures/block1_estabilidad.png")


if __name__ == "__main__":
    main()
