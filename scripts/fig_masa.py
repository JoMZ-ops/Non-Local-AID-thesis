"""Impacto del cutoff sobre la renormalizacion de masa, ec. (10).

    delta_m = (e^2 / 2 c^2) int_0^inf dz z^{-1/2} delta_B(z)
    delta_m / m = (r0/2) int_0^inf dz z^{-1/2} delta_B(z)

En forma cerrada:  r0/(6 ell) suavizado,  r0/(2 ell) desplazado.

Panel A  delta_m/m contra r0/ell. Divergencia LINEAL al remover el cutoff.
         El cruce delta_m = m marca donde el bare mass cambia de signo.
Panel B  m_B/m = 1 - delta_m/m. Se grafica esta y no m/m_B porque no diverge:
         cruza cero limpiamente en r0/ell = 6 (suavizado) y 2 (desplazado),
         que es la asintota vertical de m/m_B y separa las dos ramas de la
         Fig. 2 del paper.

Uso:  python3 scripts/fig_masa.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid.core import make_regulator

COLOR = {"shifted": "#2a78d6", "smeared": "#eb6834"}
LABEL = {"shifted": "desplazado  ec. (4)", "smeared": "suavizado  ec. (5)"}
INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"

# r0/ell donde delta_m = m, es decir donde m_B cambia de signo.
CRUCE = {"shifted": 2.0, "smeared": 6.0}


def main():
    x = np.linspace(0.05, 20, 800)          # r0/ell
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.4))

    for kind in ("shifted", "smeared"):
        dm = np.array([make_regulator(kind, 1.0 / xi).mass_shift_over_m()
                       for xi in x])
        axA.plot(x, dm, lw=2, color=COLOR[kind], label=LABEL[kind])
        axB.plot(x, 1.0 - dm, lw=2, color=COLOR[kind], label=LABEL[kind])
        axB.plot([CRUCE[kind]], [0.0], "o", ms=7, color=COLOR[kind],
                 mec="#fcfcfb", mew=1.5, zorder=5)
        axB.annotate(f"$r_0/\\ell = {CRUCE[kind]:.0f}$",
                     xy=(CRUCE[kind], 0.0), xytext=(6, 10),
                     textcoords="offset points", color=COLOR[kind], fontsize=9)

    axA.axhline(1.0, color=MUTED, lw=1, ls=(0, (5, 4)))
    axA.annotate("$\\delta m = m$   $\\Rightarrow$   $m_B$ cambia de signo",
                 xy=(9.5, 1.12), fontsize=9, color=MUTED)
    axA.set_ylabel("$\\delta m / m$")
    axA.set_title("A.  Renormalización de masa, ec. (10)\n"
                  "divergencia lineal en $1/\\ell$", loc="left",
                  fontsize=11, color=INK, pad=10)
    axA.set_ylim(0, 4)

    axB.axhline(0, color=GRID, lw=1.4, zorder=1)
    axB.annotate("$m_B > 0$   ($\\delta m < m$)", xy=(12, 0.45),
                 fontsize=9, color=MUTED)
    axB.annotate("$m_B < 0$   ($\\delta m > m$)", xy=(12, -1.6),
                 fontsize=9, color=MUTED)
    axB.set_ylabel("$m_B / m = 1 - \\delta m/m$")
    axB.set_title("B.  Signo del bare mass\n"
                  "separa las dos ramas de la Fig. 2", loc="left",
                  fontsize=11, color=INK, pad=10)
    axB.set_ylim(-2.5, 1.2)

    for ax in (axA, axB):
        ax.set_xlabel("$r_0/\\ell$   (cutoff $\\Lambda = 1/\\ell$ creciente $\\rightarrow$)")
        ax.set_xlim(0, 20)
        ax.legend(frameon=False, fontsize=9, loc="upper right")
        ax.grid(True, color=GRID, lw=.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)

    fig.suptitle("Impacto del cutoff $\\ell$ sobre la renormalización de masa "
                 "(Polonyi 2019, ec. 10)", x=0.011, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig("figures/masa_vs_cutoff.png", dpi=170, facecolor="#fcfcfb")
    print("figures/masa_vs_cutoff.png")

    print("\nTabla: delta_m/m y m_B/m")
    print(f"{'r0/ell':>7} | {'suavizado dm/m':>15} {'m_B/m':>9} | "
          f"{'desplazado dm/m':>16} {'m_B/m':>9}")
    for xi in [0.5, 1, 2, 3, 4, 6, 8, 12, 20]:
        a = make_regulator("smeared", 1.0 / xi).mass_shift_over_m()
        b = make_regulator("shifted", 1.0 / xi).mass_shift_over_m()
        print(f"{xi:7.1f} | {a:15.4f} {1-a:9.4f} | {b:16.4f} {1-b:9.4f}")


if __name__ == "__main__":
    main()
