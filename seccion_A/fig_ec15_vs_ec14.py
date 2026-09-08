"""Ec. (15), el limite local, contra la ec. (14) completa.

    chi^r_omega = 1 + r0 omega [ (2/3) i + O(omega ell) ]                  (15)

La (15) es la (14) con el cutoff removido a frecuencia fija: se queda con el
termino de Abraham-Lorentz y tira la memoria. Su unico cero esta en
omega = 3i/2 r0, en el semiplano SUPERIOR, para cualquier ell: el runaway.

Panel A  Re y Im de chi contra omega real a r0/ell = 4 (el cutoff critico).
         La (14) toca cero en omega = 1/ell -- el modo marginal --, la (15) no
         se anula nunca sobre el eje real.
Panel B  el corchete (chi - 1)/(r0 omega) contra omega*ell. Es funcion de
         omega*ell SOLAMENTE, asi que las tres curvas colapsan: eso vuelve
         literal el "O(omega ell)" de la ec. (15), cuyo limite es la recta 2i/3.
Panel C  Im del cero dominante contra r0/ell. La (15) esta clavada en 3/2,
         siempre inestable; la (14) lo mantiene abajo hasta r0/ell = 4.

Necesita el cache de `fig_ec14_polos.py` para el panel C.

Uso:  python3 seccion_A/fig_ec15_vs_ec14.py
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid import RAIZ
from nlaid.core import Params, make_regulator
from nlaid.block1_linear import susceptibility

CUTOFFS = [2.0, 3.0, 4.0]
CRITICO = 4.0
LOCUS = RAIZ / "figures/seccionA_ec14_locus.npz"

RAMPA = ["#86b6ef", "#2a78d6", "#104281"]
INK, MUTED, GRID, ROJO = "#1a1a19", "#5c5b54", "#e4e3dd", "#c1442e"


def chi(w, x):
    ell = 1.0 / x
    return susceptibility(w, make_regulator("smeared", ell), Params(ell=ell))


def panel_eje_real(ax):
    """chi(omega) sobre el eje real, a r0/ell = 4."""
    w = np.linspace(0.02, 9.0, 700)
    c = chi(w, CRITICO)
    loc = 1.0 + (2.0 / 3.0) * 1j * w

    ax.plot(w, c.real, lw=2, color=RAMPA[2], label="Re $\\chi$, ec. (14)")
    ax.plot(w, c.imag, lw=2, ls=(0, (5, 2)), color=RAMPA[0], label="Im $\\chi$, ec. (14)")
    ax.plot(w, loc.real, lw=1.6, color=ROJO, label="Re $\\chi$, ec. (15) $\\equiv 1$")
    ax.plot(w, loc.imag, lw=1.6, ls=(0, (5, 2)), color=ROJO,
            label="Im $\\chi$, ec. (15) $= \\frac{2}{3}\\omega$")

    ax.plot([CRITICO], [0.0], "o", ms=7, color=INK, mec="#fcfcfb", mew=1.4, zorder=6)
    ax.annotate("$\\chi = 0$ en $\\omega = 1/\\ell$:\nmodo marginal (la ec. 15\nno se anula jamás)",
                xy=(CRITICO, 0.0), xytext=(4.6, 1.9), fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=.9))
    ax.axhline(0, color=GRID, lw=1.2, zorder=0)
    ax.set_xlabel("$\\omega$   ($r_0 = 1$)")
    ax.set_ylabel("$\\chi^r_\\omega$")
    ax.set_title("A.  Sobre el eje real, a $r_0/\\ell = 4$\n"
                 "la ec. (15) sube sin límite; la (14) se dobla y toca cero",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, 9)
    ax.set_ylim(-1.2, 6.2)
    ax.legend(frameon=False, fontsize=8.5, loc="upper left")


def panel_corchete(ax):
    """(chi - 1)/(r0 omega) contra omega*ell: universal, y su limite es la ec. (15)."""
    # Grosor decreciente: las curvas coinciden a 1e-9, asi que superpuestas al
    # mismo ancho se verian como UNA sola y el colapso -- que es el mensaje --
    # no se leeria.
    for x, c, lw in zip(CUTOFFS, RAMPA, (6.0, 3.0, 1.4)):
        w = np.linspace(0.01, 3.0, 400) * x        # omega*ell de 0.01 a 3
        g = (chi(w, x) - 1.0) / w                  # r0 = 1
        ax.plot(w / x, g.imag, lw=lw, color=c, label=f"$r_0/\\ell$ = {x:.0f}")
        ax.plot(w / x, g.real, lw=lw, ls=(0, (4, 2)), color=c)

    ax.axhline(2.0 / 3.0, color=ROJO, lw=1.8)
    ax.annotate("ec. (15): $2i/3$, sin parte real", xy=(1.05, 0.72),
                fontsize=9, color=ROJO)
    ax.axhline(0, color=GRID, lw=1.2, zorder=0)
    ax.axvline(1.0, color=MUTED, lw=1, ls=(0, (2, 3)))
    ax.annotate("$\\omega\\ell = 1$", xy=(1.03, -0.62), fontsize=8.5, color=MUTED)
    ax.annotate("continua: Im     discontinua: Re", xy=(1.55, -0.30),
                fontsize=8.5, color=MUTED)
    ax.set_xlabel("$\\omega\\,\\ell$")
    ax.set_ylabel("$(\\chi^r_\\omega - 1)\\,/\\,r_0\\omega$")
    ax.set_title("B.  El corchete de la ec. (15) es función de $\\omega\\ell$ sola\n"
                 "las tres curvas colapsan: eso es el $O(\\omega\\ell)$",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, 3)
    ax.set_ylim(-0.75, 0.95)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")


def panel_ceros(ax):
    """Im del cero dominante: la (14) contra el 3/2 fijo de la (15)."""
    if not os.path.exists(LOCUS):
        sys.exit(f"falta {LOCUS}: corra antes seccion_A/fig_ec14_polos.py")
    d = np.load(LOCUS)
    xs, z = d["xs"], d["z"]
    m = np.isfinite(z)

    ax.axhspan(0, 8, color=ROJO, alpha=.07, zorder=0)
    ax.plot(xs[m], z[m].imag, lw=2.2, color=RAMPA[2],
            label="ec. (14): cero dominante")
    ax.axhline(1.5, lw=2, color=ROJO, label="ec. (15): $\\omega = 3i/2$, fijo")
    ax.axhline(0, color=INK, lw=1.2, zorder=1)
    ax.axvline(CRITICO, color=MUTED, lw=1, ls=(0, (2, 3)))

    i = int(np.argmin(np.abs(z[m].imag)))
    ax.plot([xs[m][i]], [0.0], "o", ms=7, color=INK, mec="#fcfcfb", mew=1.4, zorder=6)
    ax.annotate(f"cruce en $r_0/\\ell$ = {xs[m][i]:.2f}\n"
                "(estable a la izquierda: el\ncutoff retiene el polo abajo)",
                xy=(xs[m][i], 0), xytext=(5.2, -1.35), fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=.9))
    ax.annotate("inestable", xy=(15.5, 0.55), fontsize=8.5, color=ROJO)

    ax.set_xlabel("$r_0/\\ell$   (cutoff $\\Lambda = 1/\\ell$ creciente $\\rightarrow$)")
    ax.set_ylabel("Im $\\omega$ del cero dominante")
    ax.set_title("C.  La ec. (15) es inestable para todo $\\ell$\n"
                 "la (14) solo por encima de $r_0/\\ell = 4$",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(2, 20)
    ax.set_ylim(-2.0, 7.6)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def verifica():
    """El colapso del panel B, y el cero de la (15) en 3i/2."""
    print("\ncorchete (chi-1)/(r0 omega) a omega*ell = 0.5, por cutoff:")
    for x in CUTOFFS:
        g = (chi(0.5 * x, x) - 1.0) / (0.5 * x)
        print(f"  r0/ell={x:4.1f}   {g.real:+.9f}{g.imag:+.9f}i")
    w0 = 1.5j
    print(f"cero de la ec. (15):  1 + (2/3) i r0 (3i/2) = {1 + (2/3)*1j*w0:.1e}")


def main():
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.7))
    panel_eje_real(axes[0])
    panel_corchete(axes[1])
    panel_ceros(axes[2])

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

    fig.suptitle("Ec. (15) frente a la ec. (14): el límite local se queda con "
                 "Abraham-Lorentz y tira la memoria", x=0.008, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salida = RAIZ / "figures/seccionA_ec15_vs_ec14.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)
    verifica()


if __name__ == "__main__":
    main()
