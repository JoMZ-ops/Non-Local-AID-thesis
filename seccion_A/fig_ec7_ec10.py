"""Seccion II.A del paper, por dentro: el nucleo de la ec. (7) y el de la ec. (10).

    xddot = 4 r_0B int_{-inf}^{0} du delta_B'(u^2) (x - x' + u xdot')      (7)
    delta_m = (e^2 / 2 c^2) int_0^inf dz z^{-1/2} delta_B(z)              (10)

Ninguna de las dos es una curva: son integrales. Lo que se grafica es su
INTEGRANDO, que es donde se ve el mecanismo.

Panel A  nucleo 4 r0 delta_B'(u^2) de la ec. (7) contra u/r0, para tres cutoffs.
         Cambia de signo en u = -2 ell y decae como e^{-|u|/ell}: la memoria es
         infinita pero pesa solo a escala del cutoff, y al remover el cutoff se
         angosta como ell y crece como ell^-4.
Panel B  integrando z^{-1/2} delta_B(z) de la ec. (10), adimensionalizado. Su
         area es 1/(3 ell), de donde delta_m/m = r0/(6 ell).

Uso:  python3 seccion_A/fig_ec7_ec10.py
"""

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import quad

from nlaid import RAIZ
from nlaid.core import make_regulator

CUTOFFS = [2.0, 3.0, 4.0]                        # r0/ell
RAMPA = ["#86b6ef", "#2a78d6", "#104281"]        # magnitud, no identidad
INK, MUTED, GRID, NARANJA = "#1a1a19", "#5c5b54", "#e4e3dd", "#eb6834"


def panel_ec7(ax):
    """Nucleo 4 r0 delta_B'(u^2) contra u/r0 (r0 = 1)."""
    u = np.linspace(-3.0, 0.0, 1200)
    for x, c in zip(CUTOFFS, RAMPA):
        ell = 1.0 / x
        k = 4.0 * make_regulator("smeared", ell).d_delta(u ** 2)
        ax.plot(u, k, lw=2, color=c, label=f"$r_0/\\ell$ = {x:.0f}")
        ax.plot([-2 * ell], [0.0], "o", ms=5, color=c, mec="#fcfcfb", mew=1.2, zorder=5)

    ax.axhline(0, color=GRID, lw=1.4, zorder=0)
    ax.axvline(-1.0 / CUTOFFS[1], color=NARANJA, lw=1.4, ls=(0, (4, 3)))
    ax.annotate("ec. (4): toda la memoria\nconcentrada en $u = -\\ell$",
                xy=(-0.42, 4.5), fontsize=8.5, color=NARANJA, ha="right")
    ax.annotate("cambio de signo en $u = -2\\ell$;\ncola negativa, $\\sim$2% del pico",
                xy=(-0.5, 0), xytext=(-1.95, -3.2), fontsize=8.5, color=MUTED,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=.8))
    ax.set_xlabel("$u / r_0$        ($u = s' - s \\leq 0$)")
    ax.set_ylabel("$4\\,r_0\\,\\delta_B'(u^2)$   (escala symlog)", color=INK)
    ax.set_title("A.  Núcleo de la ec. (7)\nal remover el cutoff se angosta "
                 "como $\\ell$ y crece como $\\ell^{-4}$",
                 loc="left", fontsize=11, color=INK, pad=10)
    # Symlog: el nucleo cambia de signo y sus dos lobulos difieren en un factor
    # 40. En lineal el lobulo negativo -- el que apaga la autointeraccion mas
    # alla de 2 ell -- se ve como una raya sobre el cero.
    ax.set_yscale("symlog", linthresh=1.0, linscale=0.6)
    ax.set_xlim(-2.05, 0.05)
    ax.set_ylim(-6, 200)
    ax.legend(frameon=False, fontsize=9, loc="upper left")


def panel_ec10(ax):
    """Integrando de la ec. (10). Adimensional: ell * z^{-1/2} delta_B(z) vs z/ell^2.

    La forma es universal en z/ell^2 -- el regulador suavizado no tiene otra
    escala --, asi que una sola curva sirve para todo cutoff. Lo que cambia con
    ell es solo el factor global 1/ell, que es la divergencia de la ec. (10).
    """
    reg = make_regulator("smeared", 1.0)         # ell = 1: z/ell^2 = z
    z = np.linspace(1e-12, 60.0, 2000)
    g = reg.delta(z) / np.sqrt(z)                # = sqrt(z) e^{-sqrt z} / 12
    acum = np.concatenate([[0.0], np.cumsum(np.diff(z) * (g[1:] + g[:-1]) / 2)])

    ax.fill_between(z, 0, g, color=RAMPA[0], alpha=.35, lw=0)
    ax.plot(z, g, lw=2, color=RAMPA[2])
    ax.plot([1.0], [float(g.max())], "o", ms=5, color=RAMPA[2],
            mec="#fcfcfb", mew=1.2, zorder=5)
    ax.annotate("máximo en $z = \\ell^2$", xy=(1.0, float(g.max())),
                xytext=(6, 0.0335), fontsize=8.5, color=MUTED,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=.8))
    ax.annotate("área $= 1/3$\n$\\Rightarrow\\ \\delta m/m = r_0/6\\ell$",
                xy=(31, 0.0225), fontsize=10, color=INK)

    # Eje derecho: que fraccion del contratermino se ha acumulado hasta z.
    ax2 = ax.twinx()
    ax2.plot(z, acum / (1.0 / 3.0), lw=1.8, ls=(0, (5, 2)), color=NARANJA)
    ax2.set_ylabel("fracción acumulada  $\\int_0^Z / \\int_0^\\infty$",
                   color=NARANJA)
    ax2.set_ylim(0, 1.05)
    ax2.tick_params(colors=NARANJA, labelsize=9)
    for sp in ("top", "left"):
        ax2.spines[sp].set_visible(False)
    ax2.spines["right"].set_color(NARANJA)
    ax2.spines["bottom"].set_color(GRID)

    ax.set_xlabel("$z / \\ell^2$        ($z = x^2$, intervalo invariante)")
    ax.set_ylabel("$\\ell\\,z^{-1/2}\\delta_B(z)$", color=RAMPA[2])
    ax.set_title("B.  Integrando de la ec. (10)\n98% del contratérmino viene de "
                 "$z < 60\\,\\ell^2$: es de corto alcance",
                 loc="left", fontsize=11, color=INK, pad=10)
    ax.set_xlim(0, 60)
    ax.set_ylim(0, 0.038)


def verifica():
    """Tres vias al mismo numero: cuadratura, forma cerrada y el modulo."""
    ell = 0.7
    reg = make_regulator("smeared", ell)
    num = quad(lambda z: reg.delta(z) / np.sqrt(z), 0, np.inf, limit=400)[0]
    print(f"\nec. (10) con ec. (5), ell = {ell}")
    print(f"  cuadratura            {num:.12f}")
    print(f"  cerrada  1/(3 ell)    {1 / (3 * ell):.12f}")
    print(f"  reg.moment_inv_sqrt   {reg.moment_inv_sqrt():.12f}")
    print(f"  delta_m/m = r0/(6ell) {reg.mass_shift_over_m():.12f}"
          f"   (esperado {1 / (6 * ell):.12f})")
    # La lectura erronea del extractor de texto, dz*sqrt(z), da 4*ell: CRECE al
    # remover el cutoff. Es el argumento dimensional de que la (10) lleva 1/sqrt(z).
    mal = quad(lambda z: reg.delta(z) * np.sqrt(z), 0, np.inf, limit=400)[0]
    print(f"  lectura erronea dz*sqrt(z) = {mal:.6f}  (= 4 ell, diverge al reves)")


def main():
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.6))
    panel_ec7(axA)
    panel_ec10(axB)

    for ax in (axA, axB):
        ax.grid(True, color=GRID, lw=.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK)   # el color de los ejes y va por panel:
                                        # en B codifica cual curva es cual.

    fig.suptitle("Sección II.A: los integrandos de las ecs. (7) y (10) "
                 "(regulador suavizado, ec. 5)", x=0.011, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salida = RAIZ / "figures/seccionA_ec7_ec10.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)
    verifica()


if __name__ == "__main__":
    main()
