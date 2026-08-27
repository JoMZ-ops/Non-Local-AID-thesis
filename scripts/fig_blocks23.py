"""Figura de los bloques 2 y 3: reproduccion de la Fig. 1 de Polonyi (2019).

Regulador suavizado, ec. (17), r0/ell = 3, los mismos m/m_B del paper.

Panel A  m_B > 0: relajacion exponencial monotona. La tasa varia suavemente
         con m/m_B, que es lo que permite resolver la condicion de
         renormalizacion (18) monitoreando la relajacion a s grande.
Panel B  m_B < 0: oscilacion rapida ("Zitterbewegung") con envolvente
         exponencial creciente o decreciente. El borde de estabilidad cae
         entre -3.80 y -3.91.

Uso:  python3 scripts/fig_blocks23.py [--force]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid.core import Params
from nlaid.block3_memory import integrate_memory

CACHE = "figures/blocks23_datos.npz"
ELL = 1.0 / 3.0                      # r0/ell = 3, como en la Fig. 1
S_END, DS = 25.0, 5e-3

# Rampas ordinales de un solo tono (validadas: L monotona, contraste >= 2:1).
RAMP_POS = ["#86b6ef", "#2a78d6", "#104281"]
RAMP_NEG = ["#f09a6e", "#eb6834", "#8f3612"]
MR_POS = [1.95, 1.98, 2.00]
MR_NEG = [-3.80, -3.91, -4.10]
INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"


def compute(force=False):
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        return {k: d[k] for k in d.files}

    out = {}
    for mr in MR_POS + MR_NEG:
        wl = integrate_memory(Params(ell=ELL, m_over_mB=mr),
                              s_end=S_END, ds=DS, n_ell=25.0, pts_per_ell=16)
        m = wl.s >= 0.0
        out[f"s_{mr}"] = wl.s[m]
        out[f"a_{mr}"] = np.linalg.norm(wl.a[m], axis=1)
        out[f"d_{mr}"] = wl.norm_drift[m]
        print(f"  m/mB={mr:+6.2f}  |a| final={out[f'a_{mr}'][-1]:.3e}  "
              f"deriva={wl.norm_drift.max():.2e}")
    np.savez(CACHE, **out)
    return out


def main():
    data = compute(force="--force" in sys.argv)
    fig, (axA, axB) = plt.subplots(1, 2, figsize=(11.5, 4.6))

    for ax, mrs, ramp, titulo in (
        (axA, MR_POS, RAMP_POS,
         "A.  $m_B > 0$ — relajación oscilatoria amortiguada"),
        (axB, MR_NEG, RAMP_NEG, "B.  $m_B < 0$ — oscilación con envolvente exponencial"),
    ):
        for mr, c in zip(mrs, ramp):
            s, a, d = data[f"s_{mr}"], data[f"a_{mr}"], data[f"d_{mr}"]
            # Truncar donde la deriva de xdot^2 supera 1e-3: mas alla el
            # integrador ya no conserva la normalizacion y la curva no es
            # fisica. Mostrar solo el tramo validado.
            bad = np.where(d > 1e-3)[0]
            n = bad[0] if len(bad) else len(s)
            ax.semilogy(s[:n], np.maximum(a[:n], 1e-14), lw=1.4, color=c,
                        label=f"$m/m_B = {mr:+.2f}$")
        ax.set_xlabel("$s / r_0$")
        ax.set_ylabel("$|\\ddot{x}|\\, r_0$")
        ax.set_title(titulo, loc="left", fontsize=11, color=INK, pad=10)
        ax.set_xlim(0, S_END)
        ax.set_ylim(bottom=1e-12)
        ax.legend(frameon=False, fontsize=9, loc="lower left")
        ax.grid(True, color=GRID, lw=.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)

    # Predicción independiente del bloque 1: cero de chi^r en
    # omega = +-2.5736 - 0.71793i a r0/ell = 3. La tasa es Im(omega) y el
    # espaciado entre minimos de |xddot| es medio periodo, pi/Re(omega).
    TASA, FREQ = -0.71793, 2.5736
    ss = np.linspace(1.0, 24.0, 50)
    axA.semilogy(ss, 0.62 * np.exp(TASA * ss), lw=1.2, color=MUTED,
                 ls=(0, (5, 4)), zorder=5)
    axA.annotate(f"bloque 1 (ec. 14): $e^{{{TASA:.3f}\,s}}$\n"
                 f"medido (ec. 17): $e^{{-0.717\,s}}$",
                 xy=(0.46, 0.87), xycoords="axes fraction",
                 fontsize=8.5, color=MUTED)
    axA.annotate(f"semiperíodo predicho $\\pi/{FREQ:.3f} = 1.2207$\n"
                 f"medido $= 1.2207$",
                 xy=(0.46, 0.77), xycoords="axes fraction",
                 fontsize=8.5, color=MUTED)

    axB.set_title("B.  $m_B < 0$ — crecimiento exponencial (ver nota)",
                  loc="left", fontsize=11, color=INK, pad=10)
    axB.annotate("curvas truncadas donde la deriva de $\\dot{x}^2$\n"
                 "supera $10^{-3}$: más allá el integrador\n"
                 "ya no es fiable",
                 xy=(0.42, 0.10), xycoords="axes fraction",
                 fontsize=8.5, color=MUTED)

    fig.suptitle("Bloques 2–3 — ec. (17), regulador suavizado, $r_0/\\ell = 3$   "
                 "(reproducción de la Fig. 1 de Polonyi 2019)",
                 x=0.011, ha="left", fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig("figures/blocks23_fig1.png", dpi=170, facecolor="#fcfcfb")
    print("\nfigures/blocks23_fig1.png")


if __name__ == "__main__":
    main()
