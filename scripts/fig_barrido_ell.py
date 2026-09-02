"""Barrido en ell: como el cutoff afecta a la dinamica y a la masa desnuda.

No hay fisica nueva aqui. Es un bucle sobre `integrate_memory` (ec. 17) mas
las formas cerradas de `core`. La prediccion linealizada del bloque 1 se
superpone como envolvente, lo que convierte la figura en una validacion
cruzada en vez de una ilustracion.

Uso:  python3 scripts/fig_barrido_ell.py [--force]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid.core import Params, make_regulator
from nlaid.block1_linear import dominant_pole
from nlaid.block3_memory import integrate_memory
from nlaid.block4_renorm import linearized_rate

CACHE = "figures/barrido_ell.npz"
ELLS = np.array([2.0, 1.0, 0.5, 0.3, 0.2])      # r0/ell = 0.5, 1, 2, 3.33, 5
M_OVER_MB = 0.5
S_END, DS = 12.0, 5e-3

# Rampa ordinal de un solo tono (validada: L monotona, contraste 2.06:1).
RAMPA = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"


def compute(force=False):
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        if np.array_equal(d["ells"], ELLS) and float(d["m_over_mB"]) == M_OVER_MB:
            return {k: d[k] for k in d.files}

    out = {"ells": ELLS, "m_over_mB": M_OVER_MB}
    for ell in ELLS:
        wl = integrate_memory(Params(ell=ell, m_over_mB=M_OVER_MB),
                              s_end=S_END, ds=DS, n_ell=25.0)
        m = wl.s >= 0.0
        out[f"s_{ell}"] = wl.s[m]
        out[f"a_{ell}"] = np.linalg.norm(wl.a[m], axis=1)
        out[f"v_{ell}"] = wl.v[m][:, 1]                 # velocidad espacial
        out[f"tasa_{ell}"] = np.array([linearized_rate(Params(ell=ell))])
        print(f"  ell={ell:5.2f}  r0/ell={1/ell:5.2f}   "
              f"tasa linealizada={out[f'tasa_{ell}'][0]:+.4f}   "
              f"deriva={wl.norm_drift.max():.1e}", flush=True)
    np.savez(CACHE, **out)
    return out


def main():
    d = compute(force="--force" in sys.argv)
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(14, 4.3))

    for ell, c in zip(ELLS, RAMPA):
        s, a, v = d[f"s_{ell}"], d[f"a_{ell}"], d[f"v_{ell}"]
        lab = f"$r_0/\\ell$ = {1/ell:.2f}"
        axA.semilogy(s, np.maximum(a, 1e-14), lw=1.5, color=c, label=lab)
        axB.plot(s, v, lw=1.5, color=c, label=lab)

        # Envolvente predicha por el bloque 1 (espectral), independiente.
        tasa = d[f"tasa_{ell}"][0]
        if np.isfinite(tasa) and tasa < 0:
            ss = s[s > 1.0]
            i0 = np.searchsorted(s, 1.5)
            axA.semilogy(ss, a[i0] * np.exp(tasa * (ss - s[i0])),
                         lw=1, ls=(0, (4, 3)), color=c, alpha=.75)

    axA.set_xlabel("$s / r_0$"); axA.set_ylabel("$|\\ddot{x}|\\, r_0$")
    axA.set_title("A.  Aceleración: donde se ve $\\ell$\n"
                  "(punteado: envolvente del bloque 1, ec. 14)",
                  loc="left", fontsize=10.5, color=INK, pad=8)
    axA.set_ylim(1e-12, 5)

    axB.set_xlabel("$s / r_0$"); axB.set_ylabel("$\\dot{x}^1$")
    axB.set_title("B.  Velocidad: el pulso externo y su relajación",
                  loc="left", fontsize=10.5, color=INK, pad=8)

    # Panel C: masa desnuda, forma cerrada de la ec. (10). Sin calculo nuevo.
    x = np.linspace(0.05, 8, 400)
    for kind, ls, nom in (("smeared", "-", "suavizado ec. (5)"),
                          ("shifted", "--", "desplazado ec. (4)")):
        y = [1.0 - make_regulator(kind, 1/xi).mass_shift_over_m() for xi in x]
        axC.plot(x, y, lw=2, ls=ls, color=MUTED, label=nom)
    for ell, c in zip(ELLS, RAMPA):
        xi = 1/ell
        axC.plot([xi], [1 - make_regulator("smeared", ell).mass_shift_over_m()],
                 "o", ms=8, color=c, mec="#fcfcfb", mew=1.5, zorder=5)
    axC.axhline(0, color=GRID, lw=1.4)
    axC.annotate("$m_B > 0$", xy=(5.6, .45), fontsize=9, color=MUTED)
    axC.annotate("$m_B < 0$", xy=(5.6, -.55), fontsize=9, color=MUTED)
    axC.set_xlabel("$r_0/\\ell$"); axC.set_ylabel("$m_B/m = 1 - \\delta m/m$")
    axC.set_title("C.  Masa desnuda, ec. (10)\n(puntos: los $\\ell$ del barrido)",
                  loc="left", fontsize=10.5, color=INK, pad=8)
    axC.set_ylim(-1.2, 1.1); axC.set_xlim(0, 8)
    axC.legend(frameon=False, fontsize=8.5, loc="lower left")

    for ax in (axA, axB):
        ax.legend(frameon=False, fontsize=8.5)
    for ax in (axA, axB, axC):
        ax.grid(True, color=GRID, lw=.8, zorder=0); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)

    fig.suptitle(f"Barrido en el cutoff $\\ell$  —  ec. (17), $m/m_B$ = {M_OVER_MB}",
                 x=0.008, ha="left", fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig("figures/barrido_ell.png", dpi=170, facecolor="#fcfcfb")
    print("\nfigures/barrido_ell.png")


if __name__ == "__main__":
    main()
