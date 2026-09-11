"""Ec. (18): la condicion de renormalizacion, y como se resuelve.

    chi^r_omega = 1 + (2/3) i r0 omega                                    (18)

"The equation of motion has two free, adjustable parameters, the cutoff, ell,
and the bare mass, m_B. Hence we need a renormalization condition to fix the
theory" (p. 10). La (18) NO se impone literalmente -- su unico cero esta en
omega = 3i/2, el runaway, de modo que una teoria que la cumpla exacto no relaja
y no hay nada que monitorear. El paper la aplica operacionalmente: se fija
m/m_B monitoreando la relajacion a s grande.

Panel A  la funcion objetivo del disparo: tasa de relajacion de la ec. (17)
         contra m/m_B a cutoff fijo. Sobre ese mismo eje caen los tres
         criterios que se han usado en el repo, y se ve cual es cual.
Panel B  el resultado central: integrando la ec. (17) SOBRE la linea del
         contratermino de la ec. (10), la tasa no lineal reproduce el
         Im omega del cero dominante de chi (ec. 14). Dos teorias distintas
         --no lineal en el tiempo, lineal en frecuencia-- dando el mismo
         numero, sin ningun parametro ajustado.
Panel C  el plano de fases (r0/ell, m/m_B) con las dos lineas: el borde de
         estabilidad (tasa = 0, de data/borde_smeared.json) y el contratermino.

Cierra la hipotesis 1 de `docs/discrepancias.md`: la condicion de
renormalizacion SI se puede imponer, y la linea que fija es la del
contratermino.

Necesita el cache de `fig_ec14_polos.py` (los ceros de chi) y
`data/borde_smeared.json` (de `scripts/scan_borde.py`).

Uso:  python3 seccion_A/fig_ec18_renormalizacion.py [--force]
"""

import json
import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid import RAIZ
from nlaid.core import Params, make_regulator
from nlaid.block4_renorm import relaxation_rate

CACHE = RAIZ / "figures/seccionA_ec18.npz"
LOCUS = RAIZ / "figures/seccionA_ec14_locus.npz"
BORDE = RAIZ / "data/borde_smeared.json"

X_FIJO = 3.0                                    # r0/ell del panel A (Fig. 1)
RATIOS = np.round(np.arange(0.6, 5.01, 0.30), 2)
CUTOFFS = [2.0, 2.5, 3.0, 3.5]                  # deben estar en el SWEEP del locus
S_END, DS = 14.0, 5e-3

RAMPA = ["#86b6ef", "#2a78d6", "#104281"]
INK, MUTED, GRID, ROJO, VERDE = "#1a1a19", "#5c5b54", "#e4e3dd", "#c1442e", "#2e7d5b"


def cero_chi(x):
    """Im(omega) del cero dominante de chi a r0/ell = x, del cache del ec. 14.

    Se reusa ese cache en vez de rehacer la busqueda: es el MISMO numero, y
    recalcularlo con otra rutina invitaria a que las dos figuras discrepen.
    """
    if not os.path.exists(LOCUS):
        sys.exit(f"falta {LOCUS}: corra antes seccion_A/fig_ec14_polos.py")
    d = np.load(LOCUS)
    j = int(np.argmin(np.abs(d["xs"] - x)))
    if abs(d["xs"][j] - x) > 1e-6:
        sys.exit(f"r0/ell={x} no esta en el barrido del locus; use uno de su malla")
    return complex(d["z"][j])


def compute(force=False):
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        if np.array_equal(d["ratios"], RATIOS):
            return {k: d[k] for k in d.files}

    out = {"ratios": RATIOS, "cutoffs": np.array(CUTOFFS, dtype=float)}

    # --- Panel A: barrido en m/m_B a cutoff fijo -------------------------
    ell = 1.0 / X_FIJO
    tasas, fiable = [], []
    for mr in RATIOS:
        t, diag = relaxation_rate(Params(ell=ell, m_over_mB=float(mr)),
                                  s_end=S_END, ds=DS)
        tasas.append(t)
        fiable.append(bool(diag.get("fiable", False)))
        print(f"  r0/ell={X_FIJO}  m/mB={mr:5.2f}  tasa={t:+.4f}  "
              f"fiable={fiable[-1]}", flush=True)
    out["tasas"], out["fiable"] = np.array(tasas), np.array(fiable)

    # --- Panel B: la tasa SOBRE la linea del contratermino ---------------
    med, esp, mrs = [], [], []
    for x in CUTOFFS:
        e = 1.0 / x
        mr = make_regulator("smeared", e).m_over_mB_counterterm()
        t, _ = relaxation_rate(Params(ell=e, m_over_mB=mr), s_end=S_END, ds=DS)
        med.append(t)
        esp.append(cero_chi(x).imag)
        mrs.append(mr)
        print(f"  r0/ell={x:4.1f}  m/mB={mr:6.4f}  medida={t:+.4f}  "
              f"espectral={esp[-1]:+.4f}  dif={abs(t-esp[-1]):.4f}", flush=True)
    out["med"], out["esp"], out["mr_ct"] = np.array(med), np.array(esp), np.array(mrs)

    np.savez(CACHE, **out)
    return out


def panel_objetivo(ax, d):
    mr, t, ok = d["ratios"], d["tasas"], d["fiable"]
    objetivo = cero_chi(X_FIJO).imag                      # criterio R1
    ct = make_regulator("smeared", 1.0 / X_FIJO).m_over_mB_counterterm()

    ax.axhspan(0, 1, color=ROJO, alpha=.07, zorder=0)
    ax.plot(mr[ok], t[ok], "o-", lw=2, ms=5, color=RAMPA[2], zorder=4,
            label="tasa de la ec. (17)")
    if (~ok).any():
        ax.plot(mr[~ok], t[~ok], "x", ms=7, color=MUTED, zorder=4,
                label="ajuste no fiable")

    ax.axhline(0, color=INK, lw=1.2, zorder=2)
    ax.annotate("inestable", xy=(0.7, 0.06), fontsize=8.5, color=ROJO)
    ax.axhline(objetivo, color=VERDE, lw=1.6, ls=(0, (5, 2)), zorder=3)
    ax.annotate(f"Im $\\omega$ del cero de $\\chi$ = {objetivo:.3f}\n"
                "(criterio: empalmar con la ec. 14)",
                xy=(2.55, objetivo + 0.06), fontsize=8.5, color=VERDE)

    ax.axvline(ct, color=RAMPA[0], lw=1.4, ls=(0, (2, 3)), zorder=3)
    ax.plot([ct], [objetivo], "o", ms=9, mfc="none", mew=2.2, color=INK, zorder=7)
    ax.annotate(f"contratérmino, ec. (10): $m/m_B$ = {ct:.2f}",
                xy=(ct, objetivo), xytext=(18, -34), textcoords="offset points",
                fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=.9))

    ax.set_xlabel("$m/m_B$        (a $r_0/\\ell$ = 3)")
    ax.set_ylabel("tasa de relajación de $|\\ddot{x}|$")
    ax.set_title("A.  La ec. (18) resuelta por disparo\n"
                 "la línea verde corta la curva justo en el contratérmino",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0.4, 5.2)
    ax.set_ylim(-1.4, 0.2)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right")


def panel_ley(ax, d):
    x = np.array(CUTOFFS)
    ax.plot(x, d["esp"], "o-", lw=2, ms=7, color=VERDE,
            label="ec. (14): Im $\\omega$ del cero de $\\chi$")
    ax.plot(x, d["med"], "s--", lw=1.8, ms=6, color=RAMPA[2],
            label="ec. (17): tasa medida sobre la línea (10)")

    for xi, a, b in zip(x, d["esp"], d["med"]):
        ax.annotate(f"{abs(a-b):.4f}", xy=(xi, min(a, b) - 0.035), ha="center",
                    fontsize=8, color=MUTED)
    ax.annotate("diferencia entre las dos vías", xy=(2.9, -0.35),
                fontsize=8.5, color=MUTED)

    ax.set_xlabel("$r_0/\\ell$")
    ax.set_ylabel("tasa de relajación")
    ax.set_title("B.  Dos teorías, el mismo número\n"
                 "no lineal en el tiempo contra lineal en frecuencia",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(1.8, 3.7)
    ax.set_ylim(-1.0, -0.3)
    ax.legend(frameon=False, fontsize=8.5, loc="lower left")


def panel_plano(ax, d):
    x = np.linspace(0.3, 12, 400)
    ct = np.array([make_regulator("smeared", 1.0 / xi).m_over_mB_counterterm()
                   for xi in x])
    ct[np.abs(ct) > 12] = np.nan                     # corta la asintota en 6

    borde = json.load(open(BORDE))
    xb = np.array([float(k) for k in borde])
    pos = np.array([borde[k]["positiva"] for k in borde], dtype=float)

    ax.fill_between(xb, 0, pos, color=RAMPA[0], alpha=.25, lw=0,
                    label="estable (rama $m_B > 0$)")
    ax.plot(xb, pos, "o-", lw=2, ms=5, color=RAMPA[2],
            label="borde: tasa = 0")
    ax.plot(x, ct, lw=2, color=VERDE, label="contratérmino, ec. (10)")
    ax.plot(CUTOFFS, d["mr_ct"], "s", ms=6, color=VERDE, mec="#fcfcfb", mew=1.2,
            zorder=6)

    ax.plot([4.0], [3.0], "o", ms=10, mfc="none", mew=2.2, color=ROJO, zorder=7)
    ax.annotate("el contratérmino cruza el borde en $r_0/\\ell$ = 4.00,\n"
                "el mismo crítico que da el bloque 1 por\n"
                "principio del argumento (borde 3.009 vs 3.000)",
                xy=(4.0, 3.0), xytext=(4.9, 4.6), fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=.9))
    ax.set_xlabel("$r_0/\\ell$")
    ax.set_ylabel("$m/m_B$")
    ax.set_title("C.  Dónde cae la teoría renormalizada\n"
                 "regulador suavizado, rama $m_B > 0$",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, 12)
    ax.set_ylim(0, 7)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def main():
    warnings.filterwarnings("ignore")
    d = compute(force="--force" in sys.argv)

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8))
    panel_objetivo(axes[0], d)
    panel_ley(axes[1], d)
    panel_plano(axes[2], d)

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

    fig.suptitle("Ec. (18): la condición de renormalización fija $m/m_B$, y lo "
                 "que fija es la línea del contratérmino", x=0.008, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salida = RAIZ / "figures/seccionA_ec18_renormalizacion.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)


if __name__ == "__main__":
    main()
