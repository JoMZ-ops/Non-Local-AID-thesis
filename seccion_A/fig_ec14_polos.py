"""Ec. (14) en el plano complejo: donde estan los polos y que deciden.

    chi^r_omega = 1 + r0 [ (2/3) i omega
                  - (2/omega^2) int_{-inf}^{0} du delta_B(u^2) N(omega u)/u^2 ]   (14)

La respuesta lineal es F^r_omega = 1/[(omega + i eps)^2 chi^r_omega], asi que los
POLOS de la respuesta son los CEROS de chi. Un cero con Im omega > 0 es un modo
e^{-i omega s} que crece: inestabilidad. El criterio del paper (p. 8) es que chi
sea analitica y sin ceros en el semiplano superior.

Panel A  |chi| sobre el plano complejo a r0/ell = 3: estable, todos los ceros
         por debajo del eje real.
Panel B  lo mismo a r0/ell = 5: el par dominante ya cruzo. Inestable.
Panel C  lugar de raices del par dominante al barrer r0/ell. Cruza el eje real
         en r0/ell = 4.00 (con omega ell = 1), los dos ceros colisionan sobre el
         eje imaginario cerca de r0/ell = 5.8 y la rama desciende hacia 3i/2:
         el polo runaway de Abraham-Lorentz de la ec. (15), que el cutoff
         mantiene a raya mientras r0/ell < 4.

Uso:  python3 seccion_A/fig_ec14_polos.py [--force]
"""

import os
import sys
import warnings

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.optimize import newton

from nlaid import RAIZ
from nlaid.core import Params, make_regulator
from nlaid.block1_linear import convergence_lower_bound, find_zeros_uhp, susceptibility

MAPAS = [(3.0, "A", "estable"), (5.0, "B", "inestable")]
SWEEP = np.concatenate([np.linspace(2.0, 5.4, 35), np.linspace(5.6, 20.0, 30)])
CACHE = "figures/seccionA_ec14_locus.npz"

INK, MUTED, GRID = "#1a1a19", "#5c5b54", "#e4e3dd"
AZUL, ROJO, VERDE = "#2a78d6", "#c1442e", "#2e7d5b"


def ceros(reg, pr, re_max, im_lo, im_hi, n=(200, 150)):
    """Ceros de chi en la caja. Minimos locales de |chi| + Newton complejo.

    # ponytail: duplica el patron de block1_linear.find_zeros_uhp, que por
    # contrato solo mira Im omega > 0. Aqui hace falta el semiplano INFERIOR
    # (es donde viven los ceros del caso estable, que es justo lo que hay que
    # mostrar). La parte superior se contrasta contra esa funcion en `verifica`.
    """
    R, I = np.meshgrid(np.linspace(-re_max, re_max, n[0]),
                       np.linspace(im_lo, im_hi, n[1]), indexing="ij")
    W = R + 1j * I
    A = np.abs(campo(W, reg, pr))
    f = lambda w: susceptibility(np.asarray(w), reg, pr)
    out = []
    for i in range(1, A.shape[0] - 1):
        for j in range(1, A.shape[1] - 1):
            if A[i, j] > A[i - 1:i + 2, j - 1:j + 2].min():
                continue
            try:
                r = newton(f, W[i, j], tol=1e-13, maxiter=80)
            except Exception:                                     # noqa: BLE001
                continue
            if (abs(f(r)) < 1e-8 and im_lo < r.imag < im_hi
                    and abs(r.real) < re_max
                    and not any(abs(r - z) < 1e-5 for z in out)):
                out.append(complex(r))
    return sorted(out, key=lambda z: -z.imag)


def campo(W, reg, pr):
    """chi sobre una malla 2D, por columnas: la version plana pide ~100 MB."""
    return np.stack([susceptibility(col, reg, pr) for col in W])


def locus(force=False):
    """Cero dominante (Re >= 0) contra r0/ell. Se cachea: domina el costo."""
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        if np.array_equal(d["xs"], SWEEP):
            return d["z"]
    z = []
    for x in SWEEP:
        ell = 1.0 / x
        reg, pr = make_regulator("smeared", ell), Params(ell=ell)
        # La representacion integral converge para Im omega > -1/ell; se deja
        # un margen para no evaluar sobre el borde mismo.
        zs = [w for w in ceros(reg, pr, 12.0, 0.85 * convergence_lower_bound(reg), 9.0)
              if w.real >= -1e-9]
        z.append(zs[0] if zs else np.nan + 0j)
        print(f"  r0/ell={x:5.2f}  dominante = {z[-1]}", flush=True)
    z = np.array(z)
    np.savez(CACHE, xs=SWEEP, z=z)
    return z


VENTANA = (8.0, -2.5, 5.0)          # (re_max, im_lo, im_hi), la MISMA en A y B:
                                    # comparar dos mapas con ejes distintos no
                                    # dice nada sobre quien cruzo el eje.


def panel_mapa(ax, x, letra, veredicto):
    ell = 1.0 / x
    reg, pr = make_regulator("smeared", ell), Params(ell=ell)
    re_max, im_lo, im_hi = VENTANA
    assert im_lo > convergence_lower_bound(reg)     # la integral debe converger

    R, I = np.meshgrid(np.linspace(-re_max, re_max, 320),
                       np.linspace(im_lo, im_hi, 260), indexing="ij")
    A = np.log10(np.abs(campo(R + 1j * I, reg, pr)) + 1e-12)

    ax.contourf(R, I, A, levels=np.linspace(-2.0, 1.6, 22), cmap="Blues_r")
    ax.contour(R, I, A, levels=[-1.0, -0.5, 0.0], colors="#ffffff",
               linewidths=.5, alpha=.5)
    ax.axhspan(0, im_hi, color=ROJO, alpha=.07, zorder=2)
    ax.axhline(0, color=INK, lw=1.2, zorder=3)
    ax.annotate("Im $\\omega > 0$: modos que crecen", xy=(-re_max + 0.4, im_hi - 0.65),
                fontsize=8.5, color=ROJO, zorder=6,
                bbox=dict(fc="#fcfcfb", ec="none", alpha=.8, pad=2))
    if letra == "A":
        ax.annotate("azul oscuro: $|\\chi^r_\\omega|$ pequeño,\ncerca de un cero",
                    xy=(-re_max + 0.4, im_lo + 0.35), fontsize=8, color=MUTED,
                    zorder=6, bbox=dict(fc="#fcfcfb", ec="none", alpha=.8, pad=2))

    # Se descartan los ceros pegados al marco: un minimo local sobre una malla
    # recortada por el borde no es un cero, es el borde.
    dibujados = set()
    for w in ceros(reg, pr, re_max, im_lo, im_hi):
        if abs(w.real) > re_max - 0.5 or w.imag < im_lo + 0.3:
            continue
        arriba = w.imag > 0
        dibujados.add(arriba)
        ax.plot([w.real], [w.imag], "x", ms=9, mew=2.2,
                color=ROJO if arriba else VERDE, zorder=7)
    for arriba, color, etiqueta in ((False, VERDE, "ceros con Im $\\omega < 0$"),
                                    (True, ROJO, "ceros con Im $\\omega > 0$")):
        if arriba in dibujados:
            ax.plot([], [], "x", ms=8, mew=2, color=color, label=etiqueta)

    ax.set_xlabel("Re $\\omega$   ($r_0 = 1$)")
    ax.set_ylabel("Im $\\omega$")
    ax.set_title(f"{letra}.  $|\\chi^r_\\omega|$ (log) a $r_0/\\ell$ = {x:.0f} — "
                 f"{veredicto}\n"
                 + ("ningún cero cruzó el eje real" if veredicto == "estable"
                    else "el par dominante está arriba: runaway"),
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(-re_max, re_max)
    ax.set_ylim(im_lo, im_hi)
    ax.legend(frameon=False, fontsize=8.5, loc="lower right", labelcolor=INK)


def panel_locus(ax, z):
    m = np.isfinite(z)
    xs, zz = SWEEP[m], z[m]
    sc = ax.scatter(zz.real, zz.imag, c=xs, cmap="viridis", s=26, zorder=5)
    ax.scatter(-zz.real, zz.imag, c=xs, cmap="viridis", s=26, alpha=.35, zorder=4)
    cb = ax.figure.colorbar(sc, ax=ax, pad=0.02)
    cb.set_label("$r_0/\\ell$", color=INK)
    cb.ax.tick_params(colors=MUTED, labelsize=8)

    ax.axhspan(0, 9, color=ROJO, alpha=.07, zorder=0)
    ax.axhline(0, color=INK, lw=1.2, zorder=1)
    ax.axhline(1.5, color=MUTED, lw=1, ls=(0, (5, 4)), zorder=1)
    ax.annotate("$3i/2$: polo runaway de\nAbraham-Lorentz, ec. (15)",
                xy=(-8.6, 1.75), fontsize=8.5, color=MUTED)

    i = int(np.argmin(np.abs(zz.imag)))
    ax.annotate(f"cruce en $r_0/\\ell$ = {xs[i]:.2f}\n$\\omega = \\pm{zz[i].real:.2f}"
                f" = \\pm 1/\\ell$", xy=(zz[i].real, 0), xytext=(-2.6, -2.5),
                fontsize=8.5, color=INK,
                arrowprops=dict(arrowstyle="->", color=INK, lw=.9))
    ax.annotate("los dos ceros colisionan\nsobre el eje imaginario",
                xy=(-0.2, 6.0), xytext=(-8.6, 6.9), fontsize=8.5, color=MUTED,
                arrowprops=dict(arrowstyle="->", color=MUTED, lw=.8))

    ax.set_xlabel("Re $\\omega$")
    ax.set_ylabel("Im $\\omega$")
    ax.set_title("C.  Lugar de raíces del par dominante\n"
                 "sube, cruza, colisiona y baja hacia el polo de A-L",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(-9, 9)
    ax.set_ylim(-3.2, 8.2)


def verifica():
    """La parte Im>0 de `ceros`, contra find_zeros_uhp de block1_linear."""
    for x in (5.0, 8.0):
        ell = 1.0 / x
        reg, pr = make_regulator("smeared", ell), Params(ell=ell)
        clave = lambda w: (round(w.real, 6), round(w.imag, 6))
        mio = sorted((w for w in ceros(reg, pr, 11.0, 1e-4, 9.0) if w.imag > 0),
                     key=clave)
        suyo = sorted(find_zeros_uhp(reg, pr, re_max=11.0, im_hi=9.0), key=clave)
        err = max((abs(a - b) for a, b in zip(mio, suyo)), default=np.inf)
        print(f"  r0/ell={x:4.1f}  este script: {len(mio)} ceros, "
              f"find_zeros_uhp: {len(suyo)}   max |dif| = {err:.2e}")


def main():
    os.chdir(RAIZ)
    warnings.filterwarnings("ignore")        # Newton avisa al tantear lejos
    z = locus(force="--force" in sys.argv)

    fig, axes = plt.subplots(1, 3, figsize=(16.5, 4.8))
    for ax, (x, letra, veredicto) in zip(axes, MAPAS):
        panel_mapa(ax, x, letra, veredicto)
    panel_locus(axes[2], z)

    for ax in axes:
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK)
        ax.yaxis.label.set_color(INK)
    axes[2].grid(True, color=GRID, lw=.8, zorder=0)
    axes[2].set_axisbelow(True)

    fig.suptitle("Ec. (14): los ceros de $\\chi^r_\\omega$ son los polos de la "
                 "respuesta, y deciden la estabilidad (regulador suavizado)",
                 x=0.008, ha="left", fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    fig.savefig("figures/seccionA_ec14_polos.png", dpi=170, facecolor="#fcfcfb")
    print("\nfigures/seccionA_ec14_polos.png")
    verifica()


if __name__ == "__main__":
    main()
