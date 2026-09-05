"""Trayectoria de la ecuacion LINEALIZADA y sus zonas de estabilidad.

Integra la ec. (6) linealizada alrededor de la carga en reposo hacia adelante
en tiempo propio, barriendo el cutoff a traves del borde de estabilidad, y
compara con el limite MARKOVIANO (local) de Abraham-Lorentz.

La distincion importa: linealizar quita la no linealidad, NO la memoria. La
ec. (14) sigue siendo no local. Lo markoviano es quedarse con el orden
dominante del desarrollo local, que es la ec. (15) y da el runaway.

Uso:  python3 scripts/fig_trayectoria_lineal.py
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np
from scipy.integrate import solve_ivp

from nlaid import RAIZ
from nlaid.core import Params, make_regulator
from nlaid.block1_linear import dispersion, integrate_linear

# Cutoffs a ambos lados del borde. El critico del regulador suavizado esta en
# r0/ell = 4.00 (block1_linear.critical_cutoff), asi que el barrido lo cruza.
CUTOFFS = [2.0, 3.0, 3.5, 4.0, 4.5, 5.0]
CRITICO = 4.00
S_END, TAU = 40.0, 2.0 / 3.0          # tau = 2 r0 / 3, tiempo de Abraham-Lorentz
XI_MAX = 0.32     # la linealizacion supone |xi| << r0: mas alla la curva no
                  # describe nada y solo tapa las demas. Se corta, no se recorta.
CACHE = "figures/trayectoria_lineal.npz"

# Rampa secuencial del repo: r0/ell es una MAGNITUD, no una identidad, asi que
# el color va de claro a oscuro y no por hues categoricos.
RAMPA = ["#86b6ef", "#5598e7", "#2a78d6", "#1c5cab", "#104281"]
INK, MUTED, GRID, ROJO = "#1a1a19", "#5c5b54", "#e4e3dd", "#c1442e"


def _colores(n):
    from matplotlib.colors import LinearSegmentedColormap
    return LinearSegmentedColormap.from_list("nlaid", RAMPA)(np.linspace(0, 1, n))


def envolvente(s, y):
    """Maximos locales de |y|. Dibujar la oscilacion cruda a 40 r0 con periodo
    ~1 satura de tinta y oculta la tasa, que es lo unico que se quiere leer."""
    a = np.maximum(np.abs(y), 1e-16)
    # Comparacion ESTRICTA a los dos lados: una meseta de valores identicos --
    # el tramo en reposo antes del pulso, donde xiddot es cero exacto -- no es
    # un maximo. Con >= producia una decena de picos falsos en el suelo, que
    # bastaban para burlar el umbral de abajo y dejaban la curva sin dibujar.
    i = np.where((a[1:-1] > a[:-2]) & (a[1:-1] > a[2:]))[0] + 1
    # Una senal monotona -- el runaway markoviano -- no tiene picos que
    # extraer: ahi la envolvente es la propia curva.
    return (s, a) if len(i) < 5 else (s[i], a[i])


def compute(force=False):
    """Trayectorias y cero dominante por cutoff. Se cachea: el barrido de
    raices sobre la malla compleja domina el costo del replot."""
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        if np.array_equal(d["cutoffs"], np.array(CUTOFFS, dtype=float)):
            return {k: d[k] for k in d.files}
    out = {"cutoffs": np.array(CUTOFFS, dtype=float)}
    for x in CUTOFFS:
        ell = 1.0 / x
        reg, pr = make_regulator("smeared", ell), Params(ell=ell)
        s, xi, _, xa = integrate_linear(reg, pr, s_end=S_END)
        z = _cero_dominante(reg, pr)
        out[f"s_{x}"], out[f"xi_{x}"], out[f"xa_{x}"] = s, xi, xa
        out[f"z_{x}"] = np.array([np.nan if z is None else z], dtype=complex)
        print(f"  r0/ell={x:4.1f}  cero dominante = {z}", flush=True)
    np.savez(CACHE, **out)
    return out


def markov(s_end=6.0, amp=0.05, ancho=1.0):
    """Limite local: xiddot = tau xidddot + k(s), ec. (15). Runaway e^{s/tau}.

    Se integra como problema de VALORES INICIALES, igual que la rama no local,
    para que la comparacion sea justa. (Con la prescripcion de Dirac -- imponer
    la condicion FINAL xiddot -> 0 -- se elimina el runaway a cambio de
    preaceleracion; no es lo que se compara aqui.)
    """
    pulso = lambda t: (amp * np.exp(-1.0 / (1.0 - (2 * t / ancho - 1) ** 2))
                       if 0.0 < t < ancho else 0.0)
    # y = (xi, xidot, xiddot);  xidddot = (xiddot - k)/tau
    f = lambda t, y: [y[1], y[2], (y[2] - pulso(t)) / TAU]
    sol = solve_ivp(f, (0.0, s_end), [0, 0, 0], max_step=1e-3, dense_output=True)
    return sol.t, sol.y[0], sol.y[2]


def main():
    os.chdir(RAIZ)
    d = compute(force="--force" in sys.argv)
    colores = _colores(len(CUTOFFS))
    fig, (axA, axB, axC) = plt.subplots(1, 3, figsize=(14.5, 4.5))

    for x, c in zip(CUTOFFS, colores):
        s, xi, xa = d[f"s_{x}"], d[f"xi_{x}"], d[f"xa_{x}"]
        z = complex(d[f"z_{x}"][0])
        lab = f"$r_0/\\ell$ = {x:.1f}"
        estable = x < CRITICO

        # Panel A: se corta donde |xi| sale del rango de validez lineal.
        j = np.argmax(np.abs(xi) > XI_MAX)
        j = len(s) if j == 0 else j
        axA.plot(s[:j], xi[:j], lw=1.9 if estable else 1.5, color=c, label=lab,
                 ls="-" if estable else (0, (5, 2)))
        if j < len(s):
            axA.plot([s[j - 1]], [xi[j - 1]], "o", ms=5, color=c,
                     mec="#fcfcfb", mew=1.2, zorder=6)

        se, ae = envolvente(s, xa)
        if np.isfinite(z.real):
            m = se > 4.0
            i0 = np.searchsorted(se, 5.0)
            if m.sum() > 3 and i0 < len(se):
                axB.semilogy(se[m], ae[i0] * np.exp(z.imag * (se[m] - se[i0])),
                             lw=5, color=c, alpha=.30, zorder=2, solid_capstyle="round")
        axB.semilogy(se, ae, lw=1.4, color=c, label=lab, zorder=3)

    axA.axhline(0, color=GRID, lw=1.2, zorder=0)
    axA.set_xlabel("$s / r_0$"); axA.set_ylabel("$\\xi(s) / r_0$")
    axA.set_title("A.  Trayectoria $\\xi(s)$\n"
                  "continua: relaja   discontinua: $r_0/\\ell \\geq$ crítico",
                  loc="left", fontsize=10.5, color=INK, pad=8)
    axA.set_xlim(-2, S_END); axA.set_ylim(-0.36, 0.36)
    axA.legend(frameon=False, fontsize=8.5, ncol=2, loc="lower right")

    axB.set_xlabel("$s / r_0$"); axB.set_ylabel("envolvente de $|\\ddot{\\xi}|$")
    axB.set_title("B.  Tasa de relajación\n"
                  "(banda gruesa: predicción espectral, ceros de $D$)",
                  loc="left", fontsize=10.5, color=INK, pad=8)
    axB.set_xlim(-2, S_END); axB.set_ylim(1e-11, 1e4)
    axB.legend(frameon=False, fontsize=8.5, ncol=2, loc="upper right")

    # --- Panel C: no local frente al limite markoviano -------------------
    sm, xim, xam = markov()
    s, xa = d["s_3.0"], d["xa_3.0"]
    axC.semilogy(*envolvente(sm, xam), lw=2, color=ROJO,
                 label="markoviano, ec. (15)")
    axC.semilogy(*envolvente(s, xa), lw=2, color=RAMPA[3],
                 label="no local, $r_0/\\ell$ = 3")
    axC.plot(sm, 0.02 * np.exp(sm / TAU), lw=1, ls=(0, (4, 3)), color=MUTED)
    axC.annotate("$e^{\\,s/\\tau}$,  $\\tau = 2r_0/3$", xy=(2.6, 1.2e3),
                 fontsize=9, color=MUTED)
    axC.set_xlabel("$s / r_0$"); axC.set_ylabel("$|\\ddot{\\xi}|$")
    axC.set_title("C.  El límite local NO se estabiliza\n"
                  "(mismo pulso externo en los dos)",
                  loc="left", fontsize=10.5, color=INK, pad=8)
    axC.set_xlim(0, 12); axC.set_ylim(1e-9, 1e6)
    axC.legend(frameon=False, fontsize=9, loc="lower right")

    for ax in (axA, axB, axC):
        ax.grid(True, color=GRID, lw=.8, zorder=0); ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)

    fig.suptitle("Ecuación linealizada: trayectoria y borde de estabilidad "
                 f"(regulador suavizado, crítico en $r_0/\\ell$ = {CRITICO:.0f})",
                 x=0.008, ha="left", fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    fig.savefig("figures/trayectoria_lineal.png", dpi=170, facecolor="#fcfcfb")
    print("\nfigures/trayectoria_lineal.png")


def _cero_dominante(reg, pr):
    """Cero de D con mayor Im, por minimos de |D| sobre una malla + Newton."""
    from scipy.optimize import newton
    R, I = np.meshgrid(np.linspace(-9, 9, 240), np.linspace(-3, 5, 200), indexing="ij")
    W = R + 1j * I
    A = np.abs(dispersion(W.ravel(), reg, pr)).reshape(W.shape)
    seeds = sorted((W[i, j] for i in range(1, A.shape[0] - 1)
                    for j in range(1, A.shape[1] - 1)
                    if A[i, j] <= A[i - 1:i + 2, j - 1:j + 2].min() and W[i, j].real > 1e-6),
                   key=lambda z: -z.imag)
    for s0 in seeds[:15]:
        try:
            r = newton(lambda w: dispersion(w, reg, pr), s0, tol=1e-12, maxiter=90)
        except Exception:                                     # noqa: BLE001
            continue
        if abs(dispersion(r, reg, pr)) < 1e-8 and r.real > 1e-6:
            return r
    return None


if __name__ == "__main__":
    main()
