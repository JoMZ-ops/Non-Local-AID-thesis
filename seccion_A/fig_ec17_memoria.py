"""Ec. (17): el regulador suavizado y su memoria infinita.

    xddot = (r0 / 3 ell^4)(m/m_B) int_{-inf}^{0} du
            (1 - sqrt(w)/2ell) e^{-sqrt(w)/ell}
            { (x-x')(xdot.xdot') + [xdot.(x'-x)] xdot' },   w = (x-x')^2   (17)

Frente a la ec. (16), que consulta UN punto del pasado -- el que cumple
ell^2 = (x-x')^2 --, la (17) integra sobre TODO el pasado con peso
e^{-sqrt(w)/ell}. "Memoria infinita" es literal: el limite inferior es -inf.
Estas tres graficas responden qué cuesta eso.

Panel A  de dónde viene la fuerza: |nucleo * V| contra u/ell a lo largo de una
         trayectoria real, con la fraccion acumulada. La memoria que pesa NO es
         fija: se alarga conforme la trayectoria relaja, porque al enderezarse
         la linea de mundo el pasado reciente aporta cada vez menos.
Panel B  qué pasa si se trunca: |xddot| con ventanas n_ell = 2, 4, 8, 24
         contra la referencia n_ell = 35.
Panel C  el error de truncamiento contra n_ell, con su ajuste exponencial. Es
         lo que justifica el n_ell = 30 por defecto de `block3_memory`.

Uso:  python3 seccion_A/fig_ec17_memoria.py [--force]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid import RAIZ
from nlaid.core import Params, make_regulator, minkowski_dot
from nlaid.block3_memory import integrate_memory

CACHE = RAIZ / "figures/seccionA_ec17.npz"
ELL = 1.0 / 3.0                 # r0/ell = 3, el caso de la Fig. 1 del paper
M_OVER_MB = 2.0                 # el valor del contratermino ahi: 1/(1 - 3/6)
S_END, DS = 15.0, 5e-3
VENTANAS = [2, 4, 8, 24]        # n_ell a comparar
N_REF = 35                      # referencia: peso residual ~ e^-35
INSTANTES = [0.5, 3.0, 8.0, 10.0, 14.0]     # s donde se mira el peso de la memoria

RAMPA = ["#86b6ef", "#5598e7", "#2a78d6", "#104281"]
INK, MUTED, GRID, ROJO = "#1a1a19", "#5c5b54", "#e4e3dd", "#c1442e"


def perfil(n_ell):
    """|xddot|(s) para una ventana de memoria dada."""
    wl = integrate_memory(Params(ell=ELL, m_over_mB=M_OVER_MB),
                          s_end=S_END, ds=DS, n_ell=float(n_ell), pts_per_ell=16)
    return wl, np.linalg.norm(wl.a, axis=1)


def peso_memoria(wl, s_now, n_ell=25.0, n=800):
    """|nucleo * V| contra u, el integrando de la ec. (17) en s = s_now.

    Se reconstruye con las mismas piezas que `rhs_memory`: delta_B' del
    regulador de `core` y el mismo V. No se reescribe el nucleo aqui.
    """
    reg = make_regulator("smeared", ELL)
    x, v, _ = wl.sample(s_now)
    u = np.linspace(-n_ell * ELL, -1e-9, n)
    xp, vp, _ = wl.sample_many(s_now + u)
    d = x - xp
    w = np.maximum(minkowski_dot(d, d), 0.0)
    v_vp, v_d = minkowski_dot(v, vp), minkowski_dot(v, d)
    V = d * v_vp[:, None] - v_d[:, None] * vp
    return u, np.abs(reg.d_delta(w)[:, None] * V).sum(axis=1)


def compute(force=False):
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        # La receta incluye INSTANTES y S_END: sin ellos, agregar un instante
        # devolvia el cache viejo y la figura no reflejaba el cambio.
        if (np.array_equal(d["ventanas"], np.array(VENTANAS))
                and np.array_equal(d["instantes"], np.array(INSTANTES))
                and float(d["s_end"]) == S_END):
            return {k: d[k] for k in d.files}

    out = {"ventanas": np.array(VENTANAS), "instantes": np.array(INSTANTES),
           "s_end": np.array(S_END)}
    wl_ref, a_ref = perfil(N_REF)
    sg = np.linspace(0.0, S_END, 1200)
    out["sg"] = sg
    out["a_ref"] = np.interp(sg, wl_ref.s, a_ref)
    print(f"  referencia n_ell={N_REF}: |a|max={out['a_ref'].max():.3e}")

    for s_now in INSTANTES:
        u, c = peso_memoria(wl_ref, s_now)
        out[f"u_{s_now}"], out[f"c_{s_now}"] = u, c

    # La prehistoria depende de n_ell, asi que cada corrida trae su propia
    # malla: se comparan interpolando sobre `sg`, no indice a indice.
    for n in VENTANAS + [3, 6, 12, 16]:
        wl, a = perfil(n)
        ai = np.interp(sg, wl.s, a)
        out[f"a_{n}"] = ai
        out[f"err_{n}"] = np.abs(ai - out["a_ref"]).max()
        print(f"  n_ell={n:3d}  err_max={out[f'err_{n}']:.3e}", flush=True)

    np.savez(CACHE, **out)
    return out


def panel_peso(ax, d):
    """Fraccion acumulada del peso, contando desde u = 0 hacia el pasado.

    Se grafica la acumulada y no el integrando: |delta_B' V| oscila con la
    trayectoria (que relaja oscilando), y normalizado por su maximo son tres
    curvas multilobuladas de las que no se lee nada. La acumulada es monotona
    y responde la pregunta directamente: cuanto pasado hace falta.
    """
    from matplotlib.colors import LinearSegmentedColormap
    # Un color POR INSTANTE. Con una lista fija, `zip` descartaba en silencio
    # los instantes sobrantes: agregar uno a INSTANTES no dibujaba nada.
    colores = LinearSegmentedColormap.from_list("nlaid", RAMPA)(
        np.linspace(0.15, 1.0, len(INSTANTES)))
    alcances = []
    for s_now, c in zip(INSTANTES, colores):
        u, w = d[f"u_{s_now}"], d[f"c_{s_now}"]
        acum = np.concatenate([[0], np.cumsum(np.diff(u) * (w[1:] + w[:-1]) / 2)])
        desde_cero = 1.0 - acum / acum[-1]          # peso acumulado en [u, 0]
        ax.plot(u / ELL, desde_cero, lw=2, color=c, label=f"$s$ = {s_now}")
        j = int(np.searchsorted(acum / acum[-1], 0.10))
        ax.plot([u[j] / ELL], [0.9], "o", ms=6, color=c, mec="#fcfcfb",
                mew=1.2, zorder=6)
        alcances.append(-u[j] / ELL)

    ax.axhline(0.9, color=MUTED, lw=1, ls=(0, (2, 3)))
    ax.annotate("90% del peso", xy=(-14.6, 0.915), fontsize=8.5, color=MUTED)
    ax.annotate("al relajar, la memoria que importa se alarga:\n"
                f"de ${alcances[0]:.1f}\\ell$ a ${alcances[-1]:.1f}\\ell$   "
                "(los puntos marcan el 90%)",
                xy=(-14.6, 0.30), fontsize=8.5, color=INK)
    ax.set_xlabel("$u/\\ell$        ($u = s' - s \\leq 0$)")
    ax.set_ylabel("fracción del peso acumulada en $[u,\\,0]$")
    ax.set_title("A.  Cuánto pasado hace falta, ec. (17)\n"
                 "el límite es $-\\infty$, pero el peso se concentra",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(-15, 0)
    ax.set_ylim(0, 1.05)
    ax.legend(frameon=False, fontsize=8.5, loc="center left")


def panel_truncado(ax, d):
    ax.semilogy(d["sg"], np.maximum(d["a_ref"], 1e-12), lw=2.6, color=INK,
                label=f"$n_\\ell$ = {N_REF} (referencia)", zorder=5)
    for n, c in zip(VENTANAS, RAMPA):
        ax.semilogy(d["sg"], np.maximum(d[f"a_{n}"], 1e-12), lw=1.6, color=c,
                    ls=(0, (5, 2)), label=f"$n_\\ell$ = {n}")
    ax.annotate("truncar a $2\\ell$ cambia la\ndinámica, no solo la precisión",
                xy=(7.5, 3e-3), fontsize=8.5, color=INK)
    ax.set_xlabel("$s / r_0$")
    ax.set_ylabel("$|\\ddot{x}|\\; r_0$")
    ax.set_title("B.  Qué cuesta truncar la memoria\n"
                 f"$r_0/\\ell$ = 3, $m/m_B$ = {M_OVER_MB} (Fig. 1 del paper)",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, S_END)
    ax.set_ylim(1e-6, 1)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right", ncol=2)


def panel_error(ax, d):
    ns = np.array(sorted(VENTANAS + [3, 6, 12, 16]))
    err = np.array([float(d[f"err_{n}"]) for n in ns])
    ax.semilogy(ns, err, "o-", lw=1.8, ms=6, color=RAMPA[3])

    m = ns >= 6                                     # el regimen asintotico
    p = np.polyfit(ns[m], np.log(err[m]), 1)
    ax.semilogy(ns, np.exp(np.polyval(p, ns)), lw=1.2, ls=(0, (4, 3)),
                color=ROJO,
                label=f"ajuste $\\propto e^{{{p[0]:.2f}\\,n_\\ell}}$")
    ax.axvline(30, color=MUTED, lw=1, ls=(0, (2, 3)))
    ax.annotate("$n_\\ell = 30$, el valor\npor defecto del bloque 3",
                xy=(30, 3e-4), xytext=(-8, 0), textcoords="offset points",
                fontsize=8.5, color=MUTED, ha="right")
    ax.set_xlabel("$n_\\ell$        (ventana de memoria, en unidades de $\\ell$)")
    ax.set_ylabel("$\\max_s\\,|\\,|\\ddot{x}|_{n_\\ell} - |\\ddot{x}|_{\\rm ref}|$")
    ax.set_title("C.  El error de truncamiento es exponencial\n"
                 "la memoria es infinita, pero converge rápido",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, 34)
    ax.legend(frameon=False, fontsize=9, loc="lower left")


def main():
    d = compute(force="--force" in sys.argv)
    fig, axes = plt.subplots(1, 3, figsize=(16.0, 4.8))
    panel_peso(axes[0], d)
    panel_truncado(axes[1], d)
    panel_error(axes[2], d)

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

    fig.suptitle("Ec. (17): regulador suavizado — la memoria es infinita, y "
                 "eso se puede medir", x=0.008, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salida = RAIZ / "figures/seccionA_ec17_memoria.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)


if __name__ == "__main__":
    main()
