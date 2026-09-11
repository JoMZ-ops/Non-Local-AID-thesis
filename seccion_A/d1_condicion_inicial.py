"""D1, hipotesis 2: ¿el caracter oscilatorio depende de la condicion inicial?

El paper (p. 10) prescribe una trayectoria y apaga la fuente:

    "We impose the initial condition that the charge is at rest, x(t) = (t,0),
     for t < t_0 and the charge follows a prescribed trajectory, x_i(t) for
     t_0 < t < t_0 + t_i. A certain external source, k_i(s), is supposed to
     generate this motion which is turned off after this initial phase and the
     invariant length, s, of the world line is measured from t_0 + t_i."

El repo usa en cambio un pulso C^infinito (ver `worldline.smooth_bump`). La
hipotesis 2 de `docs/discrepancias.md` era que esa diferencia explicaba D1: que
la Fig. 1(a) relaje de forma monotona y la nuestra oscile.

La prediccion contraria es concreta y falsable: si la oscilacion es propiedad
de la ECUACION -- del cero dominante de chi -- y no de la excitacion, entonces
tasa y periodo deben salir iguales para CUALQUIER condicion inicial, porque la
condicion inicial solo fija la amplitud y la fase de cada modo, no cuales
existen. Eso es lo que mide este script.

Se prueban tres familias:
  - pulso C^infinito, variando amplitud y ancho;
  - trayectoria prescrita: tramo de aceleracion propia CONSTANTE (movimiento
    hiperbolico, solucion exacta de xdot^2 = 1) durante t_i, y fuente apagada
    en s = 0. Es la lectura literal del parrafo del paper, con su salto de
    aceleracion en el empalme incluido;
  - la misma, con amplitudes grandes, para asomarse al regimen no lineal.

Uso:  python3 seccion_A/d1_condicion_inicial.py [--force]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid import RAIZ
from nlaid.core import Params
from nlaid.block3_memory import integrate_memory
from nlaid.worldline import WorldLine, smooth_bump

CACHE = RAIZ / "figures/seccionA_d1.npz"
LOCUS = RAIZ / "figures/seccionA_ec14_locus.npz"
ELL, M_OVER_MB = 1.0 / 3.0, 2.0        # r0/ell = 3: el caso de la Fig. 1
S_END, DS, N_ELL = 25.0, 5e-3, 25.0
SIN_FUENTE = lambda s: 0.0

# (etiqueta, tipo, parametros)
CASOS = [
    ("pulso A=0.05, w=1",   "pulso",      dict(amplitude=0.05, width=1.0)),
    ("pulso A=0.30, w=1",   "pulso",      dict(amplitude=0.30, width=1.0)),
    ("pulso A=0.90, w=1",   "pulso",      dict(amplitude=0.90, width=1.0)),
    ("pulso A=0.30, w=3",   "pulso",      dict(amplitude=0.30, width=3.0)),
    ("prescrita a=0.3, t=1", "prescrita", dict(a0=0.3, t_i=1.0)),
    ("prescrita a=0.3, t=4", "prescrita", dict(a0=0.3, t_i=4.0)),
    ("prescrita a=1.0, t=1", "prescrita", dict(a0=1.0, t_i=1.0)),
    ("prescrita a=2.0, t=2", "prescrita", dict(a0=2.0, t_i=2.0)),
]

RAMPA = ["#86b6ef", "#5598e7", "#2a78d6", "#104281"]
NARANJA = ["#f5b78f", "#f09a6e", "#eb6834", "#8f3612"]
INK, MUTED, GRID, VERDE = "#1a1a19", "#5c5b54", "#e4e3dd", "#2e7d5b"


def historia_prescrita(a0, t_i, dim=2, ds=DS, s_rest=None):
    """Reposo, luego aceleracion propia constante a0 durante t_i, hasta s = 0.

    El tramo acelerado es movimiento hiperbolico, que cumple xdot^2 = 1
    EXACTAMENTE, de modo que la prehistoria no introduce deriva propia. En el
    empalme la aceleracion salta de 0 a a0: es justo lo que describe el paper
    ("a certain external source ... turned off after this initial phase"), y lo
    que el pulso C^infinito evita a proposito.
    """
    s_rest = s_rest if s_rest is not None else max(4.0, 1.5 * N_ELL * ELL) + t_i
    wl = WorldLine(dim)
    v0 = np.zeros(dim); v0[0] = 1.0
    cero = np.zeros(dim)

    for sv in np.arange(-s_rest, -t_i, ds):          # reposo
        x = np.zeros(dim); x[0] = float(sv)
        wl.append(float(sv), x, v0, cero)

    x_emp = np.zeros(dim); x_emp[0] = -t_i           # empalme
    for sv in np.arange(-t_i, 0.0 + 0.5 * ds, ds):
        sig = float(sv) + t_i
        x = x_emp.copy()
        x[0] += np.sinh(a0 * sig) / a0
        x[1] += (np.cosh(a0 * sig) - 1.0) / a0
        v = np.zeros(dim); v[0], v[1] = np.cosh(a0 * sig), np.sinh(a0 * sig)
        a = np.zeros(dim); a[0], a[1] = a0 * np.sinh(a0 * sig), a0 * np.cosh(a0 * sig)
        wl.append(float(sv), x, v, a)
    return wl


def corre(tipo, kw):
    p = Params(ell=ELL, m_over_mB=M_OVER_MB)
    if tipo == "pulso":
        wl = integrate_memory(p, s_end=S_END, ds=DS, n_ell=N_ELL, pts_per_ell=16,
                              drive=lambda s: smooth_bump(s, **kw))
    else:
        wl = integrate_memory(p, s_end=S_END, ds=DS, n_ell=N_ELL, pts_per_ell=16,
                              history=historia_prescrita(**kw), drive=SIN_FUENTE)
    m = wl.s >= 0.0
    return wl.s[m], np.linalg.norm(wl.a[m], axis=1), float(wl.norm_drift.max())


def mide(s, a, desde=0.55):
    """Tasa y semiperiodo de la envolvente, ajustados sobre los maximos locales."""
    m = s >= desde * s.max()
    ss, aa = s[m], np.maximum(a[m], 1e-300)
    i = np.where((aa[1:-1] > aa[:-2]) & (aa[1:-1] > aa[2:]))[0] + 1
    if len(i) < 4:
        return np.nan, np.nan
    tasa = float(np.polyfit(ss[i], np.log(aa[i]), 1)[0])
    return tasa, float(np.mean(np.diff(ss[i])))


def compute(force=False):
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE, allow_pickle=True)
        if list(d["etiquetas"]) == [c[0] for c in CASOS]:
            return {k: d[k] for k in d.files}

    out = {"etiquetas": np.array([c[0] for c in CASOS])}
    tasas, periodos, derivas = [], [], []
    for etq, tipo, kw in CASOS:
        s, a, drift = corre(tipo, kw)
        tasa, per = mide(s, a)
        out[f"s_{etq}"], out[f"a_{etq}"] = s, a
        tasas.append(tasa); periodos.append(per); derivas.append(drift)
        print(f"  {etq:24s} tasa={tasa:+.4f}  semiperiodo={per:.4f}  "
              f"deriva={drift:.1e}", flush=True)
    out["tasas"], out["periodos"] = np.array(tasas), np.array(periodos)
    out["derivas"] = np.array(derivas)
    np.savez(CACHE, **out)
    return out


def espectral():
    """Cero dominante de chi a r0/ell = 3, del cache de la ec. (14)."""
    if not os.path.exists(LOCUS):
        sys.exit(f"falta {LOCUS}: corra antes seccion_A/fig_ec14_polos.py")
    d = np.load(LOCUS)
    j = int(np.argmin(np.abs(d["xs"] - 3.0)))
    return complex(d["z"][j])


def main():
    d = compute(force="--force" in sys.argv)
    z = espectral()
    etiquetas = list(d["etiquetas"])
    colores = RAMPA + NARANJA

    fig, (axA, axB) = plt.subplots(1, 2, figsize=(13.0, 5.0))

    for etq, c in zip(etiquetas, colores):
        s, a = d[f"s_{etq}"], d[f"a_{etq}"]
        axA.semilogy(s, np.maximum(a, 1e-14), lw=1.4, color=c, label=etq)
    ss = np.linspace(6, S_END, 50)
    axA.semilogy(ss, 0.02 * np.exp(z.imag * (ss - 6)), lw=6, color=VERDE,
                 alpha=.25, zorder=1, solid_capstyle="round")
    axA.annotate(f"pendiente $e^{{{z.imag:.3f}s}}$: el cero\ndominante de "
                 "$\\chi$, ec. (14)", xy=(1.0, 2e-11), fontsize=8.5, color=VERDE)
    axA.set_xlabel("$s / r_0$"); axA.set_ylabel("$|\\ddot{x}|\\, r_0$")
    axA.set_title("A.  Ocho condiciones iniciales distintas\n"
                  "pulso $C^\\infty$ y trayectoria prescrita, misma pendiente",
                  loc="left", fontsize=11, color=INK, pad=10)
    axA.set_xlim(0, S_END); axA.set_ylim(1e-13, 5)
    axA.legend(frameon=False, fontsize=7.5, loc="upper right", ncol=2)

    # Barras desde cero harian que ocho valores identicos se vieran identicos
    # por la escala y no por el resultado. Se grafica la DESVIACION contra la
    # prediccion espectral, con el eje apretado a 1e-3: ahi si se ve que la
    # dispersion es de 2e-4 y no un artefacto de dibujo.
    t, per = d["tasas"], d["periodos"]
    y = np.arange(len(etiquetas))
    axB.axvline(0, color=VERDE, lw=2, zorder=2)
    axB.plot(t - z.imag, y, "o", ms=9, color="none", zorder=3)
    for k, c in enumerate(colores):
        axB.plot([t[k] - z.imag], [y[k]], "o", ms=9, color=c, mec="#fcfcfb",
                 mew=1.2, zorder=4)
        axB.annotate(f"$T/2$ = {per[k]:.4f}", xy=(6.5e-4, y[k]), fontsize=8,
                     color=INK, va="center")
    axB.annotate(f"predicción espectral\nIm $\\omega$ = {z.imag:.4f}",
                 xy=(0, -0.75), ha="center", fontsize=8.5, color=VERDE)
    axB.set_yticks(y); axB.set_yticklabels(etiquetas, fontsize=8)
    axB.invert_yaxis()
    axB.set_xlabel("tasa medida $-$ predicción espectral")
    axB.set_title(f"B.  Desviación contra la ec. (14)\n"
                  f"semiperíodo predicho: $\\pi/{z.real:.4f}$ = "
                  f"{np.pi/z.real:.4f}",
                  loc="left", fontsize=11, color=INK, pad=10)
    axB.set_xlim(-6e-4, 1.5e-3)
    axB.ticklabel_format(axis="x", style="sci", scilimits=(0, 0))

    for ax in (axA, axB):
        ax.grid(True, color=GRID, lw=.8, zorder=0)
        ax.set_axisbelow(True)
        for sp in ("top", "right"):
            ax.spines[sp].set_visible(False)
        for sp in ("left", "bottom"):
            ax.spines[sp].set_color(GRID)
        ax.tick_params(colors=MUTED, labelsize=9)
        ax.xaxis.label.set_color(INK); ax.yaxis.label.set_color(INK)

    fig.suptitle("D1, hipótesis 2: la condición inicial NO cambia el carácter "
                 "de la relajación", x=0.008, ha="left", fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.92))
    salida = RAIZ / "figures/seccionA_d1_condicion_inicial.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)

    print(f"\n  espectral (ec. 14):  tasa = {z.imag:.4f}   "
          f"semiperiodo = {np.pi/z.real:.4f}")
    print(f"  medido, las {len(t)} condiciones: tasa en "
          f"[{np.nanmin(t):.4f}, {np.nanmax(t):.4f}]   "
          f"semiperiodo en [{np.nanmin(per):.4f}, {np.nanmax(per):.4f}]")
    print(f"  dispersion relativa: tasa {np.nanstd(t)/abs(np.nanmean(t)):.2%}   "
          f"semiperiodo {np.nanstd(per)/np.nanmean(per):.2%}")


if __name__ == "__main__":
    main()
