"""Ec. (16): el regulador desplazado y su retardo finito.

    xddot = r0 (m/m_B) / A^2 [ (B+1)/A * V1 + V2 ]                        (16)
    A = xdot'.(x-x'),  B = xddot'.(x'-x)
    V1 = (x-x')(xdot.xdot') - [xdot.(x-x')] xdot'
    V2 = (x-x')(xdot.xddot') + [xdot.(x'-x)] xddot'

con el punto retardado x' fijado por  ell^2 = (x - x')^2  (p. 9 del paper).
Esa condicion -- un intervalo INVARIANTE, no un retardo de coordenada -- es
todo el contenido geometrico del regulador desplazado, y es lo que el panel A
dibuja.

Panel A  la construccion del punto retardado en el plano (x^1, x^0): cono de
         luz, hiperbola de distancia invariante ell y la interseccion con la
         linea de mundo. El punto retardado esta DENTRO del cono, a intervalo
         temporal ell: ahi esta la supresion de la interaccion superluminica.
Panel B  el retardo, contra s. El invariante vale ell siempre (por
         construccion); lo que varia es su reparto en (Delta t, Delta x), y el
         retardo PROPIO Delta s se aparta de ell solo mientras hay aceleracion.
         Por eso la (16) es una DDE de retardo dependiente del estado.
Panel C  envolvente de |xddot| en CUATRO puntos del plano (r0/ell, m/m_B).
         El cutoff solo no decide: a r0/ell = 0.8 la trayectoria relaja con
         m/m_B = 1 y crece con el m/m_B = 5/3 que fija el contratermino de la
         ec. (10). Sobre esa linea del contratermino -- la punteada de la
         Fig. 2 -- el borde cae donde lo pone el bloque 1, r0/ell = 0.669, y
         la tasa medida se contrasta contra Im omega del cero de chi.

Nota de comparabilidad: `susceptibility` (ec. 14) NO lee m_over_mB, mientras
que `integrate_delay` (ec. 16) multiplica el lado derecho por el. Comparar la
(16) a m/m_B = 1 contra el critico linealizado es comparar dos puntos
distintos del plano de fases; por eso el panel C corre las dos columnas.

Uso:  python3 seccion_A/fig_ec16_retardo.py [--force]
"""

import os
import sys

import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
import numpy as np

from nlaid import RAIZ
from nlaid.block1_linear import find_zeros_uhp
from nlaid.core import Params, make_regulator
from nlaid.block2_delay import integrate_delay
from nlaid.worldline import smooth_bump

CACHE = RAIZ / "figures/seccionA_ec16.npz"
CUTOFFS = [0.4, 0.8]                     # a ambos lados del critico linealizado
CRITICO_LINEAL = 0.669
ELL_GEO = 1.0                            # r0/ell = 1 para los paneles A y B
S_END = 25.0
DERIVA_MAX = 1e-3

RAMPA = ["#86b6ef", "#5598e7", "#1c5cab", "#104281"]
INK, MUTED, GRID, ROJO = "#1a1a19", "#5c5b54", "#e4e3dd", "#c1442e"

# Pulso fuerte para la geometria (hace visible la curvatura de la linea de
# mundo) y el pulso por defecto del bloque 2 para el barrido, que es el que se
# puede comparar contra la teoria linealizada.
PULSO_GEO = lambda s: smooth_bump(s, amplitude=0.9, s0=0.0, width=1.5)


def compute(force=False):
    if os.path.exists(CACHE) and not force:
        d = np.load(CACHE)
        if np.array_equal(d["cutoffs"], np.array(CUTOFFS, dtype=float)):
            return {k: d[k] for k in d.files}

    out = {"cutoffs": np.array(CUTOFFS, dtype=float)}

    # --- geometria: una sola corrida, r0/ell = 1 -------------------------
    wl = integrate_delay(Params(ell=ELL_GEO), s_end=12.0, ds=5e-3, drive=PULSO_GEO)
    sg = np.linspace(-1.5, 9.0, 220)
    ret = np.array([wl.retarded_point(q, wl.sample(q)[0], ELL_GEO)[0] for q in sg])
    xq = np.array([wl.sample(q)[0] for q in sg])
    xr = np.array([wl.sample(q)[0] for q in ret])
    out["geo_s"], out["geo_ret"] = sg, ret
    out["geo_x"], out["geo_xr"] = xq, xr
    m = wl.s >= -2.0
    out["geo_linea"] = wl.x[m]
    print(f"  geometria: max |(x-x')^2 - ell^2| = "
          f"{np.abs(((xq-xr)[:,0]**2 - (xq-xr)[:,1]**2) - ELL_GEO**2).max():.2e}")

    # --- barrido: dos cutoffs x dos valores de m/m_B ---------------------
    for x in CUTOFFS:
        ell = 1.0 / x
        contraterm = make_regulator("shifted", ell).m_over_mB_counterterm()
        out[f"mr_{x}"] = contraterm
        for etq, mr in (("libre", 1.0), ("contraterm", contraterm)):
            wl = integrate_delay(Params(ell=ell, m_over_mB=mr),
                                 s_end=S_END, ds=min(0.05 * ell, 5e-3))
            m = wl.s >= 0.0
            out[f"s_{x}_{etq}"] = wl.s[m]
            out[f"a_{x}_{etq}"] = np.linalg.norm(wl.a[m], axis=1)
            out[f"d_{x}_{etq}"] = wl.norm_drift[m]
            print(f"  r0/ell={x:4.2f}  m/m_B={mr:6.3f} ({etq:10s})  "
                  f"|a|fin={out[f'a_{x}_{etq}'][-1]:.2e}  "
                  f"deriva max={wl.norm_drift.max():.1e}", flush=True)

    # Ceros de chi (ec. 14) para el mismo regulador: la via independiente
    # contra la que se contrasta la tasa medida.
    for x in CUTOFFS:
        ell = 1.0 / x
        zs = find_zeros_uhp(make_regulator("shifted", ell), Params(ell=ell),
                            re_max=45.0, im_hi=6.0, grid=240)
        z = max(zs, key=lambda w: w.imag) if zs else np.nan + 0j
        out[f"z_{x}"] = np.array([z], dtype=complex)
        print(f"  r0/ell={x:4.2f}  cero dominante de chi = {z}")

    np.savez(CACHE, **out)
    return out


def envolvente(s, y):
    """Máximos locales de |y|: la oscilación cruda a omega ~ 8 satura de tinta."""
    a = np.maximum(np.abs(y), 1e-16)
    i = np.where((a[1:-1] > a[:-2]) & (a[1:-1] > a[2:]))[0] + 1
    return (s, a) if len(i) < 5 else (s[i], a[i])


def panel_geometria(ax, d):
    linea, s, xq, xr = d["geo_linea"], d["geo_s"], d["geo_x"], d["geo_xr"]
    ax.plot(linea[:, 1], linea[:, 0], lw=2.4, color=INK, zorder=5,
            label="línea de mundo $x(s)$")

    # En cada s elegido: el cono de luz pasado, la hiperbola de distancia
    # invariante ell, y el punto retardado donde esa hiperbola corta la linea
    # de mundo. Dibujar el SEGMENTO x -> x' no serviria: es una cuerda de la
    # propia linea de mundo y se confunde con ella.
    dx = np.linspace(-1.9, 1.9, 300)
    t = np.linspace(0, 2.1, 100)
    for q in (0.0, 2.0, 4.0):
        i = int(np.argmin(np.abs(s - q)))
        x0, xp = xq[i], xr[i]
        ax.plot(x0[1] - t, x0[0] - t, lw=.9, ls=(0, (4, 3)), color=MUTED, zorder=3)
        ax.plot(x0[1] + t, x0[0] - t, lw=.9, ls=(0, (4, 3)), color=MUTED, zorder=3)
        ax.plot(x0[1] + dx, x0[0] - np.sqrt(ELL_GEO ** 2 + dx ** 2), lw=1.5,
                color=ROJO, zorder=4)
        ax.plot([x0[1]], [x0[0]], "o", ms=6, color=INK, mec="#fcfcfb",
                mew=1.2, zorder=8)
        ax.plot([xp[1]], [xp[0]], "o", ms=6, color=RAMPA[1], mec="#fcfcfb",
                mew=1.2, zorder=8)

    ax.plot([], [], "o", ms=6, color=INK, label="$x(s)$")
    ax.plot([], [], "o", ms=6, color=RAMPA[1], label="$x'$: la hipérbola corta $x(s)$")
    ax.plot([], [], lw=1.5, color=ROJO, label="$(x-x')^2 = \\ell^2$")
    ax.plot([], [], lw=.9, ls=(0, (4, 3)), color=MUTED, label="cono de luz")

    ax.annotate("$x'$ queda DENTRO del cono,\na intervalo temporal $\\ell$:\n"
                "sin interacción superlumínica",
                xy=(0.55, -1.55), fontsize=8.5, color=INK)
    ax.set_xlabel("$x^1$")
    ax.set_ylabel("$x^0$")
    ax.set_title("A.  La condición $\\ell^2 = (x-x')^2$, ec. (16)\n"
                 "distancia invariante, no retardo de coordenada",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(-2.2, 3.0)
    ax.set_ylim(-1.8, 7.2)
    ax.legend(frameon=False, fontsize=8, loc="upper left")


def panel_retardo(ax, d):
    s, ret, xq, xr = d["geo_s"], d["geo_ret"], d["geo_x"], d["geo_xr"]
    dif = xq - xr
    ax.plot(s, (s - ret) / ELL_GEO, lw=2.2, color=RAMPA[3],
            label="retardo propio  $\\Delta s/\\ell$")
    ax.plot(s, dif[:, 0] / ELL_GEO, lw=2, color=ROJO,
            label="retardo de coordenada  $\\Delta t/\\ell$")
    ax.plot(s, dif[:, 1] / ELL_GEO, lw=2, ls=(0, (5, 2)), color=RAMPA[0],
            label="separación espacial  $\\Delta x/\\ell$")
    ax.axhline(1.0, color=MUTED, lw=1, ls=(0, (2, 3)))
    ax.axvspan(0.0, 1.5, color=RAMPA[0], alpha=.18, lw=0)
    ax.annotate("pulso externo", xy=(0.1, 1.52), fontsize=8.5, color=MUTED)
    ax.annotate("el invariante vale $\\ell$ exacto en todo $s$\n"
                "(residuo $< 10^{-9}$); lo que cambia es\ncómo se reparte",
                xy=(3.0, 0.22), fontsize=8.5, color=INK)

    ax.set_xlabel("$s / r_0$")
    ax.set_ylabel("en unidades de $\\ell$")
    ax.set_title("B.  El retardo depende del estado, a $r_0/\\ell = 1$\n"
                 "por eso la (16) es una DDE neutra, no una EDO",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(-1.5, 9)
    ax.set_ylim(-0.15, 1.75)
    ax.legend(frameon=False, fontsize=8.5, loc="upper right")


def panel_trayectorias(ax, d):
    # Color = cutoff, estilo = m/m_B. Codificar las dos variables en el mismo
    # canal (cuatro azules) las volvia indistinguibles.
    estilos = {("libre", 0.4): (RAMPA[0], (0, (5, 2))),
               ("contraterm", 0.4): (RAMPA[0], "-"),
               ("libre", 0.8): (RAMPA[3], (0, (5, 2))),
               ("contraterm", 0.8): (RAMPA[3], "-")}
    for x in CUTOFFS:
        for etq in ("libre", "contraterm"):
            s, a, dr = d[f"s_{x}_{etq}"], d[f"a_{x}_{etq}"], d[f"d_{x}_{etq}"]
            mala = np.where(dr > DERIVA_MAX)[0]
            n = mala[0] if len(mala) else len(s)
            se, ae = envolvente(s[:n], a[:n])
            c, ls = estilos[(etq, x)]
            mr = 1.0 if etq == "libre" else float(d[f"mr_{x}"])
            ax.semilogy(se, np.maximum(ae, 1e-12), lw=1.7, color=c, ls=ls,
                        label=f"$r_0/\\ell$={x:.1f},  $m/m_B$={mr:.2f}")
            if n < len(s):
                ax.plot([se[-1]], [max(ae[-1], 1e-12)], "o", ms=5, color=c,
                        mec="#fcfcfb", mew=1.2, zorder=6)

    # Tasa espectral del cero dominante de chi, ec. (14), como referencia.
    z = complex(d[f"z_{CUTOFFS[1]}"][0])
    if np.isfinite(z.real):
        ss = np.linspace(9, 22, 50)
        s8, a8 = envolvente(d["s_0.8_contraterm"], d["a_0.8_contraterm"])
        i0 = np.searchsorted(s8, 9.0)
        ax.semilogy(ss, a8[i0] * np.exp(z.imag * (ss - 9.0)), lw=5,
                    color=RAMPA[3], alpha=.25, zorder=2, solid_capstyle="round")
        ax.annotate(f"banda: $e^{{\\,{z.imag:.3f}\\,s}}$, del cero de $\\chi$\n"
                    f"(ec. 14, mismo regulador)",
                    xy=(11.5, 8e-7), fontsize=8.5, color=INK)

    ax.annotate("● se corta donde la deriva de $\\dot{x}^2$ pasa de $10^{-3}$",
                xy=(0.6, 2e-9), fontsize=8.5, color=MUTED)
    ax.set_xlabel("$s / r_0$")
    ax.set_ylabel("envolvente de $|\\ddot{x}|\\; r_0$")
    ax.set_title("C.  El cutoff solo no decide: hacen falta los dos ejes\n"
                 f"sobre la línea del contratérmino el borde cae en "
                 f"$r_0/\\ell$ = {CRITICO_LINEAL:.3f}",
                 loc="left", fontsize=10.5, color=INK, pad=8)
    ax.set_xlim(0, S_END)
    ax.set_ylim(1e-10, 1e3)
    ax.legend(frameon=False, fontsize=8, loc="upper left")


def main():
    d = compute(force="--force" in sys.argv)

    fig, axes = plt.subplots(1, 3, figsize=(16.0, 5.0))
    panel_geometria(axes[0], d)
    panel_retardo(axes[1], d)
    panel_trayectorias(axes[2], d)

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
    axes[0].set_aspect("equal", adjustable="box")

    fig.suptitle("Ec. (16): regulador desplazado — el retardo finito y su "
                 "geometría invariante", x=0.008, ha="left",
                 fontsize=12.5, color=INK)
    fig.tight_layout(rect=(0, 0, 1, 0.93))
    salida = RAIZ / "figures/seccionA_ec16_retardo.png"
    fig.savefig(salida, dpi=170, facecolor="#fcfcfb")
    print(salida)


if __name__ == "__main__":
    main()
