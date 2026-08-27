"""
BLOQUE 3 -- Regulador suavizado, ec. (17): memoria infinita.

    xddot = (r0/3 ell^4)(m/m_B) int_{-inf}^{0} du (1 - sqrt(w)/2ell) e^{-sqrt(w)/ell} V

con  w = (x - x')^2,  x' = x(s + u),  y
     V = (x - x')(xdot . xdot') + [xdot . (x' - x)] xdot'.

Naturaleza: sistema no lineal de Volterra de SEGUNDA especie (la incognita
queda despejada fuera del operador causal), con nucleo dependiente de la
solucion. El nucleo es REGULAR en s' -> s: la exponencial tiende a 1 y V se
anula linealmente, de modo que no hay singularidad tipo Abel y la cuadratura
compuesta converge a orden pleno.

Truncamiento: el nucleo decae como exp(-sqrt(w)/ell) en distancia INVARIANTE.
`n_ell` fija la ventana en unidades de ell; con n_ell = 30 el peso residual
es O(e^-30) ~ 1e-13. `pts_per_ell` fija la resolucion de la cuadratura, de
modo que ambos errores se controlan por separado.
"""

from __future__ import annotations

import numpy as np

from .core import Params, minkowski_dot
from .worldline import WorldLine, rest_history, smooth_bump, unit_normal

__all__ = ["rhs_memory", "integrate_memory"]


def rhs_memory(x, v, wl: WorldLine, s_now: float, params: Params,
               n_ell: float = 30.0, pts_per_ell: int = 16):
    """Lado derecho de la ec. (17) por cuadratura de Simpson sobre la memoria.

    La cuadratura se parametriza por PUNTOS POR ell, no por numero total: asi
    el paso du = ell/pts_per_ell queda fijo al variar la ventana n_ell, y el
    error de truncamiento se puede estudiar sin contaminarlo con el error de
    cuadratura. (Con n_quad fijo, ampliar la ventana engrosaba el paso y los
    dos errores se mezclaban.)
    """
    ell = params.ell
    window = n_ell * ell
    n = int(np.ceil(n_ell * pts_per_ell))
    n = n + 1 if n % 2 == 0 else n                       # Simpson necesita impar
    sp = np.linspace(s_now - window, s_now, n)
    du = sp[1] - sp[0]

    xp, vp, _ = wl.sample_many(sp)

    d = x - xp                                            # (x - x'), (n, dim)
    w = minkowski_dot(d, d)
    w = np.maximum(w, 0.0)                                # separacion temporal
    sq = np.sqrt(w)

    kernel = (1.0 - sq / (2.0 * ell)) * np.exp(-sq / ell)

    v_vp = minkowski_dot(v, vp)
    v_d = minkowski_dot(v, d)
    V = d * v_vp[:, None] - v_d[:, None] * vp             # [xdot.(x'-x)] = -v_d

    integrand = kernel[:, None] * V

    # Simpson compuesta
    wts = np.ones(n); wts[1:-1:2] = 4.0; wts[2:-1:2] = 2.0
    integral = (du / 3.0) * np.einsum("i,ij->j", wts, integrand)

    return params.r0 * params.m_over_mB / (3.0 * ell ** 4) * integral


def integrate_memory(
    params: Params,
    s_end: float = 30.0,
    ds: float = 2e-3,
    history: WorldLine | None = None,
    n_corrector: int = 2,
    n_ell: float = 30.0,
    pts_per_ell: int = 16,
    drive=None,
):
    """Integra la ec. (17) hacia adelante en tiempo propio.

    Misma marcha explicita que el bloque 2 -- la estructura de segunda especie
    lo permite -- con un predictor-corrector para el acoplamiento entre x(s) y
    el nucleo, que depende de x(s) a traves de w = (x - x')^2.

    La prehistoria debe cubrir al menos n_ell*ell hacia atras; `rest_then_kick`
    lo garantiza con su segmento en reposo.
    """
    ell = params.ell
    if drive is None:
        drive = lambda s: smooth_bump(s, amplitude=0.3, s0=0.0, width=1.0)

    wl = history if history is not None else rest_history(
        dim=params.dim, ds=ds, s_rest=max(4.0, 1.5 * n_ell * ell)
    )
    if wl.s_max - wl.s_min < n_ell * ell:
        raise ValueError("la prehistoria es mas corta que la ventana de memoria")

    n_steps = int(np.ceil(s_end / ds))
    for _ in range(n_steps):
        s_n = wl.s_max
        x_n, v_n, a_n = wl._x[-1], wl._v[-1], wl._a[-1]
        s_new = s_n + ds

        x_p = x_n + ds * v_n + 0.5 * ds ** 2 * a_n
        v_p = v_n + ds * a_n

        a_new = a_n
        for _ in range(n_corrector):
            a_new = rhs_memory(x_p, v_p, wl, s_new, params,
                               n_ell=n_ell, pts_per_ell=pts_per_ell)
            k = drive(s_new)
            if k:
                a_new = a_new + k * unit_normal(v_p)
            v_p = v_n + 0.5 * ds * (a_n + a_new)
            x_p = x_n + 0.5 * ds * (v_n + v_p)

        if not np.all(np.isfinite(a_new)):
            raise FloatingPointError(f"divergencia en s={s_new:.4f}")
        wl.append(s_new, x_p, v_p, a_new)

    return wl
