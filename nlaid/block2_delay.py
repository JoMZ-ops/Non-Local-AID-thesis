"""
BLOQUE 2 -- Regulador desplazado, ec. (16): retardo finito.

    xddot = r0 (m/m_B) / A^2 [ (B + 1)/A * V1 + V2 ]

con  A  = xdot' . (x - x')
     B  = xddot' . (x' - x)
     V1 = (x - x')(xdot . xdot') - [xdot . (x - x')] xdot'
     V2 = (x - x')(xdot . xddot') + [xdot . (x' - x)] xddot'

y el punto retardado x' fijado por  ell^2 = (x - x')^2.

Naturaleza: ecuacion diferencial con retardo de tipo NEUTRO y retardo
dependiente del estado. xddot' es historia ya calculada (el retardo esta
acotado por debajo por ~ell), pero la UBICACION del punto retardado depende
de x(s), que es la incognita del paso. De ahi el predictor-corrector.
"""

from __future__ import annotations

import numpy as np

from .core import Params, minkowski_dot
from .worldline import WorldLine, rest_history, smooth_bump, unit_normal

BLOWUP = 1e6      # |xddot| por encima del cual la trayectoria ya escapo

__all__ = ["rhs_delay", "integrate_delay"]


def rhs_delay(x, v, xp, vp, ap, params: Params):
    """Lado derecho de la ec. (16). Devuelve (xddot, diagnosticos)."""
    d = x - xp                                   # (x - x')
    A = minkowski_dot(vp, d)                     # xdot' . (x - x')
    if not np.isfinite(A) or abs(A) < 1e-14:
        raise FloatingPointError(f"denominador degenerado en la ec. (16): A={A}")
    B = minkowski_dot(ap, -d)                    # xddot' . (x' - x)

    v_vp = minkowski_dot(v, vp)                  # xdot . xdot'
    v_d = minkowski_dot(v, d)                    # xdot . (x - x')
    v_ap = minkowski_dot(v, ap)                  # xdot . xddot'

    V1 = d * v_vp - v_d * vp
    V2 = d * v_ap - v_d * ap                     # [xdot.(x'-x)] = -v_d

    acc = params.r0 * params.m_over_mB / A ** 2 * ((B + 1.0) / A * V1 + V2)
    return acc, {"A": A, "B": B}


def integrate_delay(
    params: Params,
    s_end: float = 30.0,
    ds: float = 2e-3,
    history: WorldLine | None = None,
    n_corrector: int = 2,
    drive=None,
):
    """Integra la ec. (16) hacia adelante en tiempo propio.

    Predictor-corrector: predice (x, xdot), relocaliza el punto retardado con
    la prediccion, evalua el lado derecho y corrige con la regla trapezoidal.
    `n_corrector` iteraciones bastan porque el retardo esta acotado por debajo
    (~ell), de modo que la implicitud es debil.

    Requiere ds << ell: el punto retardado debe caer dentro de la historia ya
    registrada, no en el paso que se esta calculando.

    `drive(s)` es la aceleracion propia externa -- la fuente k^mu de la ec. (2)
    que "diagnostica la dinamica". Se aplica a lo largo del normal unitario a
    xdot, de modo que preserva xdot^2 = 1. Por defecto un pulso suave de
    soporte compacto en s in [0, 1].
    """
    ell = params.ell
    if ds > 0.2 * ell:
        raise ValueError(
            f"ds={ds} demasiado grande frente a ell={ell}; use ds <= {0.2*ell:.2e}"
        )

    if drive is None:
        drive = lambda s: smooth_bump(s, amplitude=0.3, s0=0.0, width=1.0)

    wl = history if history is not None else rest_history(
        dim=params.dim, ds=ds, s_rest=max(20.0 * ell, 4.0)
    )

    n_steps = int(np.ceil(s_end / ds))
    for _ in range(n_steps):
        s_n = wl.s_max
        x_n, v_n, a_n = wl._x[-1], wl._v[-1], wl._a[-1]
        s_new = s_n + ds

        # predictor
        x_p = x_n + ds * v_n + 0.5 * ds ** 2 * a_n
        v_p = v_n + ds * a_n

        a_new = a_n
        for _ in range(n_corrector):
            _, xp, vp, ap = wl.retarded_point(s_new, x_p, ell)
            a_new, _ = rhs_delay(x_p, v_p, xp, vp, ap, params)
            k = drive(s_new)
            if k:
                a_new = a_new + k * unit_normal(v_p)
            v_p = v_n + 0.5 * ds * (a_n + a_new)
            x_p = x_n + 0.5 * ds * (v_n + v_p)

        if not np.all(np.isfinite(a_new)):
            raise FloatingPointError(f"divergencia en s={s_new:.4f}")
        wl.append(s_new, x_p, v_p, a_new)

        # Corte temprano por escape. Una vez que |xddot| supera BLOWUP la
        # trayectoria ya escapo y seguir integrando no anade informacion:
        # solo acumula error y, en la ec. (16), dispara la busqueda del punto
        # retardado (que debe retroceder cada vez mas para alcanzar el cono).
        # Se devuelve la historia truncada; los diagnosticos de amplitud y
        # deriva bastan para clasificarla como inestable.
        if np.linalg.norm(a_new) > BLOWUP:
            break

    return wl
