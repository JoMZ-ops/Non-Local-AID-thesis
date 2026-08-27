"""
Historia de la linea de mundo: almacenamiento denso, interpolacion y punto retardado.

Pieza COMPARTIDA por los bloques 2 (ec. 16) y 3 (ec. 17). Ambos necesitan
exactamente lo mismo -- consultar x, xdot, xddot en un tiempo propio pasado
arbitrario -- y ambos arrancan de la misma prehistoria prescrita del paper
(p. 9): carga en reposo, luego una trayectoria impuesta, luego evolucion libre.

La normalizacion xdot^2 = 1 NO se impone: los lados derechos de las ecs. (6),
(16) y (17) son ortogonales a xdot de forma exacta (se verifica contrayendo
con xdot_mu), asi que la deriva de xdot^2 mide directamente el error de
integracion. `WorldLine.norm_drift` la reporta.

Condicion inicial
-----------------
Una carga en reposo es solucion EXACTA de las tres ecuaciones: con
xdot = xdot' = (1,0) y (x-x') proporcional a xdot, el vector
    V1 = (x-x')(xdot.xdot') - [xdot.(x-x')] xdot'
se anula identicamente, es decir no hay autofuerza para movimiento inercial,
como debe ser. Por eso la prehistoria inercial NO introduce ningun salto.

La excitacion se aplica mediante la fuente externa k^mu de la ec. (2), como
un pulso suave de soporte compacto. Prescribir una trayectoria y apagar la
fuente de golpe haria saltar xddot en s=0; en una ecuacion con retardo esa
discontinuidad se propaga a s = ell, 2 ell, ... y degrada el orden del
integrador. Un pulso C^infinito evita el problema por completo.
"""

from __future__ import annotations

import numpy as np
from scipy.optimize import brentq

from .core import minkowski_dot

__all__ = ["WorldLine", "rest_history", "unit_normal", "smooth_bump"]


class WorldLine:
    """Muestreo denso de (s, x, xdot, xddot), con interpolacion cubica.

    Se guarda como buffer creciente. La interpolacion es de Hermite cubica a
    trozos usando x y xdot (que son consistentes por construccion), lo que da
    O(ds^4) en x y O(ds^3) en xdot -- suficiente para los integradores de
    segundo orden de los bloques 2 y 3, y mas estable que un spline global.
    """

    def __init__(self, dim: int):
        self.dim = int(dim)
        self._s: list[float] = []
        self._x: list[np.ndarray] = []
        self._v: list[np.ndarray] = []
        self._a: list[np.ndarray] = []

    # -- construccion -------------------------------------------------------

    def append(self, s: float, x, v, a) -> None:
        if self._s and s <= self._s[-1]:
            raise ValueError(f"s debe crecer: {s} <= {self._s[-1]}")
        self._s.append(float(s))
        self._x.append(np.asarray(x, dtype=float).copy())
        self._v.append(np.asarray(v, dtype=float).copy())
        self._a.append(np.asarray(a, dtype=float).copy())

    def __len__(self) -> int:
        return len(self._s)

    @property
    def s(self):
        return np.asarray(self._s)

    @property
    def x(self):
        return np.asarray(self._x)

    @property
    def v(self):
        return np.asarray(self._v)

    @property
    def a(self):
        return np.asarray(self._a)

    @property
    def s_min(self) -> float:
        return self._s[0]

    @property
    def s_max(self) -> float:
        return self._s[-1]

    # -- consulta -----------------------------------------------------------

    def sample(self, sq: float):
        """(x, xdot, xddot) en el tiempo propio sq, por Hermite cubica.

        Antes de s_min se extrapola con el segmento en reposo (movimiento
        inercial), que es la prehistoria exacta del problema.
        """
        s = self.s
        if sq < s[0]:
            # prehistoria inercial: xdot constante, xddot = 0
            v0, x0 = self._v[0], self._x[0]
            return x0 + (sq - s[0]) * v0, v0.copy(), np.zeros(self.dim)
        if sq >= s[-1]:
            return self._x[-1].copy(), self._v[-1].copy(), self._a[-1].copy()

        j = int(np.searchsorted(s, sq, side="right") - 1)
        j = min(max(j, 0), len(s) - 2)
        h = s[j + 1] - s[j]
        t = (sq - s[j]) / h

        x0, x1 = self._x[j], self._x[j + 1]
        v0, v1 = self._v[j], self._v[j + 1]
        a0, a1 = self._a[j], self._a[j + 1]

        t2, t3 = t * t, t * t * t
        # Hermite en x usando xdot como tangente
        h00 = 2 * t3 - 3 * t2 + 1
        h10 = t3 - 2 * t2 + t
        h01 = -2 * t3 + 3 * t2
        h11 = t3 - t2
        xq = h00 * x0 + h10 * h * v0 + h01 * x1 + h11 * h * v1
        # Hermite en xdot usando xddot como tangente
        vq = h00 * v0 + h10 * h * a0 + h01 * v1 + h11 * h * a1
        # aceleracion: lineal (basta al orden de los integradores)
        aq = (1 - t) * a0 + t * a1
        return xq, vq, aq

    # -- punto retardado ----------------------------------------------------

    def retarded_point(self, s_now: float, x_now, ell: float, s_hint=None):
        """Resuelve (x_now - x(s'))^2 = ell^2 con s' < s_now.

        Sobre una linea de mundo temporal el intervalo invariante crece
        monotonamente al retroceder, de modo que la raiz retardada es unica.
        Se busca por biseccion de Brent sobre un intervalo acotado.
        """
        x_now = np.asarray(x_now, dtype=float)

        def gap(sp: float) -> float:
            xp = self.sample(sp)[0]
            d = x_now - xp
            return minkowski_dot(d, d) - ell ** 2

        # cota superior: justo antes de s_now, gap < 0
        hi = s_now - 1e-12
        # cota inferior: retroceder hasta que gap > 0
        span = max(ell, 1e-9)
        lo = s_now - span
        for _ in range(200):
            if gap(lo) > 0:
                break
            span *= 1.7
            lo = s_now - span
        else:
            raise RuntimeError(
                "no se encontro el punto retardado; la historia es demasiado corta"
            )

        sp = brentq(gap, lo, hi, xtol=1e-14, rtol=1e-14, maxiter=300)
        xp, vp, ap = self.sample(sp)
        return sp, xp, vp, ap

    # -- diagnostico --------------------------------------------------------

    @property
    def norm_drift(self):
        """|xdot^2 - 1| a lo largo de la historia.

        Los lados derechos de (6), (16) y (17) son exactamente ortogonales a
        xdot, de modo que xdot^2 se conserva analiticamente. Su deriva es una
        medida directa y gratuita del error de integracion.
        """
        v = self.v
        return np.abs(minkowski_dot(v, v) - 1.0)


def unit_normal(v, direction: int = 1):
    """Vector unitario espacial ortogonal a la 4-velocidad v.

    Construido proyectando el eje espacial `direction` ortogonalmente a v y
    normalizando:  n = p / sqrt(-p.p)  con  p = e - (e.v) v.
    Cumple n.v = 0 y n.n = -1, de modo que cualquier aceleracion externa
    proporcional a n preserva xdot^2 = 1 exactamente.
    """
    v = np.asarray(v, dtype=float)
    e = np.zeros_like(v)
    e[direction] = 1.0
    p = e - minkowski_dot(e, v) * v
    return p / np.sqrt(-minkowski_dot(p, p))


def smooth_bump(s, amplitude: float = 0.3, s0: float = 0.0, width: float = 1.0):
    """Pulso C^infinito de soporte compacto en [s0, s0 + width].

    exp(-1/(u(1-u))) reescalado a maximo `amplitude`. Todas sus derivadas se
    anulan en los extremos, de modo que xddot queda infinitamente suave y el
    integrador conserva su orden nominal.
    """
    u = (s - s0) / width
    if u <= 0.0 or u >= 1.0:
        return 0.0
    return amplitude * np.exp(1.0 - 1.0 / (4.0 * u * (1.0 - u)))


def rest_history(dim: int = 2, s_rest: float = 60.0, ds: float = 1e-3):
    """Prehistoria inercial exacta: carga en reposo, x = (s, 0), xddot = 0.

    Es solucion exacta de las ecs. (16) y (17), asi que no aporta ningun
    transitorio espurio. Cubre s in [-s_rest, 0]; debe ser mas larga que la
    ventana de memoria del bloque 3 o que el retardo del bloque 2.
    """
    if dim < 2:
        raise ValueError("dim debe ser >= 2")
    wl = WorldLine(dim)
    v = np.zeros(dim); v[0] = 1.0
    a = np.zeros(dim)
    for sv in np.arange(-s_rest, 0.0 + 0.5 * ds, ds):
        x = np.zeros(dim); x[0] = sv
        wl.append(float(sv), x, v, a)
    return wl


def _hermite(wl: "WorldLine", sq: np.ndarray):
    """Version vectorizada de WorldLine.sample sobre un arreglo de tiempos propios."""
    sq = np.atleast_1d(np.asarray(sq, dtype=float))
    s, X, V, A = wl.s, wl.x, wl.v, wl.a

    xq = np.empty((sq.size, wl.dim))
    vq = np.empty_like(xq)
    aq = np.empty_like(xq)

    before = sq < s[0]
    if np.any(before):
        xq[before] = X[0] + (sq[before, None] - s[0]) * V[0]
        vq[before] = V[0]
        aq[before] = 0.0

    after = sq >= s[-1]
    if np.any(after):
        xq[after], vq[after], aq[after] = X[-1], V[-1], A[-1]

    mid = ~(before | after)
    if np.any(mid):
        sm = sq[mid]
        j = np.clip(np.searchsorted(s, sm, side="right") - 1, 0, len(s) - 2)
        h = (s[j + 1] - s[j])[:, None]
        t = ((sm - s[j]) / (s[j + 1] - s[j]))[:, None]
        t2, t3 = t * t, t * t * t
        h00, h10 = 2 * t3 - 3 * t2 + 1, t3 - 2 * t2 + t
        h01, h11 = -2 * t3 + 3 * t2, t3 - t2
        xq[mid] = h00 * X[j] + h10 * h * V[j] + h01 * X[j + 1] + h11 * h * V[j + 1]
        vq[mid] = h00 * V[j] + h10 * h * A[j] + h01 * V[j + 1] + h11 * h * A[j + 1]
        aq[mid] = (1 - t) * A[j] + t * A[j + 1]

    return xq, vq, aq


WorldLine.sample_many = _hermite
