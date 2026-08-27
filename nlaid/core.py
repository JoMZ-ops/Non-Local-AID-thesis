"""
Nucleo compartido: metrica, parametros y reguladores point-splitting.

Todos los bloques (14), (16) y (17) importan de aqui. La idea es que exista
UNA sola definicion de cada regulador y UN solo objeto de parametros, para que
cambiar el regulador o el cutoff se propague a los tres analisis.

Referencia: J. Polonyi, "The Abraham-Lorentz force and electrodynamics at the
classical electron radius", arXiv:1701.04068v4.

Convenciones
------------
Signatura        (+,-,-,-)   -> lineas de mundo temporales tienen x^2 > 0
Unidades         c = 1, r0 = 1.  El cutoff ell se mide en unidades de r0,
                 de modo que el eje de la Fig. 2 del paper es r0/ell = 1/ell.
Parametrizacion  tiempo propio s, con la normalizacion xdot^2 = 1.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

__all__ = [
    "minkowski_dot",
    "Params",
    "Regulator",
    "ShiftedDelta",
    "SmearedDelta",
    "make_regulator",
]


# --------------------------------------------------------------------------
# Metrica
# --------------------------------------------------------------------------

def minkowski_dot(a, b):
    """Producto (+,-,-,-) sobre el ultimo eje. Difunde sobre los ejes previos."""
    a = np.asarray(a)
    b = np.asarray(b)
    return a[..., 0] * b[..., 0] - np.sum(a[..., 1:] * b[..., 1:], axis=-1)


# --------------------------------------------------------------------------
# Parametros
# --------------------------------------------------------------------------

@dataclass(frozen=True)
class Params:
    """Parametros fisicos compartidos por los tres bloques.

    ell         cutoff de point-splitting, en unidades de r0.
    r0          radio clasico del electron (=1 en unidades naturales).
    m_over_mB   m/m_B, el segundo eje del diagrama de fases (Fig. 2).
    dim         dimension del espacio-tiempo. dim=2 es el caso 1+1 usado
                para reproducir las figuras del paper.
    """

    ell: float
    r0: float = 1.0
    m_over_mB: float = 1.0
    dim: int = 2

    def __post_init__(self):
        if self.ell <= 0:
            raise ValueError(f"ell debe ser > 0, se recibio {self.ell}")
        if self.r0 <= 0:
            raise ValueError(f"r0 debe ser > 0, se recibio {self.r0}")
        if self.dim < 2:
            raise ValueError(f"dim debe ser >= 2, se recibio {self.dim}")

    @property
    def r0_over_ell(self) -> float:
        """Abscisa del diagrama de fases de la Fig. 2."""
        return self.r0 / self.ell


# --------------------------------------------------------------------------
# Reguladores
# --------------------------------------------------------------------------

class Regulator:
    """Interfaz de un regulador delta_B(z), z = x^2 el intervalo invariante.

    Las tres condiciones que exige el paper (p. 6) son:
      (i)   int dz delta_B(z) = 1          -> preserva el flujo radiado
      (ii)  delta_B(0) = 0                 -> separa los puntos singulares
      (iii) delta_B(z) = 0 para z < 0      -> suprime la interaccion superluminica

    `check_conditions` las verifica numericamente para cualquier subclase.
    """

    name = "abstract"

    def delta(self, z):
        """delta_B(z)."""
        raise NotImplementedError

    def d_delta(self, z):
        """delta_B'(z) = d delta_B / dz. La necesita la ec. (6)."""
        raise NotImplementedError

    def integrate_u(self, g):
        """int_{-inf}^{0} du delta_B(u^2) g(u), con g vectorizada sobre u.

        Este es el funcional que aparece en la teoria linealizada, donde
        (x - x')^2 -> u^2. Cada regulador lo implementa con la cuadratura que
        le corresponde: exacta para la delta desplazada, Gauss-Laguerre para
        la suavizada.
        """
        raise NotImplementedError

    # -- diagnosticos compartidos ------------------------------------------

    def mass_shift_over_m(self, r0: float = 1.0) -> float:
        """delta_m / m a partir de la ec. (10).

        delta_m = (e^2 / 2 c^2) int_0^inf dz z^{-1/2} delta_B(z), y como
        e^2/c^2 = m r0, resulta delta_m/m = (r0/2) int_0^inf dz z^{-1/2} delta_B(z).
        """
        return 0.5 * r0 * self.moment_inv_sqrt()

    def moment_inv_sqrt(self) -> float:
        """int_0^inf dz z^{-1/2} delta_B(z)."""
        raise NotImplementedError

    def m_over_mB_counterterm(self, r0: float = 1.0) -> float:
        """m/m_B fijado por el contraterm de la ec. (10).

        Como m = m_B + delta_m, se tiene m/m_B = 1/(1 - delta_m/m). Es la
        LINEA PUNTEADA de la Fig. 2 del paper: la prediccion de la teoria
        linealizada en el plano (r0/ell, m/m_B). Diverge cuando delta_m = m
        (r0/ell = 6 suavizado, r0/ell = 2 desplazado) y pasa a m_B < 0.
        """
        return 1.0 / (1.0 - self.mass_shift_over_m(r0))

    def moment_u2(self) -> float:
        """I2 = int_{-inf}^{0} du delta_B(u^2) u^2.

        Fija la correccion O(omega^2) de la susceptibilidad:
        chi = 1 + r0 [ (2/3) i omega - (3/4) omega^2 I2 ] + O(omega^3).
        """
        return self.integrate_u(lambda u: u ** 2)

    def check_conditions(self, atol: float = 1e-10) -> dict:
        """Verifica numericamente las tres condiciones del regulador."""
        from scipy.integrate import quad

        norm = quad(self.delta, 0.0, np.inf, limit=400)[0]
        return {
            "normalizacion": norm,
            "error_normalizacion": abs(norm - 1.0),
            "delta_en_cero": float(self.delta(0.0)),
            "soporte_negativo": float(self.delta(-1.0)),
            "cumple": (
                abs(norm - 1.0) < 1e-8
                and abs(float(self.delta(0.0))) < atol
                and abs(float(self.delta(-1.0))) < atol
            ),
        }


class ShiftedDelta(Regulator):
    """Ec. (4): delta_B(x^2) = delta(x^2 - ell^2).

    Desplaza la funcion de Green retardada fuera del cono de luz. Es la
    version mas simple; produce oscilaciones en espacio de momentos, pero
    da forma cerrada para todos los funcionales, lo que la vuelve el test
    de referencia del regulador suavizado.

    Sobre u < 0 la identidad delta(u^2 - ell^2) = [delta(u+ell) + delta(u-ell)]/(2 ell)
    deja solo la raiz u = -ell, de modo que
        int_{-inf}^{0} du delta_B(u^2) g(u) = g(-ell) / (2 ell).
    """

    name = "shifted"

    def __init__(self, ell: float):
        self.ell = float(ell)

    def delta(self, z):
        # Distribucion: sin densidad puntual. Se representa por su accion.
        return np.zeros_like(np.asarray(z, dtype=float))

    def d_delta(self, z):
        return np.zeros_like(np.asarray(z, dtype=float))

    def integrate_u(self, g):
        return g(-self.ell) / (2.0 * self.ell)

    def moment_inv_sqrt(self) -> float:
        # int_0^inf dz z^{-1/2} delta(z - ell^2) = 1/ell
        return 1.0 / self.ell

    def check_conditions(self, atol: float = 1e-10) -> dict:
        # La delta desplazada satisface las tres condiciones por construccion,
        # pero como distribucion no admite la cuadratura de la clase base.
        return {
            "normalizacion": 1.0,
            "error_normalizacion": 0.0,
            "delta_en_cero": 0.0,
            "soporte_negativo": 0.0,
            "cumple": True,
            "nota": "distribucion: condiciones exactas por construccion",
        }


class SmearedDelta(Regulator):
    """Ec. (5): delta_B(x^2) = Theta(x^2) x^2 exp(-sqrt(x^2)/ell) / (12 ell^4).

    Suaviza la singularidad para evitar las oscilaciones en espacio de
    momentos, al costo de memoria infinita. Es el regulador con el que
    Polonyi produce las Figs. 1 y 3.

    La integral sobre u se hace con Gauss-Laguerre generalizada (alpha=2):
    sustituyendo u = -ell t,
        int_{-inf}^{0} du delta_B(u^2) g(u)
            = (1 / 12 ell) int_0^inf dt t^2 e^{-t} g(-ell t),
    que la cuadratura integra con precision espectral.
    """

    name = "smeared"

    def __init__(self, ell: float, n_quad: int = 120):
        from scipy.special import roots_genlaguerre

        self.ell = float(ell)
        self.n_quad = int(n_quad)
        self._t, self._w = roots_genlaguerre(self.n_quad, 2.0)

    def delta(self, z):
        z = np.asarray(z, dtype=float)
        out = np.zeros_like(z)
        pos = z > 0
        zp = z[pos] if z.ndim else (z if pos else None)
        if z.ndim == 0:
            if pos:
                return np.asarray(
                    z * np.exp(-np.sqrt(z) / self.ell) / (12.0 * self.ell ** 4)
                )
            return out
        out[pos] = zp * np.exp(-np.sqrt(zp) / self.ell) / (12.0 * self.ell ** 4)
        return out

    def d_delta(self, z):
        """delta_B'(z) = [1 - sqrt(z)/(2 ell)] exp(-sqrt(z)/ell) / (12 ell^4).

        Es el factor que aparece explicitamente en la ec. (17).
        """
        z = np.asarray(z, dtype=float)
        out = np.zeros_like(z)
        if z.ndim == 0:
            if z > 0:
                sq = np.sqrt(z)
                return np.asarray(
                    (1.0 - sq / (2.0 * self.ell))
                    * np.exp(-sq / self.ell)
                    / (12.0 * self.ell ** 4)
                )
            return out
        pos = z > 0
        sq = np.sqrt(z[pos])
        out[pos] = (
            (1.0 - sq / (2.0 * self.ell))
            * np.exp(-sq / self.ell)
            / (12.0 * self.ell ** 4)
        )
        return out

    def integrate_u(self, g):
        vals = g(-self.ell * self._t)
        return np.sum(self._w * vals) / (12.0 * self.ell)

    def moment_inv_sqrt(self) -> float:
        # int_0^inf dz z^{-1/2} delta_B(z) = 1/(3 ell), analitico
        return 1.0 / (3.0 * self.ell)


def make_regulator(kind: str, ell: float, **kwargs) -> Regulator:
    """Fabrica: 'shifted' -> ec. (4), 'smeared' -> ec. (5)."""
    kinds = {"shifted": ShiftedDelta, "smeared": SmearedDelta}
    if kind not in kinds:
        raise ValueError(f"regulador desconocido {kind!r}; use uno de {sorted(kinds)}")
    return kinds[kind](ell, **kwargs)
