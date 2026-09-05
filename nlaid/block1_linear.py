"""
BLOQUE 1 -- Teoria linealizada: susceptibilidad chi^r_omega, ec. (14).

Objetivo: decidir estabilidad y causalidad del regulador SIN integrar ninguna
trayectoria. La dinamica armonica es estable y causal si chi^r_omega es
analitica y libre de polos en el semiplano superior (p. 8 del paper); como
F^r_omega = 1 / [(omega + i eps)^2 chi^r_omega], los polos de F^r son los
CEROS de chi^r. Buscamos entonces ceros de chi^r en Im omega > 0.

Ec. (14):
    chi^r_omega = 1 + r0 [ (2/3) i omega
                           - (2/omega^2) int_{-inf}^{0} du delta_B(u^2) N(omega u)/u^2 ]
con
    N(phi) = (1 + i phi - phi^2) e^{-i phi} - 1 + phi^2/2 - (2/3) i phi^3.

Los cuatro primeros ordenes de N se cancelan identicamente, dejando
    N(phi) = sum_{n>=4} (-i phi)^n (n-1)^2 / n!,
serie que usamos a |phi| pequeno para evitar cancelacion catastrofica.
"""

from __future__ import annotations

import numpy as np

from .core import Params, Regulator

__all__ = [
    "numerator_N",
    "dispersion",
    "integrate_linear",
    "susceptibility",
    "susceptibility_small_omega",
    "count_zeros_uhp",
    "find_zeros_uhp",
    "is_stable",
    "convergence_lower_bound",
    "dominant_pole",
    "critical_cutoff",
]

# Escala de busqueda por defecto, en unidades r0=1. El polo runaway de
# Abraham-Lorentz de la ec. (15) esta en omega = 3i/2, de modo que una caja
# hasta 12i lo cubre con holgura junto con los ceros que aporta el cutoff.
_DEFAULT_IM_HI = 12.0
_DEFAULT_RE_MAX = 40.0

_SERIES_CUTOFF = 6.0   # usar la serie para |phi| < esto
_SERIES_TERMS = 64


def numerator_N(phi):
    """N(phi), estable en todo el plano complejo.

    Serie para |phi| < _SERIES_CUTOFF (evita la cancelacion de los cuatro
    ordenes que se anulan), forma directa fuera.
    """
    phi = np.asarray(phi, dtype=complex)
    scalar = phi.ndim == 0
    phi = np.atleast_1d(phi)
    out = np.empty_like(phi)

    small = np.abs(phi) < _SERIES_CUTOFF

    if np.any(small):
        ps = phi[small]
        x = -1j * ps                      # a_n = x^n / n!
        acc = np.zeros_like(ps)
        a = np.ones_like(ps)              # a_0
        for n in range(1, _SERIES_TERMS):
            a = a * x / n
            if n >= 4:
                acc += a * (n - 1) ** 2
        out[small] = acc

    if np.any(~small):
        pl = phi[~small]
        with np.errstate(over="ignore", invalid="ignore"):
            out[~small] = (
                (1.0 + 1j * pl - pl ** 2) * np.exp(-1j * pl)
                - 1.0
                + 0.5 * pl ** 2
                - (2.0 / 3.0) * 1j * pl ** 3
            )

    return out[0] if scalar else out


def susceptibility(omega, reg: Regulator, params: Params):
    """chi^r_omega segun la ec. (14). Vectorizada sobre omega complejo.

    Validez: sustituyendo u = -v (v > 0), el integrando decae como
    v^2 exp[-v (1/ell + Im omega)], de modo que la representacion converge para
    Im(omega) > -1/ell. Cubre TODO el semiplano superior. (El paper pide
    Im omega > 0, que es condicion suficiente pero no el limite real; ver
    `convergence_lower_bound`.)
    """
    omega = np.asarray(omega, dtype=complex)
    scalar = omega.ndim == 0
    omega = np.atleast_1d(omega)

    # int_{-inf}^{0} du delta_B(u^2) N(w u) / u^2, para TODOS los omega a la vez.
    # La cuadratura se pide al regulador como nodos y pesos, de modo que el
    # integrando se evalua sobre la malla completa (n_omega x n_nodos) en una
    # sola llamada. La version escalar -- un bucle de Python sobre omega, con
    # una suma de 120 nodos de Laguerre dentro -- dominaba el costo de todo el
    # bloque 1: las busquedas de ceros recorren decenas de miles de omega.
    u, wq = reg.nodes_weights()
    phi = omega[:, None] * u[None, :]
    # errstate acotado, no global: a |omega| enorme numerator_N desborda a inf
    # y aqui sale inf/inf. Ocurre solo cuando Muller tantea lejos, y la caja
    # acotante descarta esos iterados; el resultado no se usa.
    with np.errstate(invalid="ignore"):
        integ = np.einsum("j,ij->i", wq, numerator_N(phi) / u[None, :] ** 2)

    with np.errstate(divide="ignore", invalid="ignore"):
        bracket = (2.0 / 3.0) * 1j * omega - 2.0 * integ / omega ** 2

    # omega -> 0: chi -> 1 (el termino 2 I/omega^2 tiende a (3/4) omega^2 I2 -> 0)
    bracket = np.where(np.abs(omega) < 1e-12, 0.0 + 0.0j, bracket)

    chi = 1.0 + params.r0 * bracket
    return chi[0] if scalar else chi


def susceptibility_small_omega(omega, reg: Regulator, params: Params):
    """Desarrollo chi = 1 + r0 [ (2/3) i omega - (3/4) omega^2 I2 ] + O(omega^3).

    Sirve como test independiente de `susceptibility` a omega pequeno, y
    reproduce la ec. (15) al orden dominante.
    """
    omega = np.asarray(omega, dtype=complex)
    i2 = reg.moment_u2()
    return 1.0 + params.r0 * ((2.0 / 3.0) * 1j * omega - 0.75 * omega ** 2 * i2)


def convergence_lower_bound(reg: Regulator) -> float:
    """Cota INFERIOR de Im(omega) para la que converge la representacion integral.

    Suavizado: -1/ell.  Desplazado: -inf (la cuadratura es exacta).
    En ambos casos el semiplano superior completo es accesible.
    """
    if reg.name == "shifted":
        return -np.inf
    ell = getattr(reg, "ell", None)
    return -1.0 / ell if ell else -np.inf


# --------------------------------------------------------------------------
# Busqueda de ceros en el semiplano superior
# --------------------------------------------------------------------------

def _rectangle(re_max: float, im_lo: float, im_hi: float, n: int):
    """Contorno rectangular cerrado en el semiplano superior, sentido antihorario."""
    a = np.linspace(-re_max, re_max, n) + 1j * im_lo
    b = re_max + 1j * np.linspace(im_lo, im_hi, n)
    c = np.linspace(re_max, -re_max, n) + 1j * im_hi
    d = -re_max + 1j * np.linspace(im_hi, im_lo, n)
    return np.concatenate([a, b, c, d, a[:1]])


def count_zeros_uhp(
    reg: Regulator,
    params: Params,
    re_max: float = _DEFAULT_RE_MAX,
    im_lo: float = 1e-6,
    im_hi: float = _DEFAULT_IM_HI,
    n: int = 400,
) -> int:
    """Numero de ceros de chi^r dentro del rectangulo, por principio del argumento.

    El conteo es el numero de vueltas de chi alrededor del origen al recorrer
    el contorno. Es un diagnostico global: no depende de una semilla.
    """
    contour = _rectangle(re_max, im_lo, im_hi, n)
    vals = susceptibility(contour, reg, params)
    if np.any(np.abs(vals) < 1e-14):
        raise RuntimeError("el contorno pasa por un cero de chi; mueva im_lo o re_max")

    winding = np.sum(np.diff(np.unwrap(np.angle(vals)))) / (2.0 * np.pi)
    return int(np.rint(winding))


def _muller(f, x0, x1, x2, tol=1e-13, maxiter=80, box=None):
    """Muller: root-finder complejo sin derivadas. Converge cuadraticamente.

    `box = (re_max, im_lo, im_hi)` aborta si el iterado se sale de la region,
    lo que evita que la iteracion se dispare a |omega| enorme donde el
    integrando desborda.
    """
    def outside(z):
        if box is None or not np.isfinite(z):
            return not np.isfinite(z)
        re_max, im_lo, im_hi = box
        return abs(z.real) > 2 * re_max or z.imag < im_lo - 1 or z.imag > 2 * im_hi

    f0, f1, f2 = f(x0), f(x1), f(x2)
    for _ in range(maxiter):
        q = (x2 - x1) / (x1 - x0)
        a = q * f2 - q * (1 + q) * f1 + q ** 2 * f0
        b = (2 * q + 1) * f2 - (1 + q) ** 2 * f1 + q ** 2 * f0
        c = (1 + q) * f2
        disc = np.sqrt(b ** 2 - 4 * a * c + 0j)
        den = b + disc if abs(b + disc) > abs(b - disc) else b - disc
        if den == 0:
            return None
        x3 = x2 - (x2 - x1) * 2 * c / den
        if outside(x3):
            return None
        x0, x1, x2 = x1, x2, x3
        f0, f1, f2 = f1, f2, f(x3)
        if abs(x2 - x1) < tol * max(1.0, abs(x2)):
            return x2
    return None


def find_zeros_uhp(
    reg: Regulator,
    params: Params,
    re_max: float = _DEFAULT_RE_MAX,
    im_hi: float = _DEFAULT_IM_HI,
    grid: int = 160,
    tol: float = 1e-9,
) -> list[complex]:
    """Localiza los ceros de chi^r con Im(omega) > 0.

    Estrategia: minimos locales de |chi| sobre una malla como semillas, luego
    refinamiento por Muller. Contrastar siempre contra `count_zeros_uhp`, que
    es independiente de la malla.
    """
    re = np.linspace(-re_max, re_max, grid)
    im = np.linspace(1e-4, im_hi, max(grid // 2, 20))
    RE, IM = np.meshgrid(re, im, indexing="ij")
    W = RE + 1j * IM
    A = np.abs(susceptibility(W.ravel(), reg, params)).reshape(W.shape)

    seeds = []
    for i in range(1, A.shape[0] - 1):
        for j in range(1, A.shape[1] - 1):
            if A[i, j] <= A[i - 1:i + 2, j - 1:j + 2].min():
                seeds.append(W[i, j])

    f = lambda w: susceptibility(np.asarray(w), reg, params)
    zeros: list[complex] = []
    for s in seeds:
        h = 1e-3 * max(1.0, abs(s))
        r = _muller(f, s - h, s + h, s + 1j * h, box=(re_max, 0.0, im_hi))
        if r is None or r.imag <= 1e-8 or abs(r.real) > re_max or r.imag > im_hi:
            continue
        if abs(f(r)) > 1e-7:
            continue
        if not any(abs(r - z) < 1e-6 for z in zeros):
            zeros.append(complex(r))

    return sorted(zeros, key=lambda z: (z.imag, z.real))


def is_stable(reg: Regulator, params: Params, **kwargs) -> bool:
    """True si chi^r no tiene ceros en el semiplano superior.

    Este es el criterio linealizado: la linea punteada de la Fig. 2 del paper.
    """
    return count_zeros_uhp(reg, params, **kwargs) == 0


def dominant_pole(reg: Regulator, params: Params, **kwargs):
    """Cero de chi^r con mayor Im(omega), o None si la teoria linealizada es estable.

    Es el modo que gobierna la inestabilidad: Im(omega) es la tasa de
    crecimiento de la autoaceleracion. Al remover el cutoff debe tender a
    3i/2 (en unidades r0=1), el polo runaway de Abraham-Lorentz de la ec. (15).
    """
    zeros = find_zeros_uhp(reg, params, **kwargs)
    return max(zeros, key=lambda z: z.imag) if zeros else None


def critical_cutoff(
    regulator_kind: str,
    lo: float = 0.1,
    hi: float = 40.0,
    tol: float = 1e-3,
    **kwargs,
) -> float:
    """Valor critico de r0/ell donde la teoria linealizada pierde estabilidad.

    Biseccion sobre r0/ell suponiendo monotonia: estable por debajo (cutoff
    suave), inestable por encima. El paper (p. 8) situa la transicion cuando
    el cutoff Lambda = 1/ell es comparable a 1/r0.
    """
    from .core import make_regulator

    def unstable(x: float) -> bool:
        ell = 1.0 / x
        return count_zeros_uhp(make_regulator(regulator_kind, ell),
                               Params(ell=ell), **kwargs) > 0

    if unstable(lo):
        raise ValueError(f"ya inestable en r0/ell={lo}; reduzca lo")
    if not unstable(hi):
        raise ValueError(f"aun estable en r0/ell={hi}; aumente hi")

    while hi - lo > tol:
        mid = 0.5 * (lo + hi)
        lo, hi = (lo, mid) if unstable(mid) else (mid, hi)
    return 0.5 * (lo + hi)



# --------------------------------------------------------------------------
# Rama en el TIEMPO: la misma teoria linealizada, integrada hacia adelante
# --------------------------------------------------------------------------
#
# Linealizando la ec. (6) alrededor de la carga en reposo, x = (s, xi(s)) con
# |xi| pequeno, se tiene a primer orden
#     d = x - x' = (-u, Delta),  u = s'-s <= 0,  Delta = xi(s) - xi(s+u)
#     w = d^2 = u^2 + O(xi^2)
#     V = (0, Delta + u xidot(s+u))            <- solo componente espacial
# de modo que las dos componentes colapsan a UNA ecuacion integro-diferencial
# lineal y escalar:
#
#     xiddot(s) = 4 r_0B int_{-inf}^{0} du delta_B'(u^2)
#                        [xi(s) - xi(s+u) + u xidot(s+u)]                 (*)
#
# Sustituyendo xi = e^{-i omega s} sale la relacion de dispersion, y su
# desarrollo a omega pequeno reproduce los ordenes omega^3 y omega^4 de la
# ec. (14) EXACTAMENTE, pero deja un omega^2 (1 + r_0B M0) en vez de omega^2,
# con M0 = int du delta_B(u^2) = delta_m/(m r0). Es decir (*) es la ecuacion
# DESNUDA: el termino que sobra es el contratermino de masa de la ec. (10).
# Restarlo -- el -phi^2/2 de `_f_disp` -- da D(omega) = omega^2 chi^r_omega,
# lo que se verifica numericamente en tests/test_block1.py.


def _f_disp(phi, renormalizada: bool):
    """1 - (1 + i phi) e^{-i phi}, mas phi^2/2 si se renormaliza la masa."""
    f = 1.0 - (1.0 + 1j * phi) * np.exp(-1j * phi)
    return f + 0.5 * phi ** 2 if renormalizada else f


def _nodos_kernel(reg: Regulator, ell: float, n_ell: float, ds: float):
    """Nodos u <= 0 SOBRE la malla y pesos de Simpson por delta_B'(u^2).

    Elegir u_j = -j ds evita interpolar la historia en `integrate_linear`: la
    memoria queda como una combinacion lineal de los ultimos N valores ya
    calculados. `n` par garantiza que Simpson cierre.
    """
    n = 2 * max(2, int(round(n_ell * ell / (2 * ds))))
    u = -np.arange(n, -1, -1) * ds
    w = np.ones(n + 1); w[1:-1:2] = 4.0; w[2:-1:2] = 2.0
    return u, (ds / 3.0) * w * reg.d_delta(u ** 2)


def dispersion(omega, reg: Regulator, params: Params, renormalizada: bool = True,
               n_ell: float = 60.0, pts_per_ell: int = 64):
    """D(omega) de la ec. (*). Con `renormalizada`, D = omega^2 chi^r_omega.

    Es una via INDEPENDIENTE a los mismos polos que `susceptibility`: esta parte
    de delta_B' en el tiempo, aquella de la ec. (14) ya integrada por partes.
    Que coincidan valida las dos.
    """
    ell = params.ell
    u, k = _nodos_kernel(reg, ell, n_ell, ell / pts_per_ell)
    om = np.atleast_1d(np.asarray(omega, dtype=complex))
    out = om ** 2 + 4.0 * params.r0 * params.m_over_mB * np.einsum(
        "j,ij->i", k, _f_disp(om[:, None] * u[None, :], renormalizada))
    return out[0] if np.ndim(omega) == 0 else out


def integrate_linear(reg: Regulator, params: Params, s_end: float = 40.0,
                     ds: float | None = None, n_ell: float = 30.0,
                     amp: float = 0.05, ancho: float = 1.0, blowup: float = 1e8):
    """Integra la ec. (*) renormalizada hacia adelante. Devuelve (s, xi, xidot, xiddot).

        (m_B/m) xiddot(s) = 4 r0 int du delta_B'(u^2)[xi(s) - xi(s+u) + u xidot(s+u)]
                            + k(s)

    La prehistoria xi = xidot = 0 es solucion EXACTA (carga en reposo), asi que
    no hay salto en s = 0; la excitacion entra por la fuente externa k(s) de la
    ec. (2), un pulso C^infinito de soporte compacto en [0, ancho].

    Solo tiene sentido para el regulador suavizado: el desplazado es una
    distribucion y su delta_B' no admite evaluacion puntual (esa rama es la
    ec. 16, en block2_delay).
    """
    if reg.name != "smeared":
        raise ValueError(
            f"integrate_linear necesita un delta_B' evaluable; el regulador "
            f"{reg.name!r} es distribucional. Use block2_delay para la ec. (4).")

    ell, r0 = params.ell, params.r0
    dm_sobre_m = reg.mass_shift_over_m(r0)
    if abs(1.0 - dm_sobre_m) < 1e-12:
        raise ValueError(
            f"r0/ell = {r0/ell:.4g} da delta_m = m, o sea m_B = 0: el "
            "contratermino diverge y la ecuacion renormalizada no existe.")

    ds = ds or min(ell / 40.0, 5e-3)
    u, k = _nodos_kernel(reg, ell, n_ell, ds)
    n = len(u) - 1
    K0, fac = k.sum(), 4.0 * r0 / (1.0 - dm_sobre_m)      # 1 - dm/m = m_B/m

    xi, xd, xa = (np.zeros(n + 1 + int(np.ceil(s_end / ds))) for _ in range(3))
    s = (np.arange(len(xi)) - n) * ds
    pulso = lambda t: (amp * np.exp(-1.0 / (1.0 - (2 * t / ancho - 1) ** 2))
                       if 0.0 < t < ancho else 0.0)

    for i in range(n, len(xi) - 1):
        xi_p = xi[i] + ds * xd[i] + 0.5 * ds ** 2 * xa[i]
        xd_p, kext = xd[i] + ds * xa[i], pulso(s[i + 1])
        hist = k[:-1] @ (-xi[i + 1 - n:i + 1] + u[:-1] * xd[i + 1 - n:i + 1])
        for _ in range(3):                                # corrector: xi(s) es implicito
            a_new = fac * (K0 * xi_p + hist) + kext
            xd_p = xd[i] + 0.5 * ds * (xa[i] + a_new)
            xi_p = xi[i] + 0.5 * ds * (xd[i] + xd_p)
        xi[i + 1], xd[i + 1], xa[i + 1] = xi_p, xd_p, a_new
        if abs(a_new) > blowup:
            return s[:i + 2], xi[:i + 2], xd[:i + 2], xa[:i + 2]
    return s, xi, xd, xa
