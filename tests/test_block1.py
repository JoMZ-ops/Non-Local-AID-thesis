"""Tests del nucleo y del Bloque 1 contra resultados analiticos.

Cada test contrasta el codigo con un valor cerrado derivado a mano del paper,
no con otra corrida del propio codigo.
"""

import numpy as np
import pytest

from nlaid.core import Params, ShiftedDelta, SmearedDelta
from nlaid.block1_linear import (
    numerator_N,
    susceptibility,
    susceptibility_small_omega,
    count_zeros_uhp,
    find_zeros_uhp,
)


# --- Regulador suavizado, ec. (5): momentos analiticos --------------------

@pytest.mark.parametrize("ell", [0.05, 0.2, 1.0, 3.0])
def test_smeared_normalizacion(ell):
    """Condicion (i): int dz delta_B(z) = 1."""
    reg = SmearedDelta(ell)
    assert reg.check_conditions()["cumple"]


@pytest.mark.parametrize("ell", [0.05, 0.2, 1.0, 3.0])
def test_smeared_momento_inv_sqrt(ell):
    """int_0^inf dz z^{-1/2} delta_B(z) = 1/(3 ell)  =>  delta_m/m = r0/(6 ell)."""
    from scipy.integrate import quad
    reg = SmearedDelta(ell)
    num = quad(lambda z: reg.delta(np.float64(z)) / np.sqrt(z), 0, np.inf, limit=400)[0]
    assert num == pytest.approx(1.0 / (3.0 * ell), rel=1e-8)
    assert reg.mass_shift_over_m(r0=1.0) == pytest.approx(1.0 / (6.0 * ell), rel=1e-12)


@pytest.mark.parametrize("ell", [0.05, 0.2, 1.0, 3.0])
def test_smeared_momento_u2(ell):
    """I2 = int_{-inf}^0 du delta_B(u^2) u^2 = 2 ell (analitico)."""
    assert SmearedDelta(ell).moment_u2() == pytest.approx(2.0 * ell, rel=1e-10)


@pytest.mark.parametrize("ell", [0.05, 0.2, 1.0, 3.0])
def test_shifted_momento_u2(ell):
    """Para la delta desplazada, I2 = (-ell)^2 / (2 ell) = ell/2."""
    assert ShiftedDelta(ell).moment_u2() == pytest.approx(0.5 * ell, rel=1e-14)


# --- Numerador N: serie vs forma directa ----------------------------------

def test_N_serie_coincide_con_forma_directa():
    """Las dos ramas de numerator_N deben empalmar en la frontera."""
    phi = np.array([5.9, 6.1, -5.9, -6.1, 5.9 + 0.5j, 6.1 + 0.5j])
    directa = ((1 + 1j * phi - phi**2) * np.exp(-1j * phi)
               - 1 + 0.5 * phi**2 - (2/3) * 1j * phi**3)
    assert np.allclose(numerator_N(phi), directa, rtol=1e-11, atol=1e-13)


def test_N_orden_dominante():
    """N(phi) = (3/8) phi^4 - (2/15) i phi^5 + O(phi^6)."""
    phi = 1e-3
    esperado = (3/8) * phi**4 - (2/15) * 1j * phi**5
    assert numerator_N(phi) == pytest.approx(esperado, rel=1e-6)


def test_N_sin_cancelacion_catastrofica():
    """A phi muy pequeno la forma directa pierde toda la precision; la serie no."""
    phi = 1e-4
    directa = ((1 + 1j*phi - phi**2) * np.exp(-1j*phi)
               - 1 + 0.5*phi**2 - (2/3)*1j*phi**3)
    exacto = (3/8) * phi**4
    assert abs(numerator_N(phi).real - exacto) / exacto < 1e-8
    assert abs(directa.real - exacto) / exacto > 1e-4      # la directa ya fallo


# --- Susceptibilidad -------------------------------------------------------

@pytest.mark.parametrize("kind,ell", [("smeared", 0.3), ("smeared", 1.0),
                                      ("shifted", 0.3), ("shifted", 1.0)])
def test_chi_desarrollo_omega_pequeno_es_orden_tres(kind, ell):
    """chi - chi_small debe anularse como O(w^3), no solo ser pequeno.

    Verificar el ORDEN es mas fuerte que verificar la cercania: confirma que
    el coeficiente -3/4 I2 del termino w^2 es exacto. Si estuviera mal, el
    residuo seria O(w^2) y el cociente daria 4, no 8.
    """
    reg = SmearedDelta(ell) if kind == "smeared" else ShiftedDelta(ell)
    p = Params(ell=ell)
    w1, w2 = 2e-3, 1e-3
    r1 = abs(susceptibility(w1, reg, p) - susceptibility_small_omega(w1, reg, p))
    r2 = abs(susceptibility(w2, reg, p) - susceptibility_small_omega(w2, reg, p))
    assert r1 / r2 == pytest.approx(8.0, rel=0.05)


def test_chi_recupera_abraham_lorentz():
    """Ec. (15): al remover el cutoff, chi -> 1 + (2/3) i r0 omega."""
    p = Params(ell=1e-4)
    reg = SmearedDelta(1e-4)
    w = 0.05
    assert susceptibility(w, reg, p) == pytest.approx(1 + (2/3)*1j*w, rel=1e-3)


def test_chi_shifted_forma_cerrada():
    """Para la delta desplazada la integral es exacta: chi = 1 + r0[(2/3)i w - N(-w ell)/(w^2 ell^3)]."""
    ell = 0.4
    reg, p = ShiftedDelta(ell), Params(ell=ell)
    w = np.array([0.2, 1.0, 3.0, 8.0 + 0.3j])
    cerrada = 1 + p.r0 * ((2/3)*1j*w - numerator_N(-w*ell) / (w**2 * ell**3))
    assert np.allclose(susceptibility(w, reg, p), cerrada, rtol=1e-12)


# --- Estabilidad -----------------------------------------------------------

def test_polo_runaway_de_abraham_lorentz():
    """Con cutoff muy pequeno debe aparecer el cero de chi en omega ~ 1.5 i / r0.

    De la ec. (15): 1 + (2/3) i r0 omega = 0  =>  omega = 3i/(2 r0).
    Es el polo acausal responsable de las trayectorias autoaceleradas.
    """
    ell = 0.02
    reg, p = SmearedDelta(ell), Params(ell=ell)
    ceros = find_zeros_uhp(reg, p, re_max=6.0, im_hi=6.0, grid=120)
    assert ceros, "se esperaba al menos un cero en el semiplano superior"
    assert min(abs(z - 1.5j) for z in ceros) < 0.15


def test_argumento_y_muller_coinciden():
    """El conteo global (principio del argumento) debe igualar los ceros hallados."""
    ell = 0.02
    reg, p = SmearedDelta(ell), Params(ell=ell)
    n = count_zeros_uhp(reg, p, re_max=6.0, im_hi=6.0, n=600)
    ceros = find_zeros_uhp(reg, p, re_max=6.0, im_hi=6.0, grid=120)
    assert n == len(ceros)


def test_cutoff_grande_es_estable():
    """Para ell >> r0 la teoria linealizada debe ser estable (sin ceros en UHP)."""
    ell = 20.0
    reg, p = SmearedDelta(ell), Params(ell=ell)
    assert count_zeros_uhp(reg, p, re_max=5.0, im_hi=0.9/ell, n=600) == 0


# --- Funciones que consumen las figuras ------------------------------------
# Estaban sin cobertura pese a ser las que producen los numeros graficados.

@pytest.mark.parametrize("kind,ell,esperado", [
    ("shifted", 1.0, 2.0),          # delta_m/m = r0/2ell = 0.5  -> 1/(1-0.5)
    ("shifted", 0.25, -1.0),        # r0/ell = 4  -> delta_m/m = 2
    ("smeared", 1.0, 1.2),          # delta_m/m = r0/6ell = 1/6  -> 1/(1-1/6)
    ("smeared", 1.0 / 12, -1.0),    # r0/ell = 12 -> delta_m/m = 2
])
def test_contraterm_forma_cerrada(kind, ell, esperado):
    """m/m_B = 1/(1 - delta_m/m) con delta_m/m = r0/2ell o r0/6ell.

    Digitalizando la Fig. 2(a) del paper, la ley del regulador desplazado
    reproduce su linea punteada con error medio 0.0017 sobre siete puntos
    (ver docs/discrepancias.md, hallazgo F1).
    """
    from nlaid.core import make_regulator
    assert make_regulator(kind, ell).m_over_mB_counterterm() == pytest.approx(
        esperado, rel=1e-12)


def test_contraterm_diverge_donde_el_bare_mass_se_anula():
    """delta_m = m  =>  m_B = 0  =>  m/m_B diverge.

    Ocurre en r0/ell = 2 (desplazado) y 6 (suavizado): son las asintotas
    verticales que separan las dos ramas de la Fig. 2.
    """
    from nlaid.core import make_regulator
    for kind, x_c in (("shifted", 2.0), ("smeared", 6.0)):
        assert make_regulator(kind, 1.0 / x_c).mass_shift_over_m() == pytest.approx(1.0)
        # justo a un lado y al otro el signo se invierte y la magnitud es grande
        izq = make_regulator(kind, 1.0 / (x_c * 0.999)).m_over_mB_counterterm()
        der = make_regulator(kind, 1.0 / (x_c * 1.001)).m_over_mB_counterterm()
        assert izq > 100 and der < -100, (kind, izq, der)


def test_dominant_pole_devuelve_el_cero_de_mayor_parte_imaginaria():
    from nlaid.block1_linear import dominant_pole, find_zeros_uhp
    ell = 0.02
    reg, p = SmearedDelta(ell), Params(ell=ell)
    kw = dict(re_max=6.0, im_hi=6.0, grid=120)
    z = dominant_pole(reg, p, **kw)
    assert z is not None
    assert z.imag == max(w.imag for w in find_zeros_uhp(reg, p, **kw))


def test_dominant_pole_es_None_si_la_teoria_es_estable():
    """Sin ceros en el semiplano superior no hay modo inestable que devolver."""
    from nlaid.block1_linear import dominant_pole
    ell = 20.0
    assert dominant_pole(SmearedDelta(ell), Params(ell=ell),
                         re_max=5.0, im_hi=5.0, grid=80) is None


@pytest.mark.slow
def test_critical_cutoff_separa_estable_de_inestable():
    """El valor critico debe tener conteo 0 justo debajo y >0 justo encima.

    Es el contraste que valida la biseccion sin suponer su monotonia: se
    comprueba directamente el cambio de fase en el punto devuelto.
    """
    from nlaid.core import make_regulator
    from nlaid.block1_linear import critical_cutoff, count_zeros_uhp
    kw = dict(re_max=8.0, im_hi=8.0)
    xc = critical_cutoff("shifted", lo=0.3, hi=20.0, tol=5e-3, n=800, **kw)
    for factor, esperado_inestable in ((0.95, False), (1.05, True)):
        ell = 1.0 / (xc * factor)
        n = count_zeros_uhp(make_regulator("shifted", ell), Params(ell=ell),
                            n=800, **kw)
        assert (n > 0) is esperado_inestable, (factor, n, xc)


# --- Rama en el tiempo: la ec. linealizada frente a la ec. (14) -----------

@pytest.mark.parametrize("x", [2.0, 4.0, 6.0])
@pytest.mark.parametrize("omega", [0.3 + 0.1j, 1.5 - 0.4j, 3.0 + 1.0j])
def test_dispersion_es_omega2_por_chi(x, omega):
    """D(omega) de la ec. linealizada en el tiempo == omega^2 chi^r_omega.

    Son DOS derivaciones distintas del mismo objeto: `dispersion` parte de
    delta_B'(u^2) en el dominio temporal, `susceptibility` de la ec. (14) del
    paper, ya integrada por partes. Que coincidan valida ambas ramas, y en
    particular que la resta de -phi^2/2 sea exactamente el contratermino de
    masa de la ec. (10).
    """
    from nlaid.block1_linear import dispersion
    ell = 1.0 / x
    reg, pr = SmearedDelta(ell), Params(ell=ell)
    a = dispersion(omega, reg, pr)
    b = omega ** 2 * susceptibility(omega, reg, pr)
    assert abs(a - b) / abs(b) < 1e-7, (x, omega, a, b)


def test_dispersion_sin_renormalizar_deja_el_contratermino():
    """Sin renormalizar, el coeficiente de omega^2 es 1 + delta_m/m, no 1.

    Es el enunciado preciso de que la ec. linealizada cruda es la DESNUDA.
    Se lee a omega pequeno, donde los ordenes omega^3 y omega^4 son despreciables.
    """
    from nlaid.block1_linear import dispersion
    ell = 1.0 / 3.0
    reg, pr = SmearedDelta(ell), Params(ell=ell)
    om = 1e-3
    razon = dispersion(om, reg, pr, renormalizada=False) / om ** 2
    assert abs(razon.real - (1.0 + reg.mass_shift_over_m())) < 1e-6, razon


def test_integrate_linear_relaja_por_debajo_del_critico():
    """Debajo del cutoff critico (r0/ell = 4) la aceleracion decae; encima crece."""
    from nlaid.core import make_regulator
    from nlaid.block1_linear import integrate_linear
    for x, decae in ((2.0, True), (5.0, False)):
        ell = 1.0 / x
        s, xi, xd, xa = integrate_linear(make_regulator("smeared", ell),
                                         Params(ell=ell), s_end=18.0)
        m = s > 4.0
        tasa = np.polyfit(s[m], np.log(np.abs(xa[m]) + 1e-300), 1)[0]
        assert bool(tasa < 0) is decae, (x, tasa)


def test_integrate_linear_rechaza_el_regulador_distribucional():
    """delta_B' de la ec. (4) no admite evaluacion puntual: debe fallar claro."""
    from nlaid.block1_linear import integrate_linear
    with pytest.raises(ValueError, match="distribucional"):
        integrate_linear(ShiftedDelta(0.5), Params(ell=0.5), s_end=1.0)
