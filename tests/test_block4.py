"""Tests del bloque 4: medida de la tasa de relajacion y shooting sobre m/m_B.

El test central es de consistencia entre bloques: la tasa medida sobre la
dinamica NO LINEAL (ec. 17) debe coincidir con Im(omega) del cero dominante
de chi^r (ec. 14), que se calcula por analisis espectral linealizado. Son dos
caminos independientes hacia la misma cantidad.
"""

import numpy as np
import pytest

from nlaid.core import Params
from nlaid.block4_renorm import relaxation_rate, linearized_rate, solve_ratio

ELL3 = 1.0 / 3.0          # r0/ell = 3, el caso de la Fig. 1 del paper


def test_tasa_coincide_con_la_teoria_linealizada():
    """Contraste bloque 3 vs bloque 1 a r0/ell = 3, m/m_B = 2.

    Referencia establecida: el cero dominante de chi^r esta en
    omega = +-2.5736 - 0.71793i.
    """
    tasa, diag = relaxation_rate(Params(ell=ELL3, m_over_mB=2.0),
                                 s_end=14.0, ds=5e-3)
    assert diag["fiable"], diag
    assert tasa == pytest.approx(-0.71793, abs=5e-3)


def test_linearized_rate_reproduce_el_cero_dominante():
    assert linearized_rate(Params(ell=ELL3)) == pytest.approx(-0.71793, abs=1e-4)


def test_la_tasa_crece_monotonamente_con_el_acoplamiento():
    """Mas acoplamiento -> relajacion mas lenta. Es lo que permite el shooting."""
    tasas = [relaxation_rate(Params(ell=ELL3, m_over_mB=mr),
                             s_end=10.0, ds=1e-2)[0]
             for mr in (1.0, 1.5, 2.0, 2.5)]
    assert all(np.isfinite(tasas)), tasas
    assert all(b > a for a, b in zip(tasas, tasas[1:])), tasas


def test_una_trayectoria_que_escapa_se_trunca_y_se_clasifica_inestable():
    """m/m_B = -3.91 a r0/ell = 3 escapa: |xddot| llega a 1e6 hacia s ~ 14.

    El integrador corta al superar BLOWUP en vez de seguir acumulando error,
    de modo que la deriva de xdot^2 se mantiene pequena -- el corte ocurre
    ANTES de que la normalizacion se degrade. La senal de inestabilidad es por
    tanto la historia truncada y la amplitud enorme, no la deriva.
    """
    from nlaid.block3_memory import BLOWUP, integrate_memory
    from nlaid.block4_renorm import _es_inestable

    wl = integrate_memory(Params(ell=ELL3, m_over_mB=-3.91),
                          s_end=30.0, ds=1e-2)
    assert wl.s_max < 30.0, "deberia haberse truncado por escape"
    assert np.linalg.norm(wl.a[-1]) > BLOWUP

    inestable, _ = _es_inestable(-3.91, ELL3, "smeared", s_end=30.0, ds=1e-2)
    assert inestable is True


@pytest.mark.slow
def test_shooting_recupera_un_objetivo_conocido():
    """Pedirle al solver la tasa que ya sabemos que da m/m_B = 2 debe devolver 2.

    Es un test de ida y vuelta: fija el objetivo en el valor medido y comprueba
    que el shooting reconstruye el parametro de partida.
    """
    mr, diag = solve_ratio(ELL3, -0.71793, rama="positiva",
                           regulator="smeared", s_end=14.0, ds=5e-3)
    assert np.isfinite(mr), diag
    assert mr == pytest.approx(2.0, abs=0.05)
