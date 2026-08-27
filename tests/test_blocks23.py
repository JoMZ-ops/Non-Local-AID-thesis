"""Tests de los bloques 2 (ec. 16) y 3 (ec. 17).

El test estructural clave es la ortogonalidad del lado derecho a xdot: es una
identidad exacta de las ecs. (6), (16) y (17) que se verifica contrayendo con
xdot_mu, y que el codigo debe cumplir a precision de maquina. Si una
transcripcion de la ecuacion estuviera mal, este test lo detecta de inmediato.
"""

import math

import numpy as np
import pytest

from nlaid.core import Params, minkowski_dot
from nlaid.worldline import rest_history, smooth_bump, unit_normal
from nlaid.block2_delay import rhs_delay, integrate_delay
from nlaid.block3_memory import rhs_memory, integrate_memory


# --- Utilidades de la linea de mundo --------------------------------------

def test_normal_es_unitario_y_ortogonal():
    """n.v = 0 y n.n = -1: garantiza que la fuente externa preserva xdot^2."""
    for th in [0.0, 0.3, -1.2, 2.0]:
        v = np.array([math.cosh(th), math.sinh(th)])
        n = unit_normal(v)
        assert minkowski_dot(n, v) == pytest.approx(0.0, abs=1e-14)
        assert minkowski_dot(n, n) == pytest.approx(-1.0, rel=1e-14)


def test_pulso_tiene_soporte_compacto():
    assert smooth_bump(-0.1) == 0.0
    assert smooth_bump(1.1) == 0.0
    assert smooth_bump(0.5) > 0.0


def test_prehistoria_en_reposo_normalizada():
    wl = rest_history(dim=2, s_rest=5.0, ds=1e-2)
    assert wl.norm_drift.max() < 1e-15


# --- Propiedad estructural: ortogonalidad ---------------------------------

def _estado_perturbado(seed=0, ell=0.3):
    """Estado generico fuera de reposo, para no testear un caso degenerado."""
    rng = np.random.default_rng(seed)
    wl = rest_history(dim=2, s_rest=20.0, ds=2e-3)
    p = Params(ell=ell, m_over_mB=0.5)
    wl = integrate_delay(p, s_end=2.0, ds=2e-3, history=wl)
    return wl, p


def test_rhs_ec16_ortogonal_a_xdot():
    """xddot . xdot = 0 exactamente para la ec. (16)."""
    wl, p = _estado_perturbado()
    x, v = wl._x[-1], wl._v[-1]
    _, xp, vp, ap = wl.retarded_point(wl.s_max + 1e-9, x, p.ell)
    acc, _ = rhs_delay(x, v, xp, vp, ap, p)
    assert abs(minkowski_dot(acc, v)) < 1e-12 * max(1.0, np.linalg.norm(acc))


def test_rhs_ec17_ortogonal_a_xdot():
    """xddot . xdot = 0 exactamente para la ec. (17)."""
    wl = rest_history(dim=2, s_rest=30.0, ds=2e-3)
    p = Params(ell=0.3, m_over_mB=0.5)
    wl = integrate_memory(p, s_end=2.0, ds=2e-3, history=wl, n_ell=20.0)
    acc = rhs_memory(wl._x[-1], wl._v[-1], wl, wl.s_max, p, n_ell=20.0)
    assert abs(minkowski_dot(acc, wl._v[-1])) < 1e-12 * max(1.0, np.linalg.norm(acc))


# --- El reposo es solucion exacta -----------------------------------------

def test_reposo_no_tiene_autofuerza_ec16():
    """Movimiento inercial => V1 = 0 => sin autofuerza. Sin esto la
    prehistoria introduciria un transitorio espurio."""
    wl = rest_history(dim=2, s_rest=10.0, ds=1e-3)
    p = Params(ell=0.4)
    x, v = wl._x[-1], wl._v[-1]
    _, xp, vp, ap = wl.retarded_point(wl.s_max + 1e-9, x, p.ell)
    acc, _ = rhs_delay(x, v, xp, vp, ap, p)
    assert np.max(np.abs(acc)) < 1e-12


def test_reposo_no_tiene_autofuerza_ec17():
    wl = rest_history(dim=2, s_rest=30.0, ds=1e-3)
    p = Params(ell=0.4)
    acc = rhs_memory(wl._x[-1], wl._v[-1], wl, wl.s_max, p, n_ell=20.0)
    assert np.max(np.abs(acc)) < 1e-12


def test_punto_retardado_en_reposo_es_exactamente_ell():
    """En reposo (x - x')^2 = (s - s')^2, luego el retardo vale ell."""
    wl = rest_history(dim=2, s_rest=10.0, ds=1e-3)
    ell = 0.37
    sp, *_ = wl.retarded_point(wl.s_max, wl._x[-1], ell)
    assert wl.s_max - sp == pytest.approx(ell, rel=1e-9)


# --- Convergencia ----------------------------------------------------------

def _orden(sols, ds_list, probe):
    ref = sols[ds_list[-1]].sample_many(probe)[0]
    e = [np.max(np.abs(sols[d].sample_many(probe)[0] - ref)) for d in ds_list[:-1]]
    return [math.log2(e[i] / e[i + 1]) for i in range(len(e) - 1)]


@pytest.mark.slow
def test_ec16_converge_a_segundo_orden():
    ds_list = [2e-2, 1e-2, 5e-3, 2.5e-3]
    sols = {d: integrate_delay(Params(ell=0.5, m_over_mB=0.4), s_end=4.0, ds=d)
            for d in ds_list}
    ordenes = _orden(sols, ds_list, np.linspace(1.5, 4.0, 40))
    assert min(ordenes) > 1.7, f"ordenes={ordenes}"


@pytest.mark.slow
def test_ec17_converge_a_segundo_orden():
    ds_list = [2e-2, 1e-2, 5e-3]
    sols = {d: integrate_memory(Params(ell=0.5, m_over_mB=0.4), s_end=3.0, ds=d,
                                n_ell=25.0) for d in ds_list}
    ordenes = _orden(sols, ds_list, np.linspace(1.5, 3.0, 30))
    assert min(ordenes) > 1.6, f"ordenes={ordenes}"


def test_ec17_ventana_de_memoria_converge():
    """El truncamiento del nucleo es exponencial en n_ell.

    Con pts_per_ell fijo la resolucion de la cuadratura no cambia al ampliar
    la ventana, de modo que la diferencia mide truncamiento puro. Debe caer,
    no crecer, al ampliar la ventana.
    """
    v = {n: np.linalg.norm(
        integrate_memory(Params(ell=0.5, m_over_mB=0.4), s_end=2.0, ds=1e-2,
                         n_ell=n, pts_per_ell=16)._a[-1]) for n in (10, 20, 30)}
    d1, d2 = abs(v[10] - v[20]), abs(v[20] - v[30])
    assert d2 < d1, f"el truncamiento no mejora: {d1:.2e} -> {d2:.2e}"
    assert d2 < 1e-6 * v[30]      # medido: ~4e-7 relativo entre n_ell 20 y 30


# --- Comportamiento fisico -------------------------------------------------

def test_acoplamiento_debil_es_estable():
    """El paper: dinamica estable para |m/m_B| pequeno. La aceleracion decae."""
    wl = integrate_delay(Params(ell=0.5, m_over_mB=0.2), s_end=8.0, ds=5e-3)
    m = wl.s >= 1.2
    a = np.linalg.norm(wl.a[m], axis=1)
    assert a[-1] < 0.05 * a.max()
    assert wl.norm_drift.max() < 1e-6


def test_acoplamiento_fuerte_no_relaja():
    """Contraste cualitativo con el test anterior.

    Fuera de la region de estabilidad la aceleracion NO relaja: crece un orden
    de magnitud y se mantiene, en lugar de decaer a cero. (No se afirma
    divergencia: en la ec. (16) los terminos no lineales acotan el crecimiento
    -- esa es justamente la tesis del paper -- de modo que el criterio correcto
    es ausencia de relajacion, no escape a infinito.)
    """
    wl = integrate_delay(Params(ell=0.5, m_over_mB=1.0), s_end=12.0, ds=5e-3)
    m = wl.s >= 1.2
    a = np.linalg.norm(wl.a[m], axis=1)
    assert a[-1] > 5 * a[0], "deberia no relajar"
    assert a.max() > 10 * a[0]
