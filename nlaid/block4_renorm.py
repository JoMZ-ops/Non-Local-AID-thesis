"""
BLOQUE 4 -- Condicion de renormalizacion y estructura de fases.

La ec. (18) del paper, chi^r_omega = 1 + (2/3) i r0 omega, NO puede imponerse
literalmente sobre una teoria estable: su unico cero esta en omega = 3i/2, en
el semiplano SUPERIOR, es decir el polo runaway de Abraham-Lorentz. Una teoria
que la satisfaga exactamente no tiene modo de relajacion que monitorear.

El paper la aplica operacionalmente (p. 10):
    "the renormalization condition, (18), can be satisfied by monitoring the
     relaxation for large s"                                    [m_B > 0]
    "The relaxation of the envelope is used to find the physical theory"
                                                                [m_B < 0]

Se implementan por tanto DOS condiciones bien definidas, ambas por shooting
sobre m/m_B a cutoff fijo:

  R1  "empalme con la teoria linealizada"
      tasa no lineal  =  Im(omega) del cero dominante de chi^r (ec. 14).
      Es la lectura mas fiel de "monitorear la relajacion": exige que la
      dinamica no lineal reproduzca la relajacion que predice la teoria
      linealizada ya renormalizada.

  R2  "borde de estabilidad"
      tasa no lineal  =  0.
      Es la frontera de la region sombreada de la Fig. 2, y el objeto mas
      directamente comparable con el paper: las tripletas de la Fig. 1
      bracketean justamente un cambio de signo de la tasa.

Ambas devuelven dos ramas, una con m_B > 0 y otra con m_B < 0, como describe
el paper.
"""

from __future__ import annotations

import numpy as np

from .core import Params, make_regulator
from .block1_linear import dominant_pole, find_zeros_uhp, susceptibility
from .block2_delay import integrate_delay
from .block3_memory import integrate_memory

__all__ = [
    "relaxation_rate",
    "linearized_rate",
    "solve_ratio",
    "stability_edge",
    "renormalized_ratio",
]

MAX_DRIFT = 1e-3          # deriva de xdot^2 sobre la que la medida no es fiable
AMP_FLOOR = 1e-200        # por debajo de esto la senal es ruido de redondeo
R2_MIN = 0.70             # calidad minima del ajuste exponencial
AMP_SANA = 1e-6           # por encima de esto la senal es fuerte y una
                          # envolvente plana es fisica (borde de estabilidad),
                          # no ruido de redondeo


class _NoFiable(Exception):
    """La integracion diverge o pierde la normalizacion en este m/m_B."""

    def __init__(self, mr, diag):
        super().__init__(f"medida no fiable en m/m_B={mr}")
        self.mr, self.diag = mr, diag


def relaxation_rate(
    params: Params,
    regulator: str = "smeared",
    s_end: float = 14.0,
    ds: float = 5e-3,
    fit_from: float = 0.60,
    **kw,
):
    """Tasa exponencial de la envolvente de |xddot| a s grande.

    Devuelve (tasa, diagnosticos). La tasa se ajusta sobre los MAXIMOS locales
    de |xddot|, no sobre la senal: la solucion oscila, y ajustar la senal cruda
    mezcla la envolvente con la fase.

    `diag["drift"]` reporta la deriva de xdot^2; por encima de MAX_DRIFT la
    medida se marca como no fiable en vez de devolverse en silencio.
    """
    integra = integrate_memory if regulator == "smeared" else integrate_delay
    wl = integra(params, s_end=s_end, ds=ds, **kw)

    drift = float(wl.norm_drift.max())
    # La ventana de ajuste se define sobre el s REALMENTE alcanzado, no sobre
    # el pedido: los integradores cortan antes si la trayectoria escapa, y con
    # `s_end` la ventana quedaria vacia en esos casos.
    m = wl.s >= fit_from * wl.s_max
    s, a = wl.s[m], np.linalg.norm(wl.a[m], axis=1)

    # maximos locales -> envolvente
    if len(a) > 2:
        pk = np.nonzero((a[1:-1] >= a[:-2]) & (a[1:-1] >= a[2:]))[0] + 1
    else:
        pk = np.array([], dtype=int)

    if len(pk) >= 4:
        ss, aa = s[pk], a[pk]
    else:
        ss, aa = s, a                       # sin oscilacion resoluble

    good = aa > AMP_FLOOR
    if good.sum() < 3:
        return np.nan, {"drift": drift, "fiable": False,
                        "motivo": "senal por debajo del piso numerico",
                        "a_max": float(a.max()), "r2": np.nan}

    ls, la = ss[good], np.log(aa[good])
    coef = np.polyfit(ls, la, 1)
    tasa = float(coef[0])

    # Calidad del ajuste. Sin este control, una senal que ya decayo hasta el
    # piso de redondeo produce una recta casi plana que el shooting confunde
    # con tasa = 0, es decir con el borde de estabilidad: a acoplamiento debil
    # aparecian raices espurias en el propio punto inicial del bracket.
    resid = la - np.polyval(coef, ls)
    ss_tot = float(np.sum((la - la.mean()) ** 2))
    r2 = 1.0 - float(np.sum(resid ** 2)) / ss_tot if ss_tot > 0 else 0.0

    # Criterio de fiabilidad. Hay que separar dos situaciones que producen la
    # MISMA firma superficial -- una envolvente casi plana, r2 bajo:
    #
    #   (i)  borde de estabilidad: la envolvente es plana por definicion. r2
    #        tiende a cero justo en el punto que el shooting busca, asi que
    #        exigir buen r2 rechazaria la raiz. La senal es sana (amplitud
    #        O(1e-2)) y la medida es VALIDA.
    #   (ii) senal muerta: ya decayo hasta el ruido de redondeo y el ajuste da
    #        una recta espuria casi plana. La amplitud es minuscula y la
    #        medida es INVALIDA -- este era el origen de las raices espurias
    #        que aparecian en el propio punto inicial del bracket.
    #
    # Se acepta por tanto cuando el ajuste es bueno (decaimiento exponencial
    # limpio, cualquier amplitud por encima del piso) O cuando la senal es
    # fuerte (caso i, donde el mal r2 es fisico y no numerico).
    amax = float(a.max())
    fiable = (drift < MAX_DRIFT and np.isfinite(tasa)
              and ((r2 > R2_MIN and amax > AMP_FLOOR * 1e3) or amax > AMP_SANA))
    return tasa, {
        "drift": drift,
        "fiable": fiable,
        "r2": r2,
        "n_picos": int(len(pk)),
        "a_max": float(a.max()),
        "a_final": float(a[-1]),
    }


def linearized_rate(params: Params, regulator: str = "smeared") -> float:
    """Im(omega) del cero dominante de chi^r (ec. 14).

    Si hay ceros en el semiplano superior la teoria linealizada es inestable y
    se devuelve esa tasa (positiva). Si no, se busca el cero menos amortiguado
    del semiplano inferior, que es el que gobierna la relajacion.
    """
    reg = make_regulator(regulator, params.ell)

    arriba = find_zeros_uhp(reg, params, re_max=12.0, im_hi=12.0, grid=180)
    if arriba:
        return max(z.imag for z in arriba)

    # semiplano inferior: malla + refinamiento
    from .block1_linear import _muller
    f = lambda w: susceptibility(np.asarray(w), reg, params)
    lo = max(-0.97 / params.ell, -8.0)
    re = np.linspace(-8, 8, 130)
    im = np.linspace(lo, -0.02, 70)
    RE, IM = np.meshgrid(re, im, indexing="ij")
    W = RE + 1j * IM
    A = np.abs(f(W.ravel())).reshape(W.shape)

    # Solo interesa el cero MENOS amortiguado, asi que se refinan las semillas
    # empezando por las de mayor Im y se corta en cuanto una converge: refinar
    # cada minimo local de la malla cuesta decenas de miles de evaluaciones de
    # chi y no aporta nada al resultado.
    semillas = [W[i, j]
                for i in range(1, A.shape[0] - 1)
                for j in range(1, A.shape[1] - 1)
                if A[i, j] <= A[i - 1:i + 2, j - 1:j + 2].min()]
    semillas.sort(key=lambda w: -w.imag)

    for s0 in semillas[:24]:
        r = _muller(f, s0 - 1e-3, s0 + 1e-3, s0 + 1e-3j, box=(8.0, lo, 0.0))
        if r is None or r.imag >= 0 or abs(f(r)) > 1e-7:
            continue
        return float(r.imag)
    return np.nan


def solve_ratio(
    ell: float,
    objetivo: float,
    rama: str = "positiva",
    regulator: str = "smeared",
    lo: float | None = None,
    hi: float | None = None,
    tol: float = 3e-3,
    verbose: bool = False,
    **kw,
):
    """Resuelve m/m_B para que `relaxation_rate` iguale `objetivo`.

    Usa Brent (`scipy.optimize.brentq`) en lugar de biseccion pura: la tasa es
    suave y monotona en m/m_B, de modo que Brent converge en ~6 evaluaciones
    donde la biseccion necesita ~20, y cada evaluacion cuesta una integracion
    completa. Las evaluaciones se cachean por si el bracketing las repite.

    `rama`: 'positiva' -> m_B > 0,  'negativa' -> m_B < 0.
    Devuelve (m_over_mB, diagnosticos); nan si no se logra acotar la raiz.
    """
    from scipy.optimize import brentq

    cache: dict[float, tuple[float, dict]] = {}

    def evalua(mr: float):
        k = round(float(mr), 9)
        if k not in cache:
            t, d = relaxation_rate(Params(ell=ell, m_over_mB=k), regulator, **kw)
            cache[k] = (t, d)
            if verbose:
                print(f"      m/mB={k:+8.4f}  tasa={t:+9.5f}  "
                      f"deriva={d['drift']:.1e}  fiable={d['fiable']}", flush=True)
        return cache[k]

    def g(mr: float) -> float:
        t, d = evalua(mr)
        if not d["fiable"] or not np.isfinite(t):
            raise _NoFiable(mr, d)
        return t - objetivo

    # Bracketing por expansion MULTIPLICATIVA desde acoplamiento debil.
    # Debe arrancar muy cerca de cero: el borde de estabilidad de la ec. (16)
    # cae en m/m_B ~ 0.2 a cutoff alto, de modo que un arranque en 0.3 ya
    # estaria del lado inestable y la raiz no se acotaria. La expansion
    # multiplicativa cubre el rango 0.02..24 en ~16 pasos.
    signo = 1.0 if rama == "positiva" else -1.0
    a = 0.02 * signo if lo is None else lo
    tope = 24.0 * signo if hi is None else hi

    # El extremo debil del bracket no puede fijarse a ciegas: a acoplamiento
    # muy pequeno el amortiguamiento es tan fuerte que la senal cae al ruido de
    # redondeo y la medida se rechaza. Se sube hasta el primer punto FIABLE.
    # Ese punto esta siempre del lado estable (tasa < 0): la region no fiable
    # es la de amortiguamiento fuerte, es decir la mas alejada del borde.
    ga = None
    for _ in range(16):
        try:
            ga = g(a)
            break
        except _NoFiable:
            a *= 1.55
    if ga is None:
        return np.nan, {"motivo": "sin punto fiable en el extremo debil",
                        "evaluaciones": len(cache)}
    if ga > 0:
        return np.nan, {"motivo": "el primer punto fiable ya es inestable",
                        "a": a, "g_a": ga, "evaluaciones": len(cache)}

    b, gb = a, ga
    for _ in range(18):
        nb = b * 1.55
        if abs(nb) > abs(tope):
            break
        try:
            gnb = g(nb)
        except _NoFiable:
            # mas alla la dinamica escapa: la raiz queda entre b y nb
            break
        b, gb = nb, gnb
        if ga * gb < 0:
            break
    else:
        pass

    if not (np.isfinite(gb) and ga * gb < 0):
        return np.nan, {"motivo": "raiz no acotada",
                        "g_a": ga, "g_b": gb, "a": a, "b": b,
                        "evaluaciones": len(cache)}

    try:
        mr = brentq(g, a, b, xtol=tol, rtol=1e-8, maxiter=60)
    except (_NoFiable, ValueError) as e:
        return np.nan, {"motivo": f"brentq fallo: {type(e).__name__}",
                        "evaluaciones": len(cache)}

    t, d = evalua(mr)
    return float(mr), {"tasa": t, "objetivo": objetivo,
                       "evaluaciones": len(cache), **d}


def _es_inestable(mr: float, ell: float, regulator: str, **kw):
    """Clasifica un m/m_B como inestable, estable, o indeterminado.

    Devuelve (bool | None, diagnosticos). La clasificacion NO depende de
    ajustar bien una exponencial: cerca del borde la envolvente es plana y el
    ajuste pierde sentido justo donde mas importa. Se usan en cambio tres
    senales robustas:

      - deriva de xdot^2 por encima del umbral  ->  la trayectoria escapo,
        INESTABLE (la integracion se rompe porque la solucion diverge);
      - medida fiable  ->  el signo de la tasa decide;
      - senal muerta (decaida al ruido de redondeo)  ->  fuertemente
        amortiguada, ESTABLE.
    """
    tasa, d = relaxation_rate(Params(ell=ell, m_over_mB=mr), regulator, **kw)
    if d["drift"] >= MAX_DRIFT:
        return True, d
    if d["fiable"] and np.isfinite(tasa):
        return tasa > 0.0, d
    if d.get("a_max", 0.0) <= AMP_SANA:
        return False, d
    return None, d


def stability_edge(ell: float, rama: str = "positiva",
                   regulator: str = "smeared", tol: float = 3e-3,
                   hi: float = 30.0, guess: float | None = None,
                   verbose: bool = False, **kw):
    """R2: m/m_B donde la dinamica pasa de estable a inestable.

    Es la frontera de la region sombreada de la Fig. 2. Se biseca la
    CLASIFICACION estable/inestable en lugar de resolver tasa = 0: cerca del
    borde la envolvente es plana por construccion y el ajuste exponencial deja
    de ser informativo, mientras que la clasificacion sigue siendo nitida.

    Resultado empirico: el borde cae a r_0B/ell constante, con
    r_0B = r0 (m/m_B) el radio clasico desnudo.
    """
    signo = 1.0 if rama == "positiva" else -1.0

    def clas(mr):
        u, d = _es_inestable(mr, ell, regulator, **kw)
        if verbose:
            print(f"      m/mB={mr:+8.4f}  inestable={u}  "
                  f"deriva={d['drift']:.1e}  a_max={d.get('a_max', 0):.2e}",
                  flush=True)
        return u

    # extremo estable: bajar hasta encontrar uno clasificado estable.
    # `guess` (tipicamente el borde del cutoff anterior) acorta la busqueda,
    # porque la frontera varia suavemente con el cutoff.
    a = 0.5 * guess if guess is not None else 0.05 * signo
    for _ in range(12):
        if clas(a) is False:
            break
        a *= 0.55
    else:
        return np.nan, {"motivo": "no se hallo extremo estable"}

    # extremo inestable: subir hasta encontrar uno clasificado inestable
    b = a
    for _ in range(20):
        b *= 1.6
        if abs(b) > hi:
            return np.nan, {"motivo": "sin borde por debajo del tope", "b": b}
        if clas(b) is True:
            break
    else:
        return np.nan, {"motivo": "no se hallo extremo inestable"}

    # biseccion sobre la clasificacion
    while abs(b - a) > tol * max(1.0, abs(b)):
        mid = 0.5 * (a + b)
        u = clas(mid)
        if u is None:                       # indeterminado: acercarse desde el lado estable
            a = mid
        elif u:
            b = mid
        else:
            a = mid

    mr = 0.5 * (a + b)
    tasa, d = relaxation_rate(Params(ell=ell, m_over_mB=mr), regulator, **kw)
    return float(mr), {"tasa_en_el_borde": tasa, **d}


def renormalized_ratio(ell: float, rama: str = "positiva",
                       regulator: str = "smeared", **kw):
    """R1: m/m_B donde la tasa no lineal iguala la de la teoria linealizada."""
    obj = linearized_rate(Params(ell=ell), regulator)
    if not np.isfinite(obj):
        return np.nan, {"motivo": "sin cero dominante en chi^r"}
    return solve_ratio(ell, obj, rama, regulator, **kw)
