"""Barrido del borde de estabilidad (condicion R2 del bloque 4).

Para cada cutoff resuelve por shooting el m/m_B donde la tasa de relajacion
se anula, en las dos ramas (m_B > 0 y m_B < 0). El resultado es la frontera
de la region sombreada de la Fig. 2 del paper.

Resultado empirico: la frontera cae a r_0B/ell constante, donde
r_0B = r0 (m/m_B) es el radio clasico DESNUDO. Es decir el criterio de
estabilidad compara el acoplamiento desnudo con el cutoff, no con r0.

Uso:  python3 scripts/scan_borde.py [suavizado|desplazado]
Salida: data/borde_<regulador>.json
"""

import json
import sys
import time
import warnings

import numpy as np

from nlaid.block4_renorm import stability_edge

warnings.simplefilter("ignore")

AJUSTES = {
    "smeared": dict(s_end=14.0, ds=1e-2),
    "shifted": dict(s_end=10.0, ds=5e-3),
}
CUTOFFS = [2.0, 3.0, 4.0, 6.0, 8.0, 12.0]


def main():
    reg = "smeared"
    if len(sys.argv) > 1:
        reg = {"suavizado": "smeared", "desplazado": "shifted"}.get(
            sys.argv[1], sys.argv[1])

    kw = AJUSTES[reg]
    print(f"Borde de estabilidad -- regulador {reg} -- {kw}")
    print(f"{'r0/ell':>7} {'rama +':>10} {'rama -':>10} | "
          f"{'x(m/mB)+':>9} {'x(m/mB)-':>9}  tiempo")

    res, semilla = {}, {}
    for x in CUTOFFS:
        ell, fila, t0 = 1.0 / x, {}, time.time()
        for rama in ("positiva", "negativa"):
            try:
                # La frontera cae a r_0B/ell constante, asi que el borde del
                # cutoff anterior escalado por la razon de cutoffs predice el
                # siguiente y ahorra buena parte del bracketing.
                mr, d = stability_edge(ell, rama=rama, regulator=reg,
                                       guess=semilla.get(rama), **kw)
                if np.isfinite(mr):
                    semilla[rama] = mr
            except Exception as exc:                       # noqa: BLE001
                mr, d = np.nan, {"motivo": type(exc).__name__}
            fila[rama] = float(mr) if np.isfinite(mr) else None
            if not np.isfinite(mr):
                fila[rama + "_motivo"] = str(d.get("motivo", ""))[:60]
        res[x] = fila

        f = lambda v: v if v is not None else float("nan")   # noqa: E731
        p, n = f(fila["positiva"]), f(fila["negativa"])
        print(f"{x:7.1f} {p:10.4f} {n:10.4f} | {x*p:9.3f} {x*n:9.3f}  "
              f"{time.time()-t0:.0f}s", flush=True)
        json.dump(res, open(f"data/borde_{reg}.json", "w"), indent=1)

    prods = [x * v["positiva"] for x, v in res.items() if v["positiva"]]
    prodn = [x * v["negativa"] for x, v in res.items() if v["negativa"]]
    if prods:
        print(f"\n  r_0B/ell en la rama +: {np.mean(prods):.3f} "
              f"+- {np.std(prods):.3f}")
    if prodn:
        print(f"  r_0B/ell en la rama -: {np.mean(prodn):.3f} "
              f"+- {np.std(prodn):.3f}")


if __name__ == "__main__":
    main()
