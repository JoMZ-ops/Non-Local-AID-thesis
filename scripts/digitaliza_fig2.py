"""Digitaliza la Fig. 2 del paper leyendo el PDF como VECTORES, no como pixeles.

Las figuras de Polonyi son graficos vectoriales (Mathematica) embebidos en el
PDF: el content stream trae las coordenadas exactas de cada trazo. Se
interpreta ese stream -- matriz de transformacion, pila q/Q, patron de guiones
y texto -- y se convierten los trazos a coordenadas de datos con la
calibracion de las marcas de eje.

Frente a la digitalizacion por pixeles de `data/fig2_digitalizada.npz`, esto no
tiene error de lectura: son los numeros que Mathematica escribio.

Calibracion (identica en los dos paneles, verificada sobre las marcas):
    ejes cruzan en (33.7, 67.0) del espacio del PDF
    8.7222 unidades de PDF por unidad de r0/ell
    53.9  unidades de PDF por unidad de m/m_B

Salida: data/fig2_vectorial.npz con, por panel, la region visible, las dos
ramas punteadas y las dos solidas.

El PDF NO esta en el repo. Se pasa como argumento:

    python3 scripts/digitaliza_fig2.py ~/Downloads/1701.04068v4.pdf
"""

import re
import sys
import zlib

import numpy as np

from nlaid import RAIZ

SALIDA = RAIZ / "data/fig2_vectorial.npz"
PANELES = {"a": 16, "b": 19}          # indice del stream de cada panel
X0, SX, Y0, SY = 33.7, 8.7222, 67.0, 53.9


# --------------------------------------------------------------------------
# Interprete de content streams
# --------------------------------------------------------------------------

def _mul(a, b):
    """Producto de matrices PDF [a b c d e f]."""
    return (a[0]*b[0] + a[1]*b[2], a[0]*b[1] + a[1]*b[3],
            a[2]*b[0] + a[3]*b[2], a[2]*b[1] + a[3]*b[3],
            a[4]*b[0] + a[5]*b[2] + b[4], a[4]*b[1] + a[5]*b[3] + b[5])


def _ap(m, x, y):
    return (m[0]*x + m[2]*y + m[4], m[1]*x + m[3]*y + m[5])


def parse_stream(txt):
    """Trazos del content stream, ya en coordenadas de pagina.

    Devuelve dicts {pts, dash, w, op}. El `dash` es lo que distingue la linea
    punteada de la solida, y es el dato que resuelve F1.
    """
    ctm, pila = (1, 0, 0, 1, 0, 0), []
    dash, ancho = "", 1.0
    pts, cur, paths = [], None, []

    # Las cadenas (...) van PRIMERO en la alternancia: si no, un texto como
    # "(-)" -- el signo menos de una etiqueta de eje -- deja escapar su guion
    # y el tokenizador lo toma por un numero.
    tok = re.findall(r'\((?:[^()\\]|\\.)*\)|\[[^\]]*\]|/[A-Za-z0-9]+'
                     r'|[-+]?\d*\.?\d+|[A-Za-z*\'"]+', txt)
    num = []
    for t in tok:
        if re.fullmatch(r'[-+]?\d*\.?\d+', t):
            num.append(float(t))
            continue
        if t[0] in "[/(":
            num.append(t)
            continue
        n = [v for v in num if isinstance(v, float)]
        if t == "q":
            pila.append((ctm, dash, ancho))
        elif t == "Q" and pila:
            ctm, dash, ancho = pila.pop()
        elif t == "cm" and len(n) >= 6:
            ctm = _mul(tuple(n[-6:]), ctm)
        elif t == "d":
            arr = [v for v in num if isinstance(v, str) and v.startswith("[")]
            dash = arr[-1] if arr else ""
        elif t == "w" and n:
            ancho = n[-1]
        elif t == "m" and len(n) >= 2:
            if cur and len(cur) > 1:
                pts.append(cur)
            cur = [_ap(ctm, n[-2], n[-1])]
        elif t in ("l", "c") and len(n) >= 2 and cur is not None:
            cur.append(_ap(ctm, n[-2], n[-1]))
        elif t in ("S", "s", "f", "f*", "B", "b", "n"):
            if cur and len(cur) > 1:
                pts.append(cur)
            paths += [{"pts": p, "dash": dash, "w": ancho, "op": t} for p in pts]
            pts, cur = [], None
        num = []
    return paths


def a_datos(p):
    a = np.array(p["pts"])
    return np.column_stack([(a[:, 0] - X0) / SX, (a[:, 1] - Y0) / SY])


# --------------------------------------------------------------------------

def extrae(pdf):
    d = open(pdf, "rb").read()
    streams = re.findall(rb'stream\r?\n(.*?)endstream', d, re.S)
    out = {}
    for pan, idx in PANELES.items():
        paths = parse_stream(zlib.decompress(streams[idx]).decode("latin1"))

        for p in paths:                                  # area visible
            if p["op"] == "n" and len(p["pts"]) in (4, 5):
                c = a_datos(p)
                if 19 < np.ptp(c[:, 0]) < 22:
                    out[f"{pan}_area"] = np.array(
                        [c[:, 0].min(), c[:, 0].max(), c[:, 1].min(), c[:, 1].max()])
                    break

        pun = [p for p in paths if "0.8 3.2" in p["dash"] and len(p["pts"]) > 50]
        sol = [p for p in paths if p["dash"] == "[ ]" and p["w"] == 0.40
               and len(p["pts"]) > 100]
        for k, p in enumerate(sorted(pun, key=lambda q: a_datos(q)[:, 0].min())):
            out[f"{pan}_pun{k}"] = a_datos(p)
        for k, p in enumerate(sorted(sol, key=lambda q: -a_datos(q)[:, 1].mean())):
            out[f"{pan}_sol{k}"] = a_datos(p)          # sol0 = rama m/m_B > 0
    return out


def verifica(d):
    """Las tres comprobaciones que motivaron la digitalizacion."""
    print("\n1. Region visible de cada panel (eje vertical = m/m_B):")
    for pan in PANELES:
        a = d[f"{pan}_area"]
        print(f"   ({pan})  r0/ell in [{a[0]:.2f}, {a[1]:.2f}]   "
              f"m/m_B in [{a[2]:.3f}, {a[3]:.3f}]")
    print("   La Fig. 1 usa m/m_B = 1.95..2.00 y -3.80..-4.10: FUERA del eje.")

    print("\n2. Las punteadas de los dos paneles:")
    for k in (0, 1):
        a, b = d[f"a_pun{k}"], d[f"b_pun{k}"]
        igual = a.shape == b.shape and np.abs(a - b).max() == 0.0
        print(f"   rama {k}: {a.shape[0]} vs {b.shape[0]} puntos, "
              f"{'IDENTICAS (max|dif| = 0)' if igual else 'distintas'}")

    print("\n3. La punteada contra las dos leyes del contratermino:")
    c = d["a_pun1"]
    print(f"   {'r0/ell':>7} {'punteada':>10} {'1/(1-x/2) despl':>16} "
          f"{'1/(1-x/6) suav':>15}")
    for xt in (4, 6, 8, 12, 20):
        j = int(np.argmin(np.abs(c[:, 0] - xt)))
        x, y = c[j]
        print(f"   {x:7.2f} {y:10.4f} {1/(1-x/2):16.4f} {1/(1-x/6):15.4f}")


def main():
    if len(sys.argv) < 2:
        sys.exit(__doc__.strip().splitlines()[-1].strip())
    d = extrae(sys.argv[1])
    np.savez(SALIDA, **d)
    print(f"{SALIDA}  ({len(d)} arreglos)")
    verifica(d)


if __name__ == "__main__":
    main()
