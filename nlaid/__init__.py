"""Dinamica no local de una carga puntual con regularizacion point-splitting.

Implementacion de J. Polonyi, arXiv:1701.04068v4.

Bloque 1  block1_linear   ec. (14)  susceptibilidad y estabilidad linealizada
Bloque 2  block2_delay    ec. (16)  regulador desplazado, retardo finito
Bloque 3  block3_memory   ec. (17)  regulador suavizado, memoria infinita

Todos comparten `core.Params` y los reguladores de `core`.
"""

import pathlib

from .core import Params, ShiftedDelta, SmearedDelta, make_regulator, minkowski_dot

# Raiz del repo. Los scripts y el notebook la usan para localizar figures/ y
# data/ sin depender del directorio desde el que se los invoque. Valida con la
# instalacion editable (pip install -e .), que es como se usa este repo.
RAIZ = pathlib.Path(__file__).resolve().parent.parent

__all__ = ["RAIZ", "Params", "ShiftedDelta", "SmearedDelta",
           "make_regulator", "minkowski_dot"]
