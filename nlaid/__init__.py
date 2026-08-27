"""Dinamica no local de una carga puntual con regularizacion point-splitting.

Implementacion de J. Polonyi, arXiv:1701.04068v4.

Bloque 1  block1_linear   ec. (14)  susceptibilidad y estabilidad linealizada
Bloque 2  block2_delay    ec. (16)  regulador desplazado, retardo finito
Bloque 3  block3_memory   ec. (17)  regulador suavizado, memoria infinita

Todos comparten `core.Params` y los reguladores de `core`.
"""

from .core import Params, ShiftedDelta, SmearedDelta, make_regulator, minkowski_dot

__all__ = ["Params", "ShiftedDelta", "SmearedDelta", "make_regulator", "minkowski_dot"]
