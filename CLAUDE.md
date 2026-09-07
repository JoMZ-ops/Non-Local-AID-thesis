# Memoria del proyecto — Non-Local AID

## Referencia base

J. Polonyi, *The Abraham-Lorentz force and electrodynamics at the classical
electron radius*, [arXiv:1701.04068v4](https://arxiv.org/abs/1701.04068)
(v4, 4 jun 2019, 22 pp.). El PDF **no** está commiteado; las ecuaciones que
usa el repo están transcritas literalmente abajo y comentadas en
`docs/notas_paper.md`.

Tesis central: la fuerza de Abraham-Lorentz no es un defecto de la
electrodinámica clásica sino un **efecto de cutoff UV**. Al eliminar el campo EM
se forman loops de orden `O(ħ⁰)` — clásicos —, el loop diverge, y el regulador
debe ser Lorentz-invariante (*point-splitting*). La fuerza resultante es
cutoff-independiente pero **generada** por el cutoff: estructura de anomalía.

## Ecuaciones, verbatim del PDF

### Regulador — se reemplaza `δ(x²) → δ_B(x²)` en la Green retardada (p. 6)

Tres condiciones:

| # | condición | razón |
|---|---|---|
| (i) | `∫dz δ_B(z) = 1` | preserva el flujo del campo radiado |
| (ii) | `δ_B(0) = 0` | separa los puntos singulares de `Θ(x⁰)δ_B(x²)` |
| (iii) | `δ_B(z) = 0` para `z < 0` | **suprime la interacción superlumínica** |

    δ_B(x²) = δ(x² − ℓ²)                                                    (4)
    δ_B(x²) = Θ(x²) x² exp(−√(x²)/ℓ) / (12 ℓ⁴)                              (5)

La separación no local es el **intervalo invariante** `(x−x')² = ℓ²`, no `x²−x'²`.

### Ecuación de movimiento completa (p. 6)

    ẍ = 4 r_0B ∫_{−∞}^{s} ds' δ_B'((x−x')²)
        { (x−x')(ẋ·ẋ') − [ẋ·(x−x')] ẋ' }                                    (6)

con `r_0B = e²/m_B c²`, `δ'(z) = dδ(z)/dz`.

### Sección II.A — *Linearized equation of motion* (p. 7)

    ẍ = 4 r_0B ∫_{−∞}^{0} du δ'(u²) (x − x' + u ẋ')                         (7)

    F_r = 2 r_0B ∫_{−∞}^{0} du u² δ(u²) [ x − x' + u ẋ' − u²(ẍ' − ẍ/2)
                                          + (2u³/3) x⃛ ]                    (8)
    F_s = −r_0B ∫_{−∞}^{0} du δ(u²) ( ẍ + 4u³ x⃛/3 )

    F_s = −(r_0B/2) ẍ ∫_0^∞ dz z^{−1/2} δ_B(z) + (2/3) r_0B x⃛              (9)

    δm = (e²/2c²) ∫_0^∞ dz z^{−1/2} δ_B(z)                                 (10)

    F_r ≈ (2 r_0B ℓ / ℓ_x³) ∫ dũ δ̃_B(ũ²) ũ²,   ũ = u/ℓ,  δ̃_B = ℓ²δ_B      (11)

**Errata de extracción, importante:** los extractores de texto leen la (9) y la
(10) como `∫dz √z δ_B(z)`, porque el `1/√z` está compuesto como fracción. Es
dimensionalmente imposible (daría `δm ∝ ℓ`). La lectura correcta, verificada
contra la imagen de la p. 7, es `dz/√z`. Ver `seccion_A/ec10_con_ec5.md`.

### Resto (ver `docs/notas_paper.md` para el detalle)

- ec. (14) susceptibilidad `χ^r_ω`; ec. (15) límite local (Abraham-Lorentz);
- ec. (16) regulador desplazado, retardo finito; ec. (17) suavizado, memoria
  infinita; ec. (18) condición de renormalización.

## Resultados cerrados que usa el repo

| cantidad | ec. (4) desplazado | ec. (5) suavizado |
|---|---|---|
| `∫dz z^{−1/2}δ_B(z)` | `1/ℓ` | `1/(3ℓ)` |
| `δm/m` (ec. 10) | `r₀/(2ℓ)` | `r₀/(6ℓ)` |
| `m_B = 0` en | `r₀/ℓ = 2` | `r₀/ℓ = 6` |
| `I₂ = ∫du δ_B(u²)u²` | `ℓ/2` | `2ℓ` |
| `r₀/ℓ` crítico (estabilidad linealizada) | `0.669` | `4.00` |

Momento general del suavizado: `∫₀^∞ dz z^p δ_B(z) = ℓ^{2p} Γ(2p+4)/6`.

## Convenciones del repo

- Signatura `(+,−,−,−)`; líneas de mundo temporales `x² > 0`.
- `c = 1`, `r₀ = 1`: el cutoff `ℓ` se mide en unidades del radio clásico, y los
  ejes de la Fig. 2 son `r₀/ℓ` y `m/m_B`.
- Parametrización por tiempo propio `s`, con `ẋ² = 1` (consecuencia exacta de
  las ecs. 6/16/17, no una restricción impuesta: `WorldLine.norm_drift` mide su
  deriva como error de integración).
- `dim = 2` por defecto (1+1D), que es donde el paper produce sus figuras.

## Estructura

| ruta | qué hay |
|---|---|
| `nlaid/core.py` | métrica, `Params`, reguladores `ShiftedDelta` / `SmearedDelta` |
| `nlaid/worldline.py` | historia de línea de mundo, punto retardado, fuente externa |
| `nlaid/block1_linear.py` | ec. (14), dispersión, `integrate_linear` |
| `nlaid/block2_delay.py` | ec. (16) |
| `nlaid/block3_memory.py` | ec. (17) |
| `nlaid/block4_renorm.py` | ec. (18) |
| `scripts/` | figuras; **no es un paquete instalable** |
| `seccion_A/` | integrandos de las ecs. (7) y (10), y la sustitución (5)→(10) |
| `docs/` | notas de lectura y derivaciones a mano |

Salidas a `figures/` y `data/`, localizadas vía `nlaid.RAIZ`, nunca vía el
directorio de trabajo. Instalación: `pip install -e ".[dev]"`.

## Estilo de trabajo acordado

**Modo Ponytail.** Antes de escribir código: ¿tiene que existir? ¿lo resuelve
la stdlib? ¿lo resuelve algo que ya está en `nlaid/`? ¿cabe en una línea limpia?
Nada de wrappers, clases abstractas ni comentarios redundantes que no se hayan
pedido. Si algo queda a medias por diseño, se marca con `# ponytail: ...`.

- Dependencias: `numpy`, `scipy`, `matplotlib`. **No añadir más** (en particular
  sympy: las integrales de este proyecto salen a mano con `z = ℓ²t²`).
- Documentación, docstrings y comentarios **en español**.
- Cada resultado nuevo se valida contra una vía independiente (forma cerrada vs
  cuadratura vs otro módulo), no contra otra corrida del mismo código.
- Rama de desarrollo: `Electrodynamics---Abraham-&-Lorentz-Force`.
