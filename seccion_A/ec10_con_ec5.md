# Ec. (10) con el regulador de la ec. (5)

Polonyi, arXiv:1701.04068v4, p. 7. Punto 2 del plan: sustituir `δ_B(z)` de la
ec. (5) dentro de la ec. (10) y llevarlo hasta forma cerrada.

## Las dos piezas, verbatim

    δm = (e² / 2c²) ∫₀^∞ dz z^{-1/2} δ_B(z)                                (10)

    δ_B(x²) = Θ(x²) x² exp(−√(x²)/ℓ) / (12 ℓ⁴)                              (5)

es decir, con `z = x²` el intervalo invariante de Minkowski,

    δ_B(z) = Θ(z) · z · e^{−√z/ℓ} / (12 ℓ⁴).

## Sustitución

    δm = (e²/2c²) · (1/12ℓ⁴) ∫₀^∞ dz  z^{-1/2} · z · e^{−√z/ℓ}
       = (e²/2c²) · (1/12ℓ⁴) ∫₀^∞ dz  √z · e^{−√z/ℓ}

Cambio de variable `z = ℓ²t²`  (`t = √z/ℓ ≥ 0`, `dz = 2ℓ² t dt`, `√z = ℓt`):

    ∫₀^∞ dz √z e^{−√z/ℓ} = ∫₀^∞ (ℓt)(2ℓ² t) e^{−t} dt
                          = 2ℓ³ ∫₀^∞ t² e^{−t} dt
                          = 2ℓ³ · Γ(3) = 2ℓ³ · 2 = 4ℓ³.

Por lo tanto

    δm = (e²/2c²) · 4ℓ³ / (12ℓ⁴) = (e²/2c²) · 1/(3ℓ)

**Resultado:**

    δm = e² / (6 c² ℓ)

y con `e²/c² = m r₀` (definición del radio clásico del electrón, `r₀ = e²/mc²`):

    δm / m = r₀ / (6ℓ)                                                     (★)

## Consecuencias inmediatas

- `m = m_B + δm`  ⟹  `m_B/m = 1 − r₀/(6ℓ)` y `m/m_B = 1/(1 − r₀/(6ℓ))`.
  Es la **línea punteada de la Fig. 2(b)** del paper: la predicción de la teoría
  linealizada en el plano `(r₀/ℓ, m/m_B)`.
- `m_B` se anula en **`r₀/ℓ = 6`** y cambia de signo por encima. Ahí el
  contratérmino consume la masa entera y la ecuación renormalizada deja de
  existir (`integrate_linear` lo rechaza explícitamente).
- La divergencia es **lineal en el cutoff**, `δm ∼ Λ = 1/ℓ`, no cuadrática:
  es la divergencia clásica conocida de la autoenergía electrostática.

## Momento general (sirve para las tres condiciones a la vez)

El mismo cambio de variable da, para cualquier `p > −2`,

    ∫₀^∞ dz z^{p} δ_B(z) = (2 ℓ^{2p} / 12) ∫₀^∞ t^{2p+3} e^{−t} dt
                         = ℓ^{2p} Γ(2p+4) / 6

de donde salen de un golpe:

| p | integral | valor | significado |
|---|---|---|---|
| `0` | `∫dz δ_B(z)` | `Γ(4)/6 = 1` | **condición (i)**, normalización — el `12ℓ⁴` de la ec. (5) está puesto justo para esto |
| `−1/2` | `∫dz z^{-1/2}δ_B(z)` | `Γ(3)/(6ℓ) = 1/(3ℓ)` | **la ec. (10)**, ecuación (★) |
| `+1/2` | `∫dz √z δ_B(z)` | `ℓ Γ(5)/6 = 4ℓ` | la lectura *equivocada*, ver abajo |

## Advertencia sobre la lectura de la ec. (10)

Los extractores automáticos de texto leen la ec. (10) como `∫ dz √z δ_B(z)`
(porque el `1/√z` está compuesto como fracción). Esa lectura da

    δm = (e²/2c²) · 4ℓ  ⟹  δm/m = 2 r₀ ℓ,

que **crece** al remover el cutoff (`ℓ → 0` daría `δm → 0`): sería una
autoenergía que se apaga al acercar los puntos, lo contrario de lo que pasa.
Dimensionalmente, `δ_B` tiene dimensión `1/z` por la condición (i), de modo que
`∫dz z^{-1/2}δ_B ∼ 1/longitud` y `(e²/c²)·(1/longitud) = masa` ✓, mientras que
`∫dz z^{1/2}δ_B ∼ longitud` daría `masa · longitud²` ✗.

La lectura correcta es `dz/√z`, verificada contra la imagen de la p. 7.

**Ojo con un falso amigo.** La forma correcta *también* exhibe un `√z`, porque
`z^{-1/2} · z = z^{1/2}`:

    δm = (e²/24 ℓ⁴ c²) ∫₀^∞ dz √z e^{−√z/ℓ}          ← correcta,  = e²/(6c²ℓ)

de modo que ver un `√z` en el integrando no dice nada. Lo que separa las dos
lecturas es el exponente bajo el peso exponencial: `z^{1/2}` (correcta, `4ℓ³`,
momento `p = −1/2`) frente a `z^{3/2}` (errónea, `48ℓ⁵`, momento `p = +1/2`).

## Verificación numérica

`seccion_A/fig_ec7_ec10.py` imprime tres vías independientes al mismo número
(cuadratura de `scipy.quad`, la forma cerrada `1/(3ℓ)` y `Regulator.moment_inv_sqrt`
de `nlaid/core.py`), coincidiendo a 1e-12 para `ℓ = 0.7`:

```
  cuadratura            0.476190476191
  cerrada  1/(3 ell)    0.476190476190
  reg.moment_inv_sqrt   0.476190476190
  delta_m/m = r0/(6ell) 0.238095238095   (esperado 0.238095238095)
  lectura erronea dz*sqrt(z) = 2.800000  (= 4 ell, diverge al reves)
```

## Contraste con el regulador desplazado, ec. (4)

Para `δ_B(z) = δ(z − ℓ²)` la ec. (10) es inmediata: `∫dz z^{-1/2}δ(z−ℓ²) = 1/ℓ`,
luego `δm/m = r₀/(2ℓ)`, y `m_B` se anula en `r₀/ℓ = 2`. El suavizado da el mismo
comportamiento `1/ℓ` con coeficiente tres veces menor: la delta desplazada pone
todo el peso en `z = ℓ²`, mientras que la suavizada lo reparte sobre `z ≳ ℓ²`,
donde `z^{-1/2}` pesa menos. Es reparto del peso del regulador, no física
distinta — como debe ser, ya que `δm` es un contratérmino.
