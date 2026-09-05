# La ecuación linealizada en el tiempo, y qué es (y qué no es) markoviano

## Una distinción que hay que hacer explícita

**Linealizar no es hacer markoviano.** Son dos aproximaciones independientes:

| | qué se descarta | qué queda |
|---|---|---|
| **Linealizar** | los términos O(ξ²) | la memoria **entera**: ec. (14), no local |
| **Límite local (markoviano)** | la memoria, salvo el orden dominante | ec. (15), Abraham-Lorentz |

Y dan respuestas **opuestas** sobre estabilidad: la ec. (15) tiene el polo runaway
en ω = 3i/2 para cualquier cutoff, mientras que la ec. (14) es estable siempre que
r₀/ℓ esté por debajo del crítico (4.00 suavizado, 0.669 desplazado).

## Derivación

Perturbando la carga en reposo, x = (s, ξ(s)) con |ξ| ≪ r₀, a primer orden:

    d = x - x' = (-u, Δ),   u = s' - s ≤ 0,   Δ = ξ(s) - ξ(s+u)
    w = d² = u² + O(ξ²)
    V = (0, Δ + u ξ̇(s+u))          ← sólo componente espacial

Las dos componentes de la ec. (6) colapsan a **una** ecuación integro-diferencial
lineal y escalar:

    ξ̈(s) = 4 r_0B ∫_{-∞}^{0} du δ_B'(u²) [ ξ(s) - ξ(s+u) + u ξ̇(s+u) ]      (*)

Con ξ = e^{-iωs} sale la relación de dispersión

    D(ω) = ω² + 4 r_0B ∫ du δ_B'(u²) [ 1 - (1 + iφ) e^{-iφ} ] ,   φ = ω u

## El contratérmino de masa aparece solo

Desarrollando f(φ) = 1 - (1+iφ)e^{-iφ} = -φ²/2 + iφ³/3 + φ⁴/8 + … e integrando
por partes cada momento (usando δ_B(0) = 0, condición (ii) del regulador):

    ∫ du δ_B'(u²) u² = -½ ∫ du δ_B(u²) ≡ -M₀/2        M₀ = δm/(m r₀)
    ∫ du δ_B'(u²) u³ = -M₁ = +1/2                     (por normalización)
    ∫ du δ_B'(u²) u⁴ = -(3/2) I₂

resulta

    D(ω) = ω²(1 + r_0B M₀) + (2/3) i r₀ ω³ - (3/4) r₀ ω⁴ I₂ + …

Los órdenes ω³ y ω⁴ coinciden **exactamente** con ω²χ de la ec. (14). El único
sobrante es `r_0B M₀ ω²`, y ése es el contratérmino: **el término -φ²/2 de f(φ)
es δm**. Restándolo,

    D_ren(ω) = ω² χ^r_ω        verificado numéricamente a 1e-9 a 4e-10

(`tests/test_block1.py::test_dispersion_es_omega2_por_chi`). En el tiempo, la
resta pasa δm ξ̈ al miembro izquierdo y convierte m ξ̈ en m_B ξ̈.

## Validación cruzada

`integrate_linear` marcha (*) hacia adelante con los nodos de cuadratura
**sobre la malla** (u_j = -j·ds), de modo que la memoria es una combinación
lineal de los últimos N valores y no hace falta interpolar la historia.
Prehistoria ξ = ξ̇ = 0 (solución exacta), excitación por la fuente k(s) de la
ec. (2).

Tasa y periodo de la trayectoria integrada, contra Im ω y 2π/Re ω del cero
dominante de D — dos caminos numéricos sin nada en común:

| r₀/ℓ | tasa integrada | Im ω espectral | periodo integrado | 2π/Re ω |
|---|---|---|---|---|
| 2.0 | −0.8284 | −0.8346 | 4.3345 | 4.3626 |
| 3.0 | −0.7122 | −0.7179 | 2.4328 | 2.4414 |
| 3.5 | −0.4631 | −0.4623 | 1.9361 | 1.9360 |
| **4.0** | **−0.0002** | **−0.0000** | 1.5709 | 1.5708 |
| 4.5 | +0.7977 | +0.7994 | 1.2961 | 1.2961 |
| 5.0 | +2.2443 | +2.2637 | 1.0890 | 1.0891 |

El periodo coincide a 4 cifras. La fila r₀/ℓ = 4.00 da tasa −0.0002: es el
punto marginal, y cae sobre el `critical_cutoff` = 4.00 del bloque 1, obtenido
por principio del argumento sin integrar ninguna trayectoria.

## Límites de validez

- **Sólo regulador suavizado.** δ_B' de la ec. (4) es distribucional y no admite
  evaluación puntual; esa rama es la ec. (16), en `block2_delay`.
- **r₀/ℓ = 6 no existe.** Ahí δm = m, o sea m_B = 0: el contratérmino diverge y
  la ecuación renormalizada no está definida. `integrate_linear` lo rechaza con
  un mensaje explícito en vez de dividir por cero.
- **Estabilidad significa ξ̈ → 0, no ξ → 0.** Tras el pulso la carga queda con
  velocidad constante y ξ crece linealmente: la línea de mundo se vuelve recta,
  que es lo que "relajar" quiere decir aquí.

## Observación sobre el modo marginal (no es una ley general)

En el cutoff crítico hay **un solo** modo marginal (χ = 0 con ω real). Su
frecuencia, en unidades del cutoff:

| regulador | r₀/ℓ crítico | ω_c | **ω_c · ℓ** | \|χ\| en el mínimo |
|---|---|---|---|---|
| suavizado ec. (5) | 4.000011 | 4.0008 | **1.00021** | 1.5e-4 |
| desplazado ec. (4) | 0.669430 | 10.4724 | **15.644** | 1.0e-4 |

Para el suavizado el modo marginal oscila exactamente con el período del
cutoff, ω_c ℓ = 1 (el residuo 2e-4 es compatible con la tolerancia 1e-5 de la
bisección sobre r₀/ℓ). **Para el desplazado no**: ω_c ℓ = 15.6. No es una
propiedad general del point-splitting, sino del núcleo suavizado, cuyo
decaimiento invariante e^{-√z/ℓ} fija la misma escala. No se ha demostrado el
mecanismo; queda como observación numérica.

*(Cuidado: buscar el modo marginal como cero de Re χ sobre el eje real es
incorrecto — da raíces espurias. Hay que minimizar |χ|.)*
