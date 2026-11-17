# Analog Hawking radiation in a BEC acoustic black hole — reproducibility record

> **Nota de la versión pública (2026-11-17).** Este documento se escribió con todos los archivos en una sola carpeta. En el repositorio público están reorganizados: solver en `solver/`, benchmark FP en `benchmarks/fp_v1/`, cadena de validación en `methods/`, correlaciones en `correlations/`, documentos en `docs/`; `out/`, `reference/` y `residual_anatomy/` conservan su ruta. Todos los scripts se ejecutan **desde la raíz del repositorio** (p. ej. `python3 benchmarks/fp_v1/FP_check_reproduction.py`). Los SHA-256 citados aquí corresponden a la estructura original; los de la versión pública están en `docs/provenance/RELEASE_MANIFEST.txt`.

Sesión 2026-11-14. Registro primario para reproducibilidad y base de la sección
*Methods and Numerical Validation*. **Nivel de ambición (preservar exactamente):**
Hawking **análogo** en un BEC. NO "demostración de radiación de Hawking".

Unidades: ħ = m = 1 (m = 1/2 en la GPE adimensional estándar; aquí el operador
cinético es −½∂²). g = 1, n_u = 1.

---

## 0. Matriz de afirmaciones (qué se puede decir)

| Afirmación | Nivel |
|---|---|
| Conversión de modos de Bogoliubov (norma +/−) | **Demostrado numéricamente** |
| Concordancia β y N_ω entre dos métodos (transfer vs mode-matching) | **Demostrado** (dentro de tolerancias: β 1e-4, N_ω 1e-5) |
| Espectro aproximadamente térmico en el IR | **Respaldado por datos** |
| T_IR compatible con κ_geom/2π (w=2,8) | **Respaldado** (ratios 0.99–1.04) |
| Igualdad exacta T = κ/2π | **NO demostrado** |
| Incertidumbre total = ±2% | **NO demostrado** (es sensibilidad metodológica, no intervalo de confianza) |
| Residuo w=8 (~3%) no es dispersión | **NO demostrado** (argumento de escala, no excluye cuantitativamente; además ~1.5σ) |
| El ~10% inicial era "universalmente" artefacto | **Demasiado general**; identificados efectos de ventana y tamaño finito en los casos estudiados |

---

## 1. Observables (definición exacta + normalización)

**Fondo (waterfall):** GPE estacionaria por método inverso. Se elige
n(x) = n_u + (n_d−n_u)·½(1+tanh(x/w)); v(x) = j/n (continuidad, j = n_u·v_u);
V(x) = μ − g·n − ½v² + ½[ (√n)''/√n ] (despejado de la GPE). Mach_u = 0.5
(v_u = 0.5, c_u = 1), Mach_d = 2.0 → n_d = (j/Mach_d/√g)^(2/3) = 0.39685,
j = 0.5, μ = 1.125, n_H = (j²/g)^(1/3) = 0.62996.

**Modos de Bogoliubov:** dispersión asintótica (ω−vk)² = c_s²k² + k⁴/4,
c_s = √(gn). Norma de Bogoliubov: sign(N) = sign(Ω), Ω = ω−vk. Velocidad de
grupo v_g = v + (c_s²k + k³/2)/Ω. Corriente de norma conservada (derivada y
verificada J(x)=cte a 6 cifras): J = Im(u*u') + v|u|² + Im(v*v') − v|v|².
**Normalización de flujo:** cada canal propagante a |J| = 1.

**β (elemento anómalo):** |β|² = |S[neg-out, pos-in]|² (canal entrante upstream
norma +, saliente downstream norma −). Transfer: elemento de la matriz S.
Mode-matching: amplitud saliente norma − al inyectar el canal entrante norma +.

**N_ω (emisión térmica, el observable que diverge en IR):**
N_ω = Σ_{j: η_in,j<0} |S_{i,j}|² sobre canales OUT de norma +. Mode-matching:
inyectar el canal entrante de norma NEGATIVA (downstream), sumar |amplitud|² de
los salientes de norma +. **β ≠ N_ω** (β~0.01, N_ω~0.12 en ω=0.1).

**Pseudo-unitaridad (métrica de NORMA, no de corriente):**
Σ_out η_norma |a|² = η_in. (Error pasado: usar sign(corriente/v_g) en vez de
sign(norma) → falso déficit de "flujo" del 2%; corregido.)

**κ geométrica (derivada de la métrica, NO asumida):** métrica acústica
ds² = −(c_s²−v²)dt² − 2v dt dx + dx²; ξ = ∂_t; κ² = −½(∇_μξ_ν)(∇^μξ^ν)|_H.
Resultado simbólico (sympy): κ² = (c_s'−v')² → **κ = |d(c_s−v)/dx|_H** exacto.
T_grav = κ/2π. (NO lleva factor ½ ni c_H.)

---

## 2. Parámetros numéricos por configuración

| Config | L | N (dx) | método |
|---|---|---|---|
| β / N_ω validación (ω=0.1) | 24 | 4800 (0.01) | mode-matching; transfer X0=12 |
| Barrido de banda | 24 | 4000 | mode-matching + transfer |
| Convergencia hidro w=2 | 40 | 8000 | mode-matching |
| Convergencia hidro w=8 | **80** (NO 34) | 16000 | mode-matching |
| H4 IR w=2 | 60 | 12000 | mode-matching |
| H4 IR w=8 | 100 | 20000 | mode-matching |

**CRÍTICO (H3):** w=8 con L=34 (el usado inicialmente en hydro_convergence,
L=max(22,3.2w+8)=34) daba pseudo-unit 2.8e-2 y N_ω mal por TAMAÑO FINITO. Con
L≥50 se limpia (pseudo-unit 1e-3). Usar L≥80 para w=8 a baja ω.

---

## 3. Criterios de convergencia

- **Operador FD:** |M_BdG Ψ − ωΨ|/|ωΨ| ~ dx² (5.8e-7 a N=4000). Validado homogéneo.
- **β / N_ω:** convergencia en N (dN/N < 3e-4) Y en L (dN/N < 1e-4 a L convergido).
- **Pseudo-unitaridad por frecuencia:** < 1e-2 para fiar; típico mode-matching ~1e-4.
- **Fit Planck:** ventana escalada a T_grav (ω ~ 0.5–3·Tg), NO ventana fija.
- **Extrapolación IR:** T_IR = lim_{ω→0} ω·N_ω (lineal y cuadrática); leave-one-out
  spread ~0.9%; variar ventana superior ~1%.

---

## 4. Scripts → salidas (cada cifra tiene su origen)

| Script | Qué produce |
|---|---|
| `gpe.py` + `validate.py` | solver GPE (solitón 1e-7, Thomas-Fermi 0.4%) |
| `waterfall.py` | fondo waterfall autoconsistente → `out/waterfall_state.npz` |
| `evolve_fd.py` | evolución temporal CN + mapa Δn(x,t) |
| `bdg_solver/analysis/cap/cap_audit.py` | estabilidad (artefacto de pared, CAP) |
| `bdg_channels.py` | conteo de canales (upstream 2/0−, downstream 2/2−, ω_max=0.468) |
| `bdg_first_order.py` | sistema 1er orden dim-4 (char.poly = cuártico, 2e-15) |
| `bdg_smatrix.py` | matriz S por transfer (pseudo-unit 1e-4) |
| `bdg_hawking_spectrum.py` | espectro N_ω (transfer) |
| `mode_matching.py` | JUEZ independiente: matching modal + `solve_inject` (β y N_ω) |
| `fd_bdg_operator.py` | operador BdG FD disperso (base del juez) |
| `direct_scattering.py` | intento v1 (CAP) — FALLÓ convergencia en L (escenario C) |
| `convergence_study.py` | barrido anchura transfer — murió en w≥4 (conditioning) |
| `sweep_band.py` | β MM vs transfer en toda la banda → `out/*.csv` |
| `kappa_derive.py` | κ geométrica simbólica = |d(c−v)/dx| |
| `hydro_convergence.py` | T_fit/T_grav vs w |
| `temperature_window_audit.py` (H1) | sensibilidad a la ventana → `out/temperature_window_audit.csv` |
| `h3_convergence_audit.py` (H3) | convergencia L/N → `out/h3_convergence.csv` |
| `h4_ir_extrapolation.py` (H4) | límite IR, dos estimadores → `out/h4_ir.csv` |
| `h4_1_robustness.py` (H4.1) | leave-one-out, ventana, N-conv |

---

## 5. Independencia de las pruebas (honestidad — ChatGPT)

- **Independiente de verdad:** transfer (RK4, propagación) vs mode-matching (FD,
  BVP disperso) → implementaciones distintas, misma normalización física.
  Su acuerdo en β y N_ω es corroboración cruzada fuerte.
- **NO independientes entre sí:** pseudo-unitaridad y |α|²−|β|²=1 (misma identidad).
  Los dos estimadores de T (Planck y ω·N→0) usan los MISMOS datos N_ω → pueden
  compartir error sistemático de extrapolación. H4, H4.1 usan solo mode-matching.
- **Caveat global:** ambos solvers comparten el operador BdG, las raíces
  asintóticas y las convenciones de norma → el acuerdo no excluye un error
  sistemático COMÚN.

---

## 6. Resultado final (la frase defendible)

> Se construye un flujo de BEC con horizonte sónico autoconsistente. La conversión
> de modos de Bogoliubov (β) y la emisión (N_ω) se validan con dos métodos
> numéricos independientes. El espectro es aproximadamente térmico en el régimen
> infrarrojo, y su temperatura es compatible con la gravedad superficial acústica
> geométrica (derivada de la métrica): T_IR/T_grav ≈ 1.00 (w=2) y ≈ 1.03 (w=8),
> con sensibilidad metodológica ~2%, tras identificar y controlar sesgos de
> ventana de ajuste y de tamaño finito del dominio. No se demuestra igualdad
> exacta ni se excluye cuantitativamente una corrección dispersiva.

Entorno: Python 3.9.6, numpy 2.0.2, scipy 1.13.1, sympy, matplotlib 3.9.4; macOS
Darwin, Intel i7-4870HQ.
