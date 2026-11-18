# Fase 2 — Resultado de la medición pre-registrada (2026-11-18)

**Veredicto: SIN VEREDICTO — falla el control K.** Se informa tal como salió. No se ha cambiado ningún criterio ni se ha
repetido la ejecución.

Pre-registro: `FASE2_PREREG.md` (`870cb60`) + addenda 1–3. Script: `residual_anatomy/fase2/measure_F2.py` (`4a96879`),
sin modificar. Ejecución autorizada por el usuario tras la revisión externa, con I0 en OK (HEAD `df1faf5`).
Salidas: `out/fase2_medida_T.csv`, `out/fase2_medida_dstar.csv`, `out/fase2_medida_resultado.csv`,
`out/fase2_medida_manifest.json`, `out/fase2_medida_run.log` (159 s).

## Criterios
| Criterio | Resultado | ¿Cumple? |
|---|---|---|
| Convergencia de T_up (≤ 10⁻⁴) | 1·10⁻⁸ – 1,3·10⁻⁷ en los 15 casos | sí |
| d\* en ≥ 4 de 5 anchuras | 5/5 para γ = 0, 1 y 10 | sí |
| **K**, parte 1: T_up/(κ/2π) con γ = 0 frente al benchmark (≤ 10⁻⁴) | 3,2·10⁻⁵ | sí |
| **K**, parte 2: p(0) ∈ [0,25; 0,42] | **p(0) = 0,620** (16–84 %: 0,560–0,704) | **NO** |
| P1 / P0 | no se evalúan (falla K) | — |

## Datos (todos informativos, ninguno decide nada)
T_up/(κ/2π):

| w | γ = 0 | γ = 1 | γ = 10 |
|---|---|---|---|
| 2 | 0,957224 | 0,956481 | 0,939891 |
| 3 | 0,972997 | 0,972383 | 0,968957 |
| 4 | 0,978904 | 0,978519 | 0,976891 |
| 6 | 0,983241 | 0,983117 | 0,982414 |
| 8 | 0,984777 | 0,984732 | 0,984412 |

Anchura ℓ\* de la ventana con κ̄(ℓ\*)/2π = T_up (y en unidades de FP, d = ℓ\*(κΛ²)^{1/3}/c_H):

| w | γ = 0 | γ = 1 | γ = 10 |
|---|---|---|---|
| 2 | 1,573 (d = 1,49) | 1,587 (1,51) | 1,884 (1,79) |
| 3 | 1,857 (1,54) | 1,879 (1,56) | 1,996 (1,66) |
| 4 | 2,182 (1,65) | 2,202 (1,66) | 2,286 (1,72) |
| 6 | 2,910 (1,92) | 2,921 (1,92) | 2,982 (1,96) |
| 8 | 3,694 (2,21) | 3,700 (2,21) | 3,739 (2,24) |

Exponentes: p(0) = 0,620 [0,560; 0,704]; p(1) = 0,614 [0,555; 0,694]; p(10) = 0,508 [0,424; 0,644].

## Lectura (descriptiva; no cambia el veredicto)
- **Por qué falla K:** en esta familia de perfiles tanh, la anchura ℓ\* que hace κ̄(ℓ\*)/2π = T_up **no** escala como
  κ^{−1/3}: el «d» de Finazzi y Parentani en sus unidades no es constante, sino que crece de 1,49 (w = 2) a 2,21 (w = 8).
  Ya lo anticipaban la tabla de sensibilidad a d (`FP_BENCHMARK.md` §5b: con d = 1,58 el cociente va de 1,0050 a 0,9926) y
  el §2 del pre-registro. El exponente medido mezcla el ensanchamiento dispersivo con la geometría del perfil, cuya escala
  es w ∝ 1/κ (que daría p = 1). Así, **este observable no aísla el ensanchamiento dispersivo en nuestros perfiles**, y la
  pregunta del pre-registro no se puede contestar con él.
- **Lo que sí se ve:** el término k⁶ reduce T_up, y el efecto decrece al suavizar el horizonte (−1,8 % en w = 2, −0,04 % en
  w = 8, con γ = 10). Y ℓ\* crece, sobre todo en w = 2. Con γ = 1 el efecto es diez veces menor. Es coherente con F3′ (un
  efecto pequeño y parecido a un reescalado). Que p(10) sea menor que p(0) va en la dirección de la predicción, pero **no
  es significativo** (los intervalos se solapan) y, sin control válido, **no se afirma nada**.
- **Qué no se afirma:** ni que el exponente cambie a 1/5 ni que no cambie. La fase 2 deja un instrumento validado (gates
  F0′–F3′) y una medición sin veredicto.

## Si se quisiera seguir (no hecho; requeriría un pre-registro nuevo)
Un observable que separe la geometría del perfil del ensanchamiento dispersivo: por ejemplo, perfiles con κ variable a
anchura relativa fija, o el déficit 1 − T_up/(κ/2π) frente a κ en lugar de la ventana de κ̄. Cualquiera de ellos se
diseñaría y pre-registraría desde cero, **sin** usar los números de esta tabla para elegirlo de forma que salga.
