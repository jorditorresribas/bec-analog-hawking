# Fase 2 — Pre-registro de la medición física (exponente de la anchura del horizonte ensanchado)

**Estado: PRE-REGISTRO. NO EJECUTADO.** Redactado el 2026-11-18. Diseño: `docs/FASE2_DISENO.md`.

**Condiciones previas para ejecutar** (las tres):
1. Que el **revisor externo (ChatGPT)** haya revisado el diseño y este pre-registro.
2. Que los gates técnicos pasen. **Hoy no pasan:** F0 y F1 fallan por un suelo de redondeo ∝ dx⁻⁴ de la formulación actual
   (`residual_anatomy/FASE2_GATES_RESULT.md`). Hace falta una formulación revisada con criterios F0′–F3′ en commit antes
   de ejecutarlos, y que **F0′–F3′ pasen todos**, incluido el de no degeneración.
3. Que el script de medición se escriba y se ponga en commit **antes** de ejecutarlo, con una comprobación de que este
   archivo no ha cambiado desde su commit (como I0 en E1-bis).

Si la formulación revisada cambia algo de lo de abajo (por ejemplo, las mallas), este pre-registro se sustituye por uno
nuevo, en commit, antes de ejecutar; no se edita después.

## 1. Pregunta
¿Cambia la **forma** de la dispersión UV, con la escala dispersiva fija, el exponente con que la anchura del horizonte
ensanchado depende de κ? Predicción por análisis de escalas (`docs/FASE2_DISENO.md` §3, no publicada):
d ∝ κ^{−1/3} para dispersión cuártica y d ∝ κ^{−1/5} cuando domina el término k⁶.

## 2. Qué se ha visto ya (declaración)
- De γ > 0: n_L en w = 2 con γ = 1 y 10 en unas pocas ω (depuración y gates técnicos) y casos homogéneos. **Ninguna
  temperatura, ningún κ̄, ningún d\* ni ningún ajuste de exponente.**
- De γ = 0: todo el benchmark v1.0 (`out/fp_benchmark.csv`) y la tabla de sensibilidad a d de `FP_d_sensitivity.py`.
  El control γ = 0 de abajo está, por tanto, **parcialmente visto**: con d = 1,58 el cociente va de 1,0050 (w = 2) a
  0,9926 (w = 8). Eso indica que d\* no es exactamente ∝ κ^{−1/3} en nuestros perfiles. **Por eso el control tiene su propio
  criterio, y la predicción principal se juzga frente al control medido, no frente a 1/3.**

## 3. Medición (fijada)
- Perfiles: tanh de v1.0 (Mach 0,5 → 2), **w ∈ {2, 3, 4, 6, 8}**; κ = κ_local = |(c − v)′|_H.
- **γ ∈ {0; 1; 10}**. γ = 0 es el control (debe reproducir el benchmark v1.0); γ = 10 es el caso k⁶ dominante (4γc_H² ≈ 25
  en p ≈ 1); γ = 1 es el cruce y **solo se informa**.
- T_up = T₀(n_L): exactamente el estimador del benchmark (banda ω ∈ [0,05; 0,27]·κ/2π, 9 puntos, T_ω = ω/ln(1 + 1/n_L),
  ajuste cuadrático, T₀ = término constante).
- d\*(w, γ): la semianchura de la ventana **simétrica** alrededor del horizonte para la que κ̄(d\*)/2π = T_up, con κ̄ la
  media de |(c − v)′| en la ventana (FP Ec. 20, misma implementación que `FP_d_sensitivity.py`). Si para algún w no existe
  solución (T_up fuera del rango de κ̄(d)/2π), ese (w, γ) se declara «sin d\*».
- Exponente: p(γ) por mínimos cuadrados de log d\* frente a log κ en los cinco w; error de p por *bootstrap* sobre los
  cinco puntos (1000 remuestreos, semilla 0).
- Malla y solve: los que fijen los gates F0′–F3′ revisados, con convergencia demostrada (Richardson en dx) para T_up a
  ≤ 10⁻⁴ relativo en cada (w, γ).

## 4. Criterios (fijados antes de ejecutar)
| Criterio | Condición | Papel |
|---|---|---|
| **K (control)** | γ = 0: T_up/(κ/2π) coincide con `out/fp_benchmark.csv` a ≤ 10⁻⁴ en los cinco w, y p(0) ∈ [0,25; 0,42] | si falla, la medición no tiene veredicto |
| **P1 (predicción principal)** | γ = 10: \|p(10) − 1/5\| ≤ 0,05 **y** p(0) − p(10) ≥ 0,07 | **CONFIRMA** la predicción |
| **P0 (refutación)** | γ = 10: \|p(10) − p(0)\| ≤ 0,03 | **REFUTA**: la forma de la dispersión no cambia el exponente en este alcance |
| otro caso | ni P1 ni P0 | **INCONCLUSO**: se publica igual, con la tabla |
| informativo | p(1) (cruce); d\*(w) para cada γ; T_up/(κ/2π) para cada (w, γ) | se informa sin criterio |

Además, para que haya veredicto: d\* debe existir en al menos 4 de las 5 anchuras para γ = 0 y γ = 10.

**Justificación de las tolerancias:** ±0,05 es el orden de la desviación que Finazzi y Parentani encontraron en su propio
ajuste del exponente cuártico (−0,627 frente a −2/3, un 5 %), trasladado a p. La separación mínima de 0,07 deja margen sobre
la dispersión del ajuste en cinco puntos. Ninguna de las dos se ha elegido mirando datos de γ > 0, porque no existen.

## 5. Qué se afirmará
- **CONFIRMA:** «En esta familia de perfiles, al hacer dominante un término k⁶ con la escala dispersiva fija, el exponente
  de la anchura del horizonte ensanchado cambia de ≈ 1/3 a ≈ 1/5, como predice el argumento de la fase». Es la única
  afirmación nueva: la generalización del argumento de Coutant y Parentani. No es física transplanckiana en sentido
  gravitatorio.
- **REFUTA / INCONCLUSO:** se publica igual, con la tabla, y se revisa el argumento del §3 del diseño.
