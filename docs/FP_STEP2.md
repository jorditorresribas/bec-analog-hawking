# Experimento B (Finazzi-Parentani) — PASO 2: κ̄ promediada (Ec.20) + T_IR(w)

> **⚠ ESTA "FALSACIÓN" NO APLICA A FP — RESUELTO EN PASO 4 (`FP_STEP4.md`, 2026-11-15).**
> La matemática de abajo (techo certificado, cota del valor medio, w≥4) es correcta como
> enunciado sobre T_IR_tot = T(n_tot). PERO el artículo (arXiv:1012.1556, Ap. B Ec. 42)
> define su temperatura Hawking como n_ω=|β_ω|² = el canal u UPSTREAM SOLO = nuestro n_L,
> NO la suma. Nuestro n_tot = n_L+n_R+ = la ocupación del COMPAÑERO n_{−ω} de FP, que no
> es su temperatura Hawking. Con el observable correcto (n_L), T_upstream<κ_local y
> =κ̄/2π a <1% → la simulación REPRODUCE a FP, no lo falsa. Paso 2 acotaba el observable
> equivocado. Ver `FP_STEP4.md`.

Sesión 2026-11-14 (cont.). Entregable: **construir κ̄ (horizonte ensanchado, Ec.20)
con la ventana d_ξ DERIVADA de nuestro fondo (Ec.24, sin importar d=1.58), medir
T_IR(w) con la maquinaria N_ω validada, y comparar T_IR vs κ_local/2π vs κ̄/2π.**
Construido sobre PASO 1 (`FP_SCALES.md`, escalas validadas a máquina) y sobre
`mode_matching.solve_inject` (el juez modal exacto, ya validado en H4).
Archivo: `FP_step2_kappabar.py`. CSV: `out/fp_step2_kappabar.csv`.

## 1. Definiciones (verificadas contra arXiv:1012.1556v2, HTML)
- **Ec.9 dispersión:** Ω² = c²k² + c⁴k⁴/Λ², con Λ = √2 c/ξ (ESCALA DE FRECUENCIA).
  En nuestras unidades (ℏ=m=g=1, λ=1): **Λ_FP = 2c²** (= c·k_*), idéntico a PASO 1.
- **Ec.20 κ̄** = (1/d_ξ)∫_{−d_ξ^L}^{+d_ξ^R} ∂ₓw_FP dx = [w_FP(d_ξ^R)−w_FP(−d_ξ^L)]/d_ξ.
  PROMEDIO UNIFORME del gradiente del flujo sobre una ventana que ENSANCHA el horizonte.
  FP usa w_FP = c+v (su flujo v<0); en nuestra convención (v>0, horizonte c=v) la
  combinación equivalente es **(c−v)**. Default honesto: split simétrico
  d_ξ^L=d_ξ^R=d_ξ/2 centrado en el horizonte interpolado.
- **Ec.24 d_ξ** = d·c₀/(κ_K·Λ²)^{1/3}, con **d≈1.58 AJUSTADO por FP a SUS perfiles
  (NO universal; Ec.26-27 es una regresión log-lineal).** Transferimos solo el
  ESCALADO (PASO 1 §3): aquí **NO fijamos d**; barremos d∈{1.0, 1.58, 2.0}.
- **κ_local** = |d(c−v)/dx|_H. `kappa_derive.py` demuestra SIMBÓLICAMENTE que la
  gravedad superficial covariante del acústico es EXACTAMENTE |c′−v′|_H (ratio=1,
  sin factor c_H). Es la gravedad superficial, no una candidata.

## 2. Geometría (λ=1). c_H y Λ_FP son IGUALES para toda w (sólo cambia la pendiente)
Horizonte interpolado (c=v); c_H=0.7937, Λ_FP=1.2599 en todas. Régimen **débilmente
dispersivo**: κ/Λ baja de 0.215 (w=2) a 0.054 (w=8).

| w | κ_local | κ_local/2π | κ/Λ | κ̄/κ_local (d=1.58) | techo max|f′|/κ_local |
|---|---|---|---|---|---|
| 2 | 0.2703 | 0.04302 | 0.215 | 0.9525 | 1.01434 |
| 3 | 0.1802 | 0.02868 | 0.143 | 0.9717 | 1.01434 |
| 4 | 0.1351 | 0.02151 | 0.107 | 0.9805 | 1.01434 |
| 6 | 0.0901 | 0.01434 | 0.072 | 0.9885 | 1.01434 |
| 8 | 0.0676 | 0.01075 | 0.054 | 0.9922 | 1.01434 |

- **κ̄ < κ_local siempre**, y más reducida en horizonte empinado (w chico): el
  promedio del gradiente SÍ reduce la gravedad efectiva, cualitativamente como FP.
  `FP_step2_audit.py` (B): en la **ventana simétrica centrada en x_H**, κ̄(d) es
  **monótona decreciente** desde κ_local (d→0⁺), verificada en el rango escaneado
  **d∈[10⁻³, 50]** (no una prueba analítica para todo d>0). Por tanto κ̄≤κ_local en ese
  rango. Independientemente del escaneo, κ̄ ≤ max|f′| para **TODO d** por el teorema del
  valor medio (abajo).
- **Techo CERTIFICADO:** max_x|d(c−v)/dx| = **1.014343·κ_local** (pendiente máxima
  1.43% río abajo del horizonte, invariante en w). Calculado con la **derivada
  analítica cerrada** de f=c−v + root-find de f″=0 (no `np.gradient`); `np.gradient`
  converge a él O(dx²) (error <2·10⁻⁶ ya con N=20000). **Máximo GLOBAL certificado:**
  |f′|=|n′|·B con n′∝sech²(x/w) → |f′|→0 como e^{−2|x|/w} en AMBOS infinitos; se
  escanean TODOS los puntos críticos de f″ en [x_H−30w, x_H+30w] (3 en todas las w) y
  |f′| en los bordes es 0 a máquina → el máximo hallado es global, no local. Como κ̄ es
  un promedio de |f′|, por el teorema del valor medio κ̄=|f′(x_*)| para algún x_* del
  intervalo, luego **ninguna κ̄ (cualquier centro/split/anchura) puede superarlo**.
  `FP_step2_audit.py` (A).

## 3. T_IR(w) medido INDEPENDIENTEMENTE (ω·N_ω → 0, estimador 2 de H4)
Banda baja om/Tscale∈[0.12,0.60], extrapolación lineal y cuadrática del límite IR
ω·N_ω→T. N_ω por matching modal exacto (`solve_inject`). **No se impone T_IR=κ/2π.**
Convergencia verificada (sonda aparte): T_IR/κ_local **estable a 4 cifras** al variar
caja y malla — w=2: 1.0093 con L=40,60,80 (N≤18000); w=8: 1.0383 con L=84,100,130
(N≤30000). Pseudo-unitariedad ≤2·10⁻³. Lin vs quad difieren 0.8% (robusto).

| w | T_IR (quad) | T_IR/κ_local | T_IR/κ̄(1.58) | T_IR/techo | ¿T_IR>techo? |
|---|---|---|---|---|---|
| 2 | 0.04342 | 1.0093 | 1.060 | 0.9950 | no |
| 3 | 0.02942 | 1.0259 | 1.056 | 1.0114 | **SÍ** |
| 4 | 0.02220 | 1.0320 | 1.053 | 1.0174 | **SÍ** |
| 6 | 0.01486 | 1.0366 | 1.049 | 1.0220 | **SÍ** |
| 8 | 0.01117 | 1.0383 | 1.047 | 1.0236 | **SÍ** |

## 4. Resultado (marco preciso, tras auditoría)
**La prescripción de gravedad superficial promediada (Ec.20) con el escalado de ventana
(Ec.24) NO explica el T_IR de esta familia de fondos; y para los horizontes suaves
(w≥4) NINGUNA ventana basada en el promedio del gradiente puede explicarlo.** (No se
afirma "FP falsado" en general: FP no reivindica universalidad de d≈1.58 ni de la
Ec.20, que ellos mismos declaran aproximada.)

1. **κ̄ < κ_local** en todas las ventanas examinadas (§3).
2. **T_IR por intervalos** (presupuesto §4b), NO una desigualdad uniforme:
   - **w=3,4,6,8:** T_IR/κ_local está ENTERO por encima de 1 (mín 1.012/1.018/1.023/
     1.025) → T_IR > κ_local robustamente → **no hay igualdad con una κ̄≤κ_local**
     (ventana simétrica centrada). FP predice T por DEBAJO; medimos por ENCIMA.
   - **w=2:** T_IR/κ_local ∈ [0.995, 1.009] **cruza 1** → el presupuesto NO permite
     decidir; podría existir un d>0 con κ̄=T_IR. **Queda indeterminado.**
3. **Independiente de convención para w≥4** (el resultado fuerte): T_IR supera el techo
   CERTIFICADO del gradiente (1.014343·κ_local) con TODA banda/orden (mín T_IR/techo =
   1.004/1.008/1.010). Como κ̄ ≤ max|f′| para TODO d y TODA ventana (teorema del valor
   medio), **ningún promediado del gradiente de este fondo puede alcanzar T_IR** ahí.
   Esto NO usa monotonía ni la convención de ventana. Para w=2–3 no se establece
   (el exceso cae dentro del margen del techo).
4. **La diferencia CRECE con w** (horizontes suaves/más hidrodinámicos), opuesto a la
   expectativa FP de que la dispersión importa más en horizontes empinados.

### 4b. Presupuesto de error de T_IR (respuesta a la auditoría). `FP_step2_audit.py` (C)
T_IR desde una malla densa de N_ω (om/T∈[0.08,1.1]), variando **banda** (6 ventanas) ×
**orden** de extrapolación (lin/quad/cúbica). Determinista (sin ruido estadístico);
"error" = sistemático de banda+orden. Dispersión ≈0.43% en toda w (≪ excesos de w≥4).

| w | T_IR/κ_local [min,max] | T_IR/techo [min,max] | ¿supera techo con TODA banda/orden? | Planck-fit ancho [0.15,1.1]T |
|---|---|---|---|---|
| 2 | [0.995, 1.009] | [0.981, 0.995] | **no** | 1.019·κ_local |
| 3 | [1.012, 1.026] | [0.998, 1.012] | no (a caballo) | 1.036·κ_local |
| 4 | [1.018, 1.032] | [1.004, 1.018] | **sí** | 1.042·κ_local |
| 6 | [1.023, 1.037] | [1.008, 1.022] | **sí** | 1.047·κ_local |
| 8 | [1.025, 1.038] | [1.010, 1.024] | **sí** | 1.049·κ_local |

Convergencia conjunta L,N (sonda aparte): T_IR/κ_local estable a 4 cifras (w=2: 1.0093
con L=40→80; w=8: 1.0383 con L=84→130, N≤30000). pseudo-u≤3·10⁻³, el exceso NO
correlaciona con pu (mejora con N, T_IR estable). El Planck-fit ancho da ratios
MAYORES que el límite IR → el exceso no es artefacto de baja ω.

### Corroboración independiente
El Tfit viejo (método transfer, `kappa_derive.py`) ya daba T/κ_geom = 1.062 (w=2) →
1.106 (w=8): **mismo signo y misma tendencia creciente con w**. Dos métodos
independientes (transfer y matching modal) coinciden en que la T acústica está por
ENCIMA de κ_local y el exceso crece con w. El rasgo es físico, no numérico.

## 5. Pregunta ABIERTA (honesta, no resuelta aquí)
El exceso T_IR/κ_local crece hacia el límite hidrodinámico (κ/Λ→0). Una corrección
DISPERSIVA debería DESAPARECER ahí, no crecer. Por tanto el +1..+4% **no es** una
corrección en κ/Λ: su origen (greybody/potencial de dispersión del horizonte suave,
o estructura no-térmica sub-dominante) queda PENDIENTE. Es la siguiente pregunta, no
la cierro a mano. Conecta con §7 (correlaciones/robustez, ya cerrada) — ver
`CORRELATIONS.md`.

## 6. Límites del alcance y pendientes (honestos)
- **Familia de perfiles ≠ la de FP.** FP construye fondos como componente de fondo +
  perturbación localizada y AJUSTA d_ξ^L, d_ξ^R a series de espectros; aquí usamos una
  familia parametrizada por w y aplicamos el escalado Ec.24 directamente. Por eso el
  resultado se enuncia sobre "esta familia", no como refutación de su procedimiento de
  ajuste. Reproducir sus perfiles exactos sería el control definitivo (no hecho).
- Reparto asimétrico d_ξ^L≠d_ξ^R y centro x̂ (Ec.20): no explorados. Irrelevantes para
  la falsación de w≥4 porque el TECHO certificado acota CUALQUIER ventana (cota del
  valor medio). Para w=2–3 sí podrían importar (exceso dentro del margen del techo).
- d* (d que reproduce T_IR): no existe para w=3,4,6,8 (T_IR>κ_local robusto, κ̄≤κ_local
  en d∈[10⁻³,50]); **indeterminado para w=2** (T_IR cruza κ_local). El escalado FP no
  transfiere a este régimen para los horizontes resueltos.
- Compatibilidad espectral: el ajuste de Planck en la banda [0.15,1.1]T es consistente
  (dispersión ~0.4%, sin oscilaciones visibles en ω·N). Es COMPATIBILIDAD con un ajuste
  térmico en esa banda, NO una prueba independiente de adiabaticidad ni de ausencia de
  amplificación paramétrica (el caso donde FP admite desviaciones queda sin verificar
  por separado).

Archivos: `FP_step2_kappabar.py` (producción), `FP_step2_audit.py` (controles A/B/C de
esta auditoría), `out/fp_step2_kappabar.csv`. Base: `FP_SCALES.md` (paso 1),
`mode_matching.py` (N_ω), `kappa_derive.py` (κ covariante = |c′−v′|).
