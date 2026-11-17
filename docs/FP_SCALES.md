# Experimento B (Finazzi-Parentani) — PASO 1: escalas dispersivas (derivación + validación)

Sesión 2026-11-14 (cont.). Entregable: **fijar las escalas en nuestras unidades y
validarlas con derivación + aserciones**, antes de construir κ̄ (paso 2). Referencia:
Finazzi-Parentani, PRD 83 084010 (2011), arXiv:1012.1556.
**Revisado tras auditoría de ChatGPT (2026-11-14): corregidos el mapeo Λ_FP, las dos
convenciones de ξ, el horizonte interpolado y las aserciones del script.**

## 1. Derivación

Dispersión validada (ℏ=m=g=1), con knob λ:
```
Omega^2 = (omega - v k)^2 = c^2 k^2 + (lambda/4) k^4 ,   c^2 = g n(x)
```
En λ=1 es EXACTAMENTE la dispersión de Bogoliubov: Ω²=(k²/2)²+(2gn)(k²/2)=k⁴/4+c²k²
(el término k⁴ = presión cuántica (k²/2)²). λ escala la presión cuántica
(masa efectiva m_eff=1/√λ; acoplo efectivo G=gn/√λ; coef. cinético t=√λ/2 en
`build_op_lambda`).

### 1a. Longitudes — DOS convenciones de curación (corrección auditoría §1)
- **Número de onda de cruce dispersivo** (k⁴ iguala a c²k²): de (λ/4)k_*⁴=c²k_*² →
  **k_* = 2c/√λ**,  longitud de cruce **ℓ = 1/k_* = √λ/(2c)**.
- **Curación EFECTIVA** del modelo modificado (con m_eff, G efectivos):
  **ξ_λ = √λ/(√2 c)**  →  **ξ_λ/ℓ = √2** (independiente de λ). Es la que debe usarse al
  hablar de "la curación del modelo con dispersión modificada".
- **Curación de REFERENCIA** (la del modelo original λ=1, mantenida fija):
  **ξ_ref = 1/(√2 c)**  →  **ξ_ref/ℓ = √(2/λ)**. En λ=1 coincide con ξ_λ.
- Ambas son legítimas pero NO intercambiables. Reservamos "longitud de curación" sin
  calificativo para **ξ_λ** (la del modelo efectivo); ξ_ref se nombra explícitamente.

### 1b. Mapeo a FP — Λ_FP es una ESCALA DE FRECUENCIA (corrección auditoría §2)
FP, Ec. (9): `Omega^2 = c^2 k^2 + (ħ²k⁴)/(4m²) = c^2 k^2 + c^4 k^4/Λ_FP^2`, y definen
Λ_FP como **escala de FRECUENCIA**, ξ_FP = √2 c/Λ_FP. Igualando a la nuestra (ħ=m=1):
```
c^4/Λ_FP^2 = lambda/4   =>   Λ_FP = 2 c^2 / sqrt(lambda) = c · k_*
```
→ **Λ_FP = c·k_*** (frecuencia), NO k_* (número de onda): difieren en un factor c.
Forma FP de la dispersión: **Ω² = c²k² (1 + c²k²/Λ_FP²)**. Con esto, ξ_FP = √2 c₀/Λ_FP;
con c₀ = c local da **ξ_FP = √λ/(√2 c) = ξ_λ** (la curación EFECTIVA), que coincide con
1/(√2 c) SOLO en λ=1 (corrección auditoría §2: no confundir con ξ_ref para λ≠1).
(La afirmación previa "Λ = k_*" era una confusión dimensional; k_*=1.587 en el horizonte
sigue bien como número de onda, pero NO es la Λ del paper.)

## 2. Validación numérica (contra el solver, bordes de caja + HORIZONTE interpolado)

`FP_step1_scales.py` (con aserciones). Fondo w=2. Horizonte = raíz INTERPOLADA de
c−v=0: x_H=0.46209, c_H=v_H=0.79370 (|c−v|=4.9e-7, no el punto de malla).

Escalas (λ=1):

| punto | c | v | ξ_ref | ξ_λ | ℓ | k_* | Λ_FP |
|---|---|---|---|---|---|---|---|
| borde izq. caja (x=−60) | 1.0000 | 0.5000 | 0.7071 | 0.7071 | 0.5000 | 2.0000 | 2.0000 |
| horizonte interp. (c=v) | 0.7937 | 0.7937 | 0.8909 | 0.8909 | 0.6300 | 1.5874 | 1.2599 |
| borde der. caja (x=+60) | 0.6300 | 1.2599 | 1.1225 | 1.1225 | 0.7937 | 1.2599 | 0.7937 |

(ξ_ref=ξ_λ en λ=1, como debe.) Nota: upstream/downstream son BORDES de la caja
numérica, no los límites matemáticos x→±∞.

**Aserciones (λ=0.5, 1, 2), todas PASAN:**
- ξ_λ/ℓ = √2 (exacto, indep. de λ).
- ξ_ref/ℓ = √(2/λ) (2.0 / 1.414 / 1.0).
- Λ_FP/(c·k_*) = 1 (a 1e-12).
- Forma FP c²k²(1+c²k²/Λ_FP²) = c²k²+λk⁴/4 (err ≤ 6e-17).
- Dispersión del solver vs analítica (bordes + horizonte): max err **1.8e-14** (24 raíces).
- Símbolo de `build_op_lambda`, autovalor físico ~ ω: max err **≤2.4e-15**.

## 3. Qué queda DETERMINADO y qué es CONVENCIÓN / PENDIENTE

**Determinado sin ambigüedad por nuestra dispersión:**
- k_*(x) = 2c/√λ (número de onda). Horizonte λ=1: k_* = 1.587.
- Λ_FP(x) = 2c²/√λ = c·k_* (frecuencia FP). Horizonte λ=1: Λ_FP = 1.260.
- ξ_λ = √λ/(√2 c), ξ_ref = 1/(√2 c), ℓ = √λ/(2c); ξ_λ/ℓ=√2, ξ_ref/ℓ=√(2/λ).
  Validado a máquina, bordes + horizonte.

**Convención / pendiente (NO fijar sin justificar):**
- **c₀ de referencia de FP:** su fórmula usa la velocidad de sonido de referencia del
  PERFIL FP, no automáticamente la local → elegir y justificar.
- **d≈1.58 (Ec.24):** FIT a los perfiles de FP → NO universal; transferible SOLO el
  escalado **d_ξ ~ c₀/(κ_K Λ_FP²)^{1/3}**, con Λ_FP de FRECUENCIA (no k_*).
- **Reparto d_ξ^L vs d_ξ^R y centro x̄_hor (Ec.20):** dependen del fondo, pendientes.

## 4. Siguiente paso (si se decide gastar presupuesto)
PASO 2: implementar κ̄ (Ec.20) = gradiente de (c−v) promediado sobre una ventana d_ξ
DERIVADA de nuestro fondo (sin importar d=1.58), con Λ_FP de frecuencia; medir
T_IR(w) con la maquinaria N_ω validada (H4); comparar **T_IR vs κ_local/2π vs κ̄/2π**.
Predicción FP falsable: T_IR sigue a κ̄/2π (horizonte ensanchado) mejor que a
κ_local/2π, y la diferencia importa más en horizontes empinados (w pequeño).

Archivos: `FP_step1_scales.py`.
