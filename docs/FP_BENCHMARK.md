# Benchmark Finazzi–Parentani en un agujero negro acústico BEC (1+1D)

> **Nota de la versión pública (2026-11-17).** Este documento se escribió con todos los archivos en una sola carpeta. En el repositorio público están reorganizados: solver en `solver/`, benchmark FP en `benchmarks/fp_v1/`, cadena de validación en `methods/`, correlaciones en `correlations/`, documentos en `docs/`; `out/`, `reference/` y `residual_anatomy/` conservan su ruta. Todos los scripts se ejecutan **desde la raíz del repositorio** (p. ej. `python3 benchmarks/fp_v1/FP_check_reproduction.py`). Los SHA-256 citados aquí corresponden a la estructura original; los de la versión pública están en `docs/provenance/RELEASE_MANIFEST.txt`.

### Documento autocontenido — modelo, convenciones, ecuaciones de referencia, controles numéricos, resultados y limitaciones

Sesión 2026-11-14/15. **Solver CONGELADO; benchmark CERRADO.** Objetivo: que un físico
independiente pueda reproducir y verificar el resultado sin reconstruir la conversación.
Resultado principal, en la formulación precisa (**no se afirma reproducción íntegra del
artículo**): **la simulación reproduce, PARA LA FAMILIA DE PERFILES ESTUDIADA Y BAJO LAS
CONVENCIONES ESPECIFICADAS, la relación de temperatura de horizonte ensanchado de
Finazzi–Parentani (T_0 = κ̄/2π) con una desviación inferior al 0.75 % en cinco anchuras,
usando el observable que FP define (el cuanto Hawking upstream n_ω=|β_ω|²), y el acuerdo
es robusto en banda de ajuste y resolución numérica.** Es validación numérica, no nueva
física (ver §8–9).

Documentos de apoyo (derivaciones y auditorías por paso): `FP_SCALES.md` (escalas),
`FP_STEP2.md` (κ̄ + techo del gradiente), `FP_STEP3.md` (auditoría de canales S),
`FP_STEP4.md` (correspondencia con el artículo). Este documento los resume y es la
referencia única.

---

## 1. Modelo y fondo
Condensado 1+1D, GPE con ℏ=m=g=1. Fluctuaciones = Bogoliubov–de Gennes (BdG) sobre un
flujo estacionario n(x), v(x) con corriente j=n·v constante. Perfil de horizonte
acústico (familia tanh, un parámetro de anchura w):
```
s(x) = ½(1+tanh(x/w));  n(x) = n_u + (n_d−n_u)·s(x);  v(x) = j/n(x);  c(x) = √(g·n)
```
Asintóticos FIJOS (independientes de w): upstream Mach_u=0.5 (subsónico), downstream
Mach_d=2 (supersónico). Constantes: n_u=1, v_u=0.5, j=0.5, n_d=(j/(Mach_d·√g))^{2/3}=0.39685,
μ=g·n_u+½v_u²=1.125. Horizonte acústico en c=v (raíz INTERPOLADA, |c−v|<5·10⁻⁷).
c_H=v_H=0.79370 en TODAS las w (sólo cambia la pendiente). Flujo v>0 (izq→der):
**upstream = izquierda (subsónico), downstream = derecha (supersónico).**

## 2. Convenciones y escalas (derivadas y validadas: `FP_SCALES.md`)
- **Dispersión BdG** (ℏ=m=g=1, knob λ, Bogoliubov exacto en λ=1):
  `(ω−vk)² = c²k² + (λ/4)k⁴`. Todo este benchmark usa **λ=1**.
- **Escala de frecuencia de FP** (su Ec. 9, Λ escala de FRECUENCIA): `Λ_FP = 2c²/√λ = c·k_*`
  (en λ=1, horizonte: Λ_FP=1.25992). k_* = 2c/√λ es un número de onda (≠ Λ_FP).
- **Gravedad superficial local** κ_local = |d(c−v)/dx|_H. DERIVADA simbólicamente de la
  métrica acústica (Unruh/Visser): κ covariante = |c′−v′|_H EXACTAMENTE (ratio=1, sin
  factor c_H; `kappa_derive.py`). Es la gravedad superficial, no una candidata.
- **κ̄ de FP (Ec. 20), horizonte ensanchado:** promedio uniforme del gradiente sobre una
  ventana d_ξ: `κ̄ = |[(c−v)(x_H+d_ξ^R) − (c−v)(x_H−d_ξ^L)]| / d_ξ`. Ventana (Ec. 24):
  `d_ξ = d·c_H/(κ_local·Λ_FP²)^{1/3}`, con **d≈1.58 AJUSTADO por FP a SUS perfiles (NO
  universal)**; sólo el escalado es transferible. Aquí usamos ventana SIMÉTRICA
  d_ξ^L=d_ξ^R=d_ξ/2 con d=1.58 (limitación §7).

## 3. Estructura de canales y estimador de ocupación espontánea (`FP_STEP3.md`)
Scattering BdG por **matching modal exacto** (`mode_matching.solve_inject`), operador FD
interior + forma modal asintótica impuesta en los bordes. Juez independiente del transfer
(validado: homogéneo T=1/R=0/β=0; subsónico sin canal de norma −; con horizonte β²
converge en N y L a 0.01008).

Canales (verificado, `FP_step3_audit.py`): **un único canal entrante de norma −1**
(IN: izq [+1], der [+1,−1]) → la inyección de ese modo (incL=0, amplitud 1) da la
**ocupación espontánea COMPLETA** de cada canal saliente en el vacío entrante,
n_i = |S_{i,p}|². Tres canales salientes:

| canal | norma | lado | papel | n_i (w=8) |
|---|---|---|---|---|
| (+1, L) | +1 | izq (upstream) | cuanto Hawking u | n_L = 2.814 |
| (+1, R) | +1 | der (downstream) | modo v co-propagante (pequeño) | n_R+ = 0.178 |
| (−1, R) | −1 | der (downstream) | compañero (negativo) | \|α\|² = 3.990 |

**Normalización de Bogoliubov** verificada: |α|² − n_tot = 1 (a ~10⁻³), con
n_tot=n_L+n_R+. **Termalidad:** n_L y n_tot son ambos Planckianos (ln(1+1/n)=ω/T lineal
por el origen, intercept ~4·10⁻⁴, rms ~2·10⁻⁴), con temperaturas DISTINTAS.

## 4. Correspondencia con Finazzi–Parentani (`FP_STEP4.md`; arXiv:1012.1556v2 Sec. III, Ap. B)
Verificada contra las Ecs. 39–42 (confirmada independientemente). FP tiene el flujo
invertido; "against the flow / right region" de FP = nuestro upstream.

| FP | Ec. FP | Nuestra cantidad |
|---|---|---|
| n_ω = \|β_ω\|² (cuanto Hawking u, upstream) | 42 | **n_L** |
| n^v_ω (modo v, con el flujo, ≪ n_ω) | 42 | n_R+ |
| n_{−ω} = n_ω + n^v_ω (compañero) | 42 | n_tot = n_L+n_R+ |
| T_0 = lím_{ω→0} T_ω, con n_ω ≡ 1/(e^{ω/T}−1) | 11,12 | **T_upstream = T_0(n_L)** |
| comparada con κ̄/2π ("T_0 equals κ̄/2π", Fig. 6) | 20 | κ̄/2π |

**El observable Hawking de FP es n_ω = n_L (canal u upstream SOLO), NO la suma.** FP
compara T_0 con κ̄/2π, y en el caso simétrico monótono (el nuestro) halla T_0 < κ_K/2π
(= κ_local/2π). El estimador "H4" previo sumaba n_L+n_R+ = n_{−ω} (compañero): observable
distinto (ver §6, resultado secundario).

## 5. RESULTADO PRINCIPAL (`FP_benchmark_table.py` → `out/fp_benchmark.csv`)

### 5a. Estimadores — definición ÚNICA (para no mezclar magnitudes)
Hay UN estimador de temperatura por canal, el de FP (Ec. 11-12):
`T_ω ≡ ω/ln(1+1/n)`,  `T_0 ≡ lím_{ω→0} T_ω`, obtenido por **extrapolación cuadrática**
de T_ω(ω) a ω=0 sobre la banda om/T∈[0.05,0.27] (malla fija, 9 puntos; caja
L=max(40,8w+20), N=240·L, horizonte interpolado). Dos observables, MISMO estimador:
- **T_upstream ≡ T_0(n_L)** — canal u (Hawking de FP). ES el observable de comparación.
- **T_total ≡ T_0(n_tot)** — ocupación del compañero n_{−ω} (resultado secundario §6).

Las columnas "T_up/κ_local·2π" y "T_up/(κ̄/2π)" son el **mismo T_upstream** dividido por
dos referencias distintas (κ_local/2π y κ̄/2π), no estimadores distintos. Control de
estimador (`FP_step4_upstream_robust.py`): T_0(n_L) por **extrapolación ω→0 es
independiente de la banda a 4 cifras** ([0.05,0.25]→[0.12,0.60]T); el ajuste de
PENDIENTE de ln(1+1/n_L) (asume Planck puro en la banda) da un valor ~0.1-0.4 % mayor y
dependiente de banda — por eso el benchmark usa la extrapolación ω→0 (def. exacta de FP,
Ec. 12), no la pendiente.

**Error de aproximar el límite ω→0** (guardado en el CSV, columnas `fit_rms_resid`,
`order_dispersion`): la sensibilidad de la extrapolación al orden del polinomio
(lineal/cuadrático/cúbico) es inferior al 0.005 %, y el residuo RMS del ajuste cuadrático
de T_ω(ω) es del orden de 10⁻⁹, en las cinco anchuras. Ambos controles indican que la
incertidumbre numérica asociada a la extrapolación es muy pequeña frente a la desviación
máxima del 0.75 %. **No constituyen, por sí solos, una cota matemática del error
sistemático del límite** (los tres ajustes podrían compartir parte del mismo sesgo).

### 5b. Tabla
Valores de alta precisión (CSV con todos los dígitos):

Columnas: κ_local y las /2π son GRAVEDADES (referencias); T_upstream es una TEMPERATURA;
las dos últimas son COCIENTES adimensionales del mismo T_upstream sobre dos referencias.

| w | κ_local | κ_local/2π | κ̄/2π (d=1.58) | **T_upstream** | **T_upstream/(κ̄/2π)** | T_upstream/(κ_local/2π) |
|---|---|---|---|---|---|---|
| 2 | 0.270283 | 0.043017 | 0.040972 | 0.041178 | **1.0050** | 0.9572 |
| 3 | 0.180189 | 0.028678 | 0.027866 | 0.027904 | **1.0014** | 0.9730 |
| 4 | 0.135142 | 0.021508 | 0.021089 | 0.021055 | **0.9984** | 0.9789 |
| 6 | 0.090094 | 0.014339 | 0.014175 | 0.014099 | **0.9947** | 0.9833 |
| 8 | 0.067571 | 0.010754 | 0.010670 | 0.010591 | **0.9926** | 0.9848 |

- **T_upstream ≈ κ̄/2π a ≤ 0.75 %** en las cinco anchuras (desv.: +0.50, +0.14, −0.16,
  −0.53, −0.74 %). Reproduce la relación T_0 = κ̄/2π de FP.
- T_upstream < κ_local/2π siempre, y el déficit ENCOGE con w (0.957→0.985): T_0 → κ_K/2π
  al suavizar el horizonte, como FP (T_0 → T_H para D_p grande).

## 6. Resultado secundario y por qué NO es la T de FP
T_total = T_0(n_tot) con n_tot = n_{−ω} (compañero): T_total/(κ_local/2π) = 1.009, 1.026,
1.032, 1.037, 1.038 (w=2,3,4,6,8) — por ENCIMA de κ_local y creciendo con w. Es una
temperatura de Bogoliubov legítima de la ocupación del compañero, **pero NO el observable
Hawking de FP.** (La "falsación para w≥4" explorada en `FP_STEP2.md`, basada en la cota
del valor medio κ̄ ≤ max|f′| —demostrada— con el valor max|f′|=1.014343·κ_local
—verificado numéricamente, §7—, es un enunciado correcto sobre T_total, pero T_total no es
la temperatura Hawking de FP; por eso no refuta a FP.)

## 7. Controles numéricos (preservados, no sólo cocientes)
- **Dispersión del solver vs analítica** (bordes + horizonte interpolado): max err
  1.8·10⁻¹⁴; símbolo del operador vs ω: ≤2.4·10⁻¹⁵ (`FP_step1_scales.py`).
- **Horizonte interpolado:** |c_H−v_H| < 5·10⁻⁷ (no punto de malla).
- **Pseudo-unitariedad** del scattering (|Σ η|a|² − (−1)|): ≤ 3·10⁻³ en toda la malla,
  mejora con N. Normalización Bogoliubov |α|²−n_tot=1 a ~10⁻³.
- **Convergencia L,N de T_upstream:** T_0(n_L)/κ_local estable a 4 cifras variando banda
  [0.05,0.25]→[0.12,0.60]T y N→1.5N (`FP_step4_upstream_robust.py`). Ej. w=8: 0.9848 en
  todas las bandas; N=20160→30240: 0.9851→0.9851.
- **Convergencia de T_total** (independiente): estable a 4 cifras con L=40→130, N≤30240
  (`FP_STEP2.md`/§convergencia).
- **Techo del gradiente** (`FP_step2_audit.py`), con tres niveles distinguidos:
  - **DEMOSTRADO analíticamente:** κ̄ = |Δf/d_ξ| ≤ max_{[x_L,x_R]}|f′| para CUALQUIER
    ventana (teorema del valor medio). La derivada f′=c′−v′ tiene forma CERRADA.
  - **VERIFICADO numéricamente:** el valor max_x|f′| = 1.014343·κ_local (invariante en w)
    se obtiene por búsqueda de raíces de f″=0 (forma cerrada) en un dominio finito
    [x_H−30w, x_H+30w] con |f′| en los bordes = 0 a máquina (decae e^{−2|x|/w}). Es una
    certificación numérica del máximo GLOBAL, NO una prueba analítica exhaustiva de que
    no haya otro máximo fuera del dominio escaneado (aunque la cola e^{−2|x|/w} lo
    excluye con holgura).
  - **Tolerancia/dominio:** `np.gradient` converge al valor cerrado O(dx²) (err <2·10⁻⁶
    con N=20000); dominio ±30w; bordes 0.0 a máquina.
- **Termalidad:** intercept ~4·10⁻⁴, rms ~2·10⁻⁴ en ln(1+1/n)=ω/T (Planckiano).

## 8. Qué se puede afirmar y qué NO
**Se afirma:** con estas convenciones y esta familia de perfiles, la simulación reproduce
la relación de temperatura de horizonte ensanchado de FP (T_0 = κ̄/2π) con desviación
≤ 0.75 %, usando el observable que FP define (n_ω = canal u upstream), con controles de
estabilidad en banda y resolución. Es un resultado de **validación numérica**, no de
nueva física.

**NO se afirma:** haber reproducido íntegramente el artículo, ni que la predicción valga
para cualquier perfil. En particular, la discrepancia previa (T_total > κ_local) NO era
una refutación de FP: comparaba la temperatura del compañero con la del cuanto Hawking.

## 9. Limitaciones (honestas)
1. **Familia de perfiles ≠ la de FP.** FP construye fondos (componente de fondo +
   perturbación localizada) y ajusta d_ξ^L, d_ξ^R a series de espectros; aquí: familia
   tanh de un parámetro w con la Ec. 24 aplicada directamente y ventana SIMÉTRICA d=1.58.
2. **d≈1.58 NO es universal** (ajuste de FP a sus perfiles). El acuerdo ≤0.75 % con la κ̄
   simétrica es notable; un ajuste con la κ̄ ASIMÉTRICA de FP (Ec. 20-21) sería el cierre
   cuantitativo fino (no hecho).
3. **Ecuaciones de FP leídas del HTML** (arXiv:1012.1556v2) por un modelo lector +
   verificación independiente; **cotejar Ecs. 11, 12, 20, 39–42 (signos y convenciones)
   con el PDF** antes de usar en un manuscrito.
4. **Factor greybody no factorizado** (no necesario: n_L sale directo de la inyección
   única del modo −1, columna completa). No se reconstruyó la S compleja 3×3 de la
   familia tanh (sí existe para el perfil waterfall fijo: `bdg_smatrix.py`).
5. Régimen débilmente dispersivo κ_local/Λ_FP ∈ [0.054, 0.215]; caso simétrico monótono
   sin undulaciones (donde FP admite T_0 > κ_K/2π).

## 10. Reproducibilidad
Entorno: Python 3 + numpy + scipy. ℏ=m=g=1. Todo determinista (sin ruido).
| Script | Produce |
|---|---|
| `FP_step1_scales.py` | escalas + validación dispersión (aserciones) |
| `FP_step2_kappabar.py` | κ̄(w), T_total(w), techo; `out/fp_step2_kappabar.csv` |
| `FP_step2_audit.py` | techo del gradiente (cota demostrada + max|f′| verificado), monotonía κ̄(d), presupuesto error T_total |
| `FP_step3_origin.py` | T_eff(ω) por canal, bandas; `out/fp_step3_origin.csv` |
| `FP_step3_audit.py` | auditoría S: canales IN, Bogoliubov, termalidad |
| `FP_step4_upstream_robust.py` | robustez banda/N de T_upstream |
| **`FP_benchmark_table.py`** | **tabla principal §5; `out/fp_benchmark.csv`** |
| `mode_matching.py` | solver de scattering (juez modal) — CONGELADO |
| `kappa_derive.py` | κ covariante = |c′−v′| (simbólico) |
Núcleo del resultado: ejecutar `FP_benchmark_table.py` reproduce la tabla §5.
**`out/fp_benchmark.csv`** guarda TODOS los valores sin redondear. Columnas:
`kappa_local, cH` (gravedad, velocidad del sonido en el horizonte);
`kloc_over_2pi, kbar_over_2pi` (gravedades /2π = referencias); `T_upstream, T_total`
(TEMPERATURAS, def. FP Ec.12); `Tup_over_kloc2pi, Tup_over_kbar2pi, Ttot_over_kloc2pi`
(COCIENTES adimensionales); `pseudo_u_max` (control de unitariedad); `fit_rms_resid,
order_dispersion, T_up_linear, T_up_cubic` (error de la extrapolación ω→0, §5a).
Referencia: F. Finazzi, R. Parentani, *Phys. Rev. D* 83, 084010 (2011), arXiv:1012.1556.
