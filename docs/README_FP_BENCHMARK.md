# README — Reproducibilidad del benchmark Finazzi–Parentani (para revisión externa)

> **Nota de la versión pública (2026-11-17).** Este documento se escribió con todos los archivos en una sola carpeta. En el repositorio público están reorganizados: solver en `solver/`, benchmark FP en `benchmarks/fp_v1/`, cadena de validación en `methods/`, correlaciones en `correlations/`, documentos en `docs/`; `out/`, `reference/` y `residual_anatomy/` conservan su ruta. Todos los scripts se ejecutan **desde la raíz del repositorio** (p. ej. `python3 benchmarks/fp_v1/FP_check_reproduction.py`). Los SHA-256 citados aquí corresponden a la estructura original; los de la versión pública están en `docs/provenance/RELEASE_MANIFEST.txt`.

Guía mínima para que un investigador independiente (física de BEC / Bogoliubov–de Gennes)
reproduzca y examine el benchmark **sin reconstruir el proceso de desarrollo**. El
documento científico es **`FP_BENCHMARK.md`**; este README solo explica cómo ejecutarlo
y qué mirar.

---

## 0. Qué se afirma (y qué no)
**Afirmación (acotada):** para la familia de perfiles estudiada y bajo las convenciones
especificadas, la simulación reproduce la relación de temperatura de horizonte ensanchado
de Finazzi–Parentani: **T_upstream/(κ_local/2π) = 0.957 → 0.985** en cinco anchuras de
horizonte (tiende a 1 al suavizar el horizonte), usando el observable que FP define (el
cuanto Hawking upstream n_ω = |β_ω|²). Con la κ̄ de FP (T₀ = κ̄/2π) y d = 1.58 la desviación
es ≤ 0.75 %, **condicional a d** (con d = 1 o 2 llega al 2–3 %; tabla en `FP_BENCHMARK.md` §5b).
Es **validación numérica**, NO nueva física, y **NO una reproducción íntegra del
artículo** (perfiles distintos; d≈1.58 es ajuste de FP, no universal).

Referencia: F. Finazzi & R. Parentani, *Phys. Rev. D* **83**, 084010 (2011),
arXiv:1012.1556.

---

## 1. Entorno (fijado)
- Python **3.9.6** (`.python-version`); dependencias con versión fijada en
  **`requirements.txt`**: `numpy==2.0.2`, `scipy==1.13.1`, `sympy==1.14.0`.
- Instalación reproducible: `pip install -r requirements.txt`.
- Sin GPU, sin semillas aleatorias: **todo es determinista**. En la misma máquina
  reejecutar da los mismos números bit a bit. **Entre máquinas**, la aritmética float
  puede variar ligeramente según el BLAS/LAPACK de scipy; por eso la comparación con la
  referencia es por **tolerancia**, no por igualdad textual (§5).
- Todos los comandos se ejecutan desde el directorio del proyecto (`bec-gpe-solver/`).

---

## 2. Reproducción del resultado principal (1 comando, ~30 s)
```
python3 benchmarks/fp_v1/FP_benchmark_table.py
```
Produce la **Tabla §5** de `FP_BENCHMARK.md` y escribe `out/fp_benchmark.csv` (valores sin
redondear). Salida esperada (columnas clave; comparar con el CSV para plena precisión):

```
 w   kappa_local  kloc/2pi  kbar/2pi   T_up      T_tot   T_up/kloc  T_up/kbar  T_tot/kloc
 2     0.270283   0.043017  0.040972  0.041178  0.043417  0.9572     1.0050     1.0093
 3     0.180189   0.028678  0.027866  0.027904  0.029422  0.9730     1.0014     1.0259
 4     0.135142   0.021508  0.021089  0.021055  0.022199  0.9789     0.9984     1.0321
 6     0.090094   0.014339  0.014175  0.014099  0.014865  0.9833     0.9947     1.0367
 8     0.067571   0.010754  0.010670  0.010591  0.011166  0.9848     0.9926     1.0383
```
**Lo que hay que mirar:** la columna `T_up/kbar` = T_upstream/(κ̄/2π) ∈ [0.9926, 1.0050]
→ acuerdo ≤ 0.75 % con la predicción de FP **para d = 1.58** (condicional a d; la columna
`T_up/kloc` = 0.957–0.985 es el resultado que no depende de d). `T_up` es el observable de FP (canal Hawking
upstream); `T_tot` es un resultado SECUNDARIO (temperatura de la ocupación del compañero),
que **no** es la temperatura Hawking de FP.

**Verificación automática (recomendada):** `python3 benchmarks/fp_v1/FP_check_reproduction.py` regenera la
tabla y la compara con la referencia congelada por tolerancia (§5), imprimiendo la
desviación máxima y PASS/FAIL. Es la forma de confirmar la reproducción sin cotejar
números a mano.

NOTA: al importar `mode_matching` se ejecuta su autotest (imprime un banner de validación
del solver: homogéneo T=1/R=0, subsónico sin canal de norma −, horizonte β²→0.01008
convergente en N y L). Es esperado; no forma parte de la tabla.

---

## 3. Definición exacta del observable y del estimador (lo esencial a auditar)
- **Perfil (familia tanh, 1 parámetro w):** n(x)=n_u+(n_d−n_u)·½(1+tanh(x/w)), v=j/n,
  c=√(gn); ℏ=m=g=1; Mach upstream 0.5, downstream 2; horizonte c=v (raíz interpolada).
- **Estimador (def. FP, Ecs. 11–12):** T_ω = ω/ln(1+1/n), **T₀ = lím_{ω→0} T_ω**, por
  extrapolación cuadrática de T_ω sobre la banda ω/T₀∈[0.05,0.27] (9 puntos).
- **Observable de FP (Ec. 42):** n_ω = |β_ω|² = ocupación del canal u **upstream**. En el
  código es `n_L` = suma de |amp|² de los canales salientes de norma **+** en el lado
  **upstream (L)**, al inyectar el único modo entrante de **norma −** (ver §4). El
  co-propagante downstream n_{R+} = n^v_ω de FP (≪ n_ω) **no** se suma a la temperatura
  Hawking; T_total sí lo suma (= n_{−ω} de FP, el compañero).
- **Gravedad superficial local:** κ_local = |d(c−v)/dx|_H (derivada covariante exacta,
  verificable con `python3 solver/kappa_derive.py`).
- **κ̄ (FP Ec. 20), ventana simétrica d=1.58:** promedio del gradiente de (c−v) sobre
  d_ξ = d·c_H/(κ_local·Λ_FP²)^{1/3}, Λ_FP = 2c_H² (escala de frecuencia).

---

## 4. Controles independientes (cada uno 1 comando, ~30–40 s)
| Comando | Qué verifica | Esperado |
|---|---|---|
| `python3 FP_step3_audit.py` | **Un único canal IN de norma −1** (inyección completa); normalización Bogoliubov **\|α\|²−n_tot=1**; **termalidad** Planckiana de n_L y n_tot | IN: izq[+1], der[+1,−1]; \|α\|²−n_tot=1 a ~10⁻³; intercept ~4·10⁻⁴ |
| `python3 FP_step4_upstream_robust.py` | **Robustez de T_upstream** al variar banda y resolución N | T₀(n_L)/κ_local idéntico a 4 cifras en bandas [0.05,0.25]→[0.12,0.60]T; estable en N |
| `python3 FP_step2_audit.py` | Techo del gradiente: **cota del valor medio (demostrada)** + **max\|f′\|=1.014343·κ_local (verificado numéricamente)**; monotonía κ̄(d); presupuesto de error de T_total | max\|f′\| invariante en w; np.gradient→O(dx²) |
| `python3 FP_step1_scales.py` | Dispersión del solver vs analítica (aserciones) | max err 1.8·10⁻¹⁴; todas las aserciones PASAN |
| `python3 FP_step3_origin.py` | T_eff(ω) por canal, bandas descendentes (origen del exceso de T_total) | estable en banda; `out/fp_step3_origin.csv` |

Interpretación de los controles en `FP_BENCHMARK.md` §7 (con tolerancias y dominios).

---

## 5. Datos de referencia y comparación por tolerancia
Hay que distinguir **regenerar** la tabla de **compararla** con una referencia:
- **`reference/fp_benchmark_reference.csv`** — tabla CONGELADA de referencia (no
  modificar). Es el patrón contra el que comparar.
- **`out/fp_benchmark.csv`** — tabla que genera cada ejecución de `FP_benchmark_table.py`.
- **`python3 benchmarks/fp_v1/FP_check_reproduction.py`** — regenera y **compara con la referencia por
  TOLERANCIA mixta** |b−a| ≤ ATOL + RTOL·|a| (RTOL=10⁻⁴, ATOL=10⁻¹²; el término ATOL la
  hace segura cuando |a|≈0), muy por debajo de la señal física de 0.5–0.75 %. Imprime la
  desviación máxima y PASS/FAIL (exit 0/1). **No se exige igualdad textual** entre
  máquinas. El comprobador además: valida la integridad de filas (sin w ausentes,
  adicionales ni duplicados), **rechaza NaN/inf**, compara las 9 columnas de RESULTADO
  bajo tolerancia (las de diagnóstico —pseudo-u, error de extrapolación— solo se exige
  que sean finitas, pueden variar más entre plataformas), y **abre la referencia en solo
  lectura** (no la modifica). Verificado con pruebas negativas (valor perturbado, fila
  ausente, NaN → FAIL correcto).

**Qué significa (y qué no) un PASS:** que las 9 columnas de RESULTADO se han reproducido
dentro de las tolerancias definidas **y** que los diagnósticos son finitos. NO significa
que todos los detalles numéricos coincidan entre plataformas. El umbral **10⁻⁴ es un
UMBRAL DE REGRESIÓN del comprobador** (para detectar cambios de código/entorno), **NO una
estimación del error sistemático físico** de la simulación; la incertidumbre física se
discute aparte en `FP_BENCHMARK.md` §5a/§7 (y recuérdese que esos controles tampoco
constituyen una cota matemática del error sistemático del límite).

Columnas del CSV: `kappa_local, cH` (gravedad, c en horizonte); `kloc_over_2pi,
kbar_over_2pi` (referencias); `T_upstream, T_total` (temperaturas); `Tup_over_kloc2pi,
Tup_over_kbar2pi, Ttot_over_kloc2pi` (cocientes); `pseudo_u_max` (unitariedad);
`fit_rms_resid, order_dispersion, T_up_linear, T_up_cubic` (error de la extrapolación
ω→0). Datos de los controles: `out/fp_step2_kappabar.csv`, `out/fp_step3_origin.csv`.

Sobre el error de extrapolación: residuo RMS ~10⁻⁹ y sensibilidad al orden del polinomio
<0.005 %. Son indicios de incertidumbre numérica muy pequeña frente al 0.75 %, **pero NO
una cota matemática del error sistemático del límite** (los tres ajustes podrían compartir
sesgo). Ver `FP_BENCHMARK.md` §5a.

---

## 6. Mapa de archivos
| Archivo | Papel |
|---|---|
| **`FP_BENCHMARK.md`** | **Documento científico (leer primero).** Modelo, convenciones, Ecs. FP, correspondencia de canales, resultado, controles, limitaciones. |
| `FP_SCALES.md` / `FP_STEP2.md` / `FP_STEP3.md` / `FP_STEP4.md` | Detalle y auditoría por paso (escalas; κ̄+techo; canales S; correspondencia con el artículo). |
| `FP_benchmark_table.py` | Tabla principal → CSV. |
| `FP_step{1,2,3,4}_*.py` | Controles (§4). |
| `mode_matching.py` | **Solver de scattering BdG (CONGELADO).** Matching modal exacto; autotest al importar. |
| `kappa_derive.py` | κ covariante = \|c′−v′\| (derivación simbólica). |
| `requirements.txt`, `.python-version` | Entorno fijado (§1). |
| `reference/fp_benchmark_reference.csv` | Tabla de referencia congelada (§5). |
| `FP_check_reproduction.py` | Regenera y compara con la referencia por tolerancia (§5). |
| `FP_REVIEW_MANIFEST.txt` | SHA-256 de los 28 artefactos entregados (versión congelada). Verificar integridad: `grep '^[0-9a-f]' FP_REVIEW_MANIFEST.txt \| shasum -a 256 -c`. |
| `bdg_smatrix.py` | Matriz S 3×3 completa con test S†ηS=η (sobre el perfil fijo `out/waterfall_state.npz`; referencia cruzada de la estructura de canales). |

---

## 7. Qué escrutar con especial atención (limitaciones, `FP_BENCHMARK.md` §9)
1. **Familia de perfiles ≠ la de FP** (ellos: fondo + perturbación localizada, ajuste de
   d_ξ^L≠d_ξ^R a series de espectros). Aquí: familia tanh de un parámetro w, ventana
   **simétrica** d=1.58. Un ajuste con la κ̄ asimétrica de FP (Ec. 20-21) sería el cierre
   cuantitativo fino (no hecho).
2. **d≈1.58 no es universal** (ajuste de FP a sus perfiles).
3. **Ecuaciones de FP leídas del HTML de arXiv** (v2) por un modelo + verificación;
   **cotejar Ecs. 11, 12, 20, 39–42 (signos y convenciones) con el PDF** antes de
   cualquier publicación.
4. **Factor greybody no factorizado** (no necesario: n_L sale directo de la inyección
   única del modo −, que es columna completa de la matriz S).
5. Régimen **débilmente dispersivo** κ_local/Λ_FP ∈ [0.054, 0.215], caso simétrico
   monótono sin undulaciones.

---

## 8. Procedencia
Registro general del proyecto análogo-Hawking y matriz de afirmaciones: `REPRODUCE.md`,
`PROVENANCE.txt`. Este benchmark FP se construyó sobre ese solver ya validado (concordancia
transfer vs mode-matching: β a 10⁻⁴, N_ω a 10⁻⁵). Solver **congelado** a fecha 2026-11-15.
