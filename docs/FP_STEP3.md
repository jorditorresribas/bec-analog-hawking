# Experimento B — PASO 3: origen de T_IR > κ_local/2π (REVISA la lectura de Paso 2)

Sesión 2026-11-15. Programa de ChatGPT (pruebas 1–3), **sin tocar el solver BdG**.
Objetivo: antes de postular un mecanismo, distinguir si el exceso T_IR>κ_local está
(i) en la creación de pares del canal Hawking, (ii) en la suma/mezcla de canales que
hace el estimador H4 (greybody), o (iii) en un sesgo del límite IR (banda).
Archivo: `FP_step3_origin.py`. CSV: `out/fp_step3_origin.csv`.

## 0. Definición EXACTA del estimador H4 / N_ω (lo que realmente se mide)
`mode_matching.solve_inject` resuelve el scattering BdG inyectando el modo entrante de
**norma negativa** desde la derecha (downstream), amplitud 1, y cero desde la izquierda
(`incL=0`). Devuelve, por canal propagante SALIENTE, (signo de norma, lado, |amplitud|²).
Para el perfil con horizonte hay **tres canales salientes** (sonda verificada):

| canal | norma | lado | papel | |amp|² (w=8) |
|---|---|---|---|---|
| (+1, L) | +1 | izq (subsónico) | **cuanto Hawking** (upstream) | 2.814 |
| (−1, R) | −1 | der (supersónico) | **compañero** (partner) | 3.990 |
| (+1, R) | +1 | der (supersónico) | co-propagante downstream | 0.178 |

**N_ω(H4) = Σ |amp|² sobre los de NORMA +1 = (+1,L) + (+1,R).** La pseudo-unitariedad
impone Σ η|amp|² = −1 (norma conservada: +2.81 −3.99 +0.18 = −1.00). Luego
**T_IR = lím_{ω→0} ω·N_ω**. El estimador SUMA dos canales de norma +1: el Hawking
grande (L) y el co-propagante pequeño (R, ~5.5% de N_ω).

## 1. Auditoría S (pasos 1–2 de ChatGPT) — estructura, NO un veredicto
**Advertencia de honestidad:** una versión previa de este documento concluyó que el
exceso era "mero artefacto greybody" y que "apoya a FP". Eso fue una SOBRE-CORRECCIÓN.
La auditoría de la matriz S muestra una estructura más matizada; ningún lado está
cerrado.

**(1) Un ÚNICO canal IN de norma −1** (certificado w=2,4,8, ambos bordes): IN = izq
[+1], der [+1, −1]. Luego la inyección única del modo −1 da la **ocupación espontánea
COMPLETA** de cada canal saliente en el vacío entrante: n_i = |S_{i,p}|² (p = único IN
de norma −1). Resuelto el aviso §3 de ChatGPT ("si hay más de uno, no basta").

**(2) Normalización de Bogoliubov |α|² − n_tot = 1** (a ~10⁻³ en todo ω; α = canal
saliente de norma −1). Es decir **n_tot = |β_total|²** estándar, y e^{ω/T}=|α|²/n_tot=
(1+n_tot)/n_tot. Por tanto **T_eff,tot (= el estimador H4) ES la temperatura de balance
detallado del sector de norma +1, SIN ambigüedad de greybody.** (Esto contradice la
lectura previa de que T_tot era un artefacto.)

**(3) Ambos canales son Planckianos** con temperaturas DISTINTAS (test ln(1+1/n)=ω/T
lineal por el origen: intercept ~4·10⁻⁴, rms ~2·10⁻⁴):

| w | T(n_L)/κ_local (upstream, flujo observable) | T(n_tot)/κ_local (=H4, Bogoliubov total) | κ̄_FP/κ_local |
|---|---|---|---|
| 2 | 0.961 | 1.014 | 0.9525 |
| 4 | 0.983 | 1.037 | 0.9805 |
| 8 | 0.988 | 1.043 | 0.9922 |

La diferencia entre ambas es el canal co-propagante downstream (+1,R), él mismo
Planckiano (~5% de n_tot). **Las dos definiciones naturales de "la temperatura"
CABALGAN κ_local:** la ocupación total de Bogoliubov (n_tot) está por ENCIMA y crece
con w; el flujo Hawking upstream observable (n_L) está por DEBAJO, rastrea κ̄ y encoge
con w. Cuál corresponde al observable de FP es lo que queda por decidir (§4).

## 2. Prueba 3 (banda): el exceso del estimador NO es sesgo del límite IR
T_IR_tot/κ_local es estable a 4 cifras en tres bandas descendentes:

| w | [0.12,0.60]T | [0.08,0.40]T | [0.05,0.25]T |
|---|---|---|---|
| 2 | 1.0092 | 1.0093 | 1.0093 |
| 4 | 1.0320 | 1.0321 | 1.0321 |
| 8 | 1.0383 | 1.0383 | 1.0383 |

Convergencia en N (N→1.5N) estable a 5 cifras en el ω más bajo; pseudo-u mejora con N.
→ **Hipótesis C (sesgo de banda) RECHAZADA.** El valor del estimador sumado es real; lo
que cambia la física es QUÉ canal se identifica con la ocupación.

## 3. Prueba 2 (espectral): T_eff(ω) es casi meseta, no deriva
T_eff_tot/κ_local y T_eff_L/κ_local varían <0.5% en toda la banda ω/κ_local∈[0.008,0.096]
→ temperatura efectiva genuina por canal (meseta IR), no una desviación espectral que el
ajuste térmico comprima en un número.

## 4. CONSECUENCIA para la lectura de Paso 2 (`FP_STEP2.md`) — matizada
La **matemática** de Paso 2 (techo certificado max|f′|=1.014343·κ_local, cota del valor
medio, independencia de ventana para w≥4) **sigue siendo correcta**. Lo que cambia es la
INTERPRETACIÓN del observable que acota:
- T_IR_tot **sí es** una temperatura de balance detallado legítima (n_tot=|β_total|²,
  §1.2), NO un artefacto numérico. Que supere el techo del gradiente para w≥4 significa
  que **la ocupación TOTAL de Bogoliubov no la fija ningún promedio espacial de
  d(c−v)/dx**. Eso es un resultado real sobre esa definición de temperatura.
- PERO existe una SEGUNDA temperatura legítima, la del flujo Hawking **upstream
  observable** (n_L), que está por DEBAJO de κ_local y rastrea κ̄. Para ESA definición
  no hay tensión con FP.
- Por tanto NO se sostiene "FP falsado" (Paso 2) NI "FP confirmado / exceso artefacto"
  (versión previa de Paso 3). **El número único "temperatura" es ambiguo por la
  estructura multicanal; las dos definiciones naturales cabalgan κ_local.**

## 5. Qué queda ABIERTO (lo esencial, sin cerrar)
> **RESUELTO en PASO 4 (`FP_STEP4.md`):** la lectura del artículo (Ec. 42: n_ω=|β_ω|²
> del canal u upstream; T_0 de n_ω SOLO; comparada con κ̄/2π) — verificada
> independientemente por ChatGPT contra las Ecs. 39–42 — identifica el observable de FP
> con **n_L (upstream)**, no n_tot. n_tot = n_L+n_R+ = FP's n_{−ω} (compañero). Con el
> observable correcto, T_upstream≈κ̄/2π a <0.75% (robusto en banda y N) → la simulación
> REPRODUCE a FP. El punto 1 de abajo queda así respondido; 2–3 también (el canal (+1,R)
> es el v de FP, no parte de su T; la paradoja era de n_tot=compañero).

1. **¿Qué observable usa FP?** Hay que mapear la "temperatura de la radiación Hawking"
   de FP (su definición de espectro) a una de nuestras dos: ocupación total de
   Bogoliubov (n_tot) vs flujo upstream (n_L). Hasta hacerlo, la comparación con κ̄ no
   está cerrada. (El rastreo n_L≈κ̄ a sub-% es sugerente pero podría ser coincidencia.)
2. **El canal co-propagante (+1,R):** ¿es parte del "espectro Hawking" o una emisión
   downstream separada? Determina cuál de las dos T es la pertinente.
3. **La paradoja original persiste para n_tot:** T(n_tot)/κ_local CRECE hacia el límite
   hidrodinámico (w grande). Si n_tot es el observable correcto, esto sigue sin
   explicación dispersiva (§ original de Paso 2). Si lo es n_L, la paradoja se disuelve.
4. **Siguiente paso técnico:** reconstruir la matriz S COMPLETA (3×3) de la familia tanh
   con amplitudes complejas (no solo |amp|²) para verificar S†ηS=η canal a canal
   —`bdg_smatrix.py` lo hace pero sobre el perfil waterfall fijo, no la familia w— y, con
   ella, intentar la factorización greybody de cada canal. Después, si procede, barrido
   (w,λ) (prueba 4) para separar dispersión de suavidad.

Archivos: `FP_step3_origin.py` (pruebas 1-3), `FP_step3_audit.py` (auditoría S: canales
IN, Bogoliubov, termalidad), `out/fp_step3_origin.csv`. Base: `mode_matching.py`
(solve_inject, estructura de canales), `bdg_smatrix.py` (S completa, perfil waterfall),
`FP_STEP2.md` (κ̄, techo), `FP_SCALES.md` (escalas).
