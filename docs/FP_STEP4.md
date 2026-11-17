# Experimento B — PASO 4: correspondencia con Finazzi-Parentani (CIERRE conceptual)

Sesión 2026-11-15. Paso puramente documental (sin tocar el solver), el que pedía
ChatGPT: **determinar qué cantidad llama FP "temperatura de Hawking"** y mapearla a
nuestros observables. Fuente: arXiv:1012.1556v2 (HTML), Sec. III y Apéndices A–B.
Lectura verificada con DOS consultas independientes al HTML (citas verbatim abajo);
**salvedad: los números de ecuación y el indexado (esp. Ec. 40/42) convienen
cotejarse con el PDF antes de citarlos en un paper.**

## 1. Qué llama FP "temperatura de Hawking" (citas verbatim)
- **Espectro Hawking (Ap. B, Ec. 42):** `n_ω = |β_ω|²,  n^v_ω = |B̃_ω|²,
  n_{−ω} = n_ω + n^v_ω`. FP identifica **n_ω = |β_ω|²** como "the mean occupation
  number of phonons spontaneously emitted **to the right region, i.e. against the
  flow**" (su flujo va derecha→izquierda ⇒ es el **upstream**). Es un **canal ÚNICO**.
- **El canal v es aparte y pequeño:** "In the presence of dispersion, one generally
  finds **n^v_ω ≪ n_ω**"; n^v_ω = "the number of left moving phonons" (con el flujo,
  downstream). **NO se suma a la temperatura.**
- **Temperatura (Sec. III.1, Ec. 11–12):** `n_ω ≡ 1/(exp(ω/T_ω) − 1)`,
  `T_0 ≡ lím_{ω→0} T_ω` — definida **solo a partir de n_ω** (el canal u). Equivale a
  la pendiente de ln(1+1/n_ω) vs ω (mi estimador).
- **Comparación (Sec. III.3, Fig. 6):** T_0 se ajusta a **κ̄/2π** (su Ec. 20, la que
  implementamos), NO a κ_local: "The perfect agreement establishes that **the
  temperature T_0 equals κ̄/2π**." Su κ_K = nuestra κ_local.
- **Signo (Sec. III, Fig. 4):** en el caso **simétrico monótono** (κ_p>0), T_0 "starts
  from κ_B/2π for D_p→0, asymptotes to T_H=κ_K/2π for D_p larger than critical", con
  κ_B<κ_K ⇒ **T_0 < κ_K/2π = κ_local/2π** (horizonte ensanchado REDUCE la T). Salvedad
  de FP: con undulaciones T_0 puede SUPERAR κ_K/2π — pero nuestro perfil es el simétrico
  monótono, sin undulaciones (compatibilidad espectral, Paso 3 §3).

## 2. Tabla de correspondencia (concepto FP ↔ nuestro simulador)
Nuestro flujo es v>0 (izquierda→derecha), horizonte c=v; upstream = IZQUIERDA (subsónico),
downstream = DERECHA (supersónico). FP tiene el flujo invertido; "against the flow /
right region" de FP = nuestro upstream/izquierda.

| Concepto en FP | Ec. FP | Nuestra cantidad | Verificado |
|---|---|---|---|
| Cuanto Hawking u (contra el flujo, upstream) | n_ω=|β_ω|² (42) | **canal (+1,L) = n_L** | sí (norma +, lado upstream) |
| Fonón v (con el flujo, downstream), pequeño | n^v_ω=|B̃_ω|² (42) | **canal (+1,R) = n_R+** | sí; n_R+/n_L≈0.06 ≪1 (como FP) |
| Compañero (partícula de norma −) | n_{−ω}=n_ω+n^v_ω (42) | **n_tot = n_L+n_R+** | sí; |α|²−n_tot=1 (Bogoliubov) |
| Temperatura Hawking T_0 | 1/(e^{ω/T}−1) (11-12) | **T_upstream = T(n_L)** | pendiente ln(1+1/n_L) |
| Gravedad superficial promediada | κ̄ (20) | κ̄ (FP_STEP2) | sí |
| Gravedad superficial local κ_K | κ_K (2,8) | κ_local = |c′−v′| | sí |

## 3. RESOLUCIÓN de la ambigüedad del Paso 3
El Paso 3 dejó dos temperaturas Planckianas que cabalgan κ_local (T_upstream=T(n_L)<κ_local
y T_total=T(n_tot)>κ_local) y preguntaba cuál es la de FP. **El artículo lo decide sin
ambigüedad:** la temperatura Hawking de FP es **T_0 de n_ω=|β_ω|² = el canal u upstream
SOLO = nuestro n_L = T_upstream.** El estimador H4 sumaba n_L+n_R+ = **n_{−ω}
(ocupación del COMPAÑERO)**, que FP NO llama temperatura Hawking.

**Consecuencia — nuestra simulación REPRODUCE a FP:**

| w | T_upstream/κ_local (= FP T_0/T_H) | κ̄/κ_local (FP Ec.20) | T_upstream/κ̄ |
|---|---|---|---|
| 2 | 0.961 | 0.953 | 1.008 |
| 4 | 0.983 | 0.981 | 1.002 |
| 8 | 0.988 | 0.992 | 0.996 |

- T_upstream está **por debajo de κ_local** y **= κ̄/2π a <1%** — exactamente el resultado
  central de FP (T_0 = κ̄/2π < κ_K/2π en el caso simétrico monótono).
- El déficit **encoge con w** (T_upstream/κ_local: 0.961→0.988): T_0 tiende a κ_K/2π al
  suavizar, como FP (T_0→T_H para D_p grande). Dirección física correcta.

## 4. Qué queda de los pasos 2 y 3 (honesto)
- **Paso 2 ("falsación para w≥4"):** acotaba T_total = T(n_tot) = temperatura de la
  ocupación del COMPAÑERO n_{−ω}, que **NO es el observable Hawking de FP**. Por tanto
  esa "falsación" **no aplica a FP**. (La matemática del techo sigue siendo correcta como
  enunciado sobre n_tot, pero n_tot no es la T de Hawking.)
- **Paso 3 v1 ("apoya a FP"):** la INTUICIÓN (n_L es el canal correcto) era correcta; se
  sobre-afirmó sin base y luego se sobre-corrigió a "ambiguo". **El artículo cierra la
  ambigüedad a favor de n_L.** Con el observable de FP, la simulación es CONSISTENTE con
  FP (T_0=κ̄/2π). 
- **La "paradoja" del exceso creciente hacia el límite hidrodinámico** era un artefacto
  de mirar n_tot (compañero). En n_L (Hawking) el déficit encoge hacia hidrodinámico,
  sin paradoja.

## 4b. Robustez del observable upstream (ChatGPT §5.4). `FP_step4_upstream_robust.py`
T_0 = lím_{ω→0} T_ω del canal u (def. FP Ec. 11-12), variando banda y resolución:

| w | T_0/κ_local (bandas [0.05,0.25]→[0.12,0.60]T) | conv. N (N→1.5N) | **T_upstream/(κ̄/2π)** |
|---|---|---|---|
| 2 | 0.9572 (idéntica a 4 cifras) | 0.9576→0.9576 | **1.0050** |
| 4 | 0.9789 (idéntica a 4 cifras) | 0.9793→0.9792 | **0.9984** |
| 8 | 0.9848 (idéntica a 4 cifras) | 0.9851→0.9851 | **0.9926** |

El T_0 extrapolado a ω→0 (la def. de FP, Ec.12) es **independiente de la banda** a 4
cifras y estable en N. **T_upstream ≈ κ̄/2π a <0.75%** en las tres anchuras. El acuerdo
con FP es robusto, no un accidente de banda ni de resolución.

## 5. Caveats y posible cierre cuantitativo (opcional)
- Correspondencia basada en extracción del HTML por un modelo lector; **cotejar Ec. 11,
  12, 20, 42 con el PDF** antes de usarla en un manuscrito.
- κ̄ de FP usa reparto ASIMÉTRICO d_ξ^L≠d_ξ^R ajustado (Ec. 21); nuestra κ̄ usó ventana
  simétrica con d=1.58. El acuerdo T_upstream≈κ̄ a <1% ya es notable; un ajuste con la
  κ̄ asimétrica de FP (Ec. 20-21) sería el cierre cuantitativo fino, si se decide gastar.
- No se ha reconstruido la S compleja 3×3 de la familia tanh; no hace falta para esta
  correspondencia (n_L=|β_ω|² sale directo de la inyección única del modo de norma −1,
  que es columna completa, Paso 3 §1).

Correspondencia verificada INDEPENDIENTEMENTE por ChatGPT contra el artículo original
(Ecs. 39–42, Sec. III, Ap. B): identificación de canales confirmada como correcta, no
solo plausible. Observable principal de comparación con FP = **T_upstream/(κ̄/2π)**;
T_total = T(n_tot) se conserva como resultado secundario (temperatura de la ocupación
del compañero n_{−ω}), sin confundirlo con la T Hawking de FP.

Archivos: `FP_step4_upstream_robust.py` (robustez §4b); resto documental. Base:
`FP_STEP3.md` (auditoría S, canales, n_L/n_tot), `FP_STEP2.md` (κ̄, Ec.20), `FP_SCALES.md`
(escalas, Λ_FP). Fuente: arXiv:1012.1556v2 Sec. III + Ap. A–B.
