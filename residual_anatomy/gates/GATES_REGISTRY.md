# Registro de gates — realización Padé disperso (Pasos 2–3)

> **★ VALIDACIÓN DEL OPERADOR CERRADA Y APROBADA DOCUMENTALMENTE por ChatGPT (2026-11-15)** dentro de su
> alcance documentado → `OPERATOR_VALIDATION_CLOSURE.md` (2 correcciones finales aplicadas: C1 "principal
> sospechosa, no causa exclusiva"; B1 simetría "para la realización ensayada, con limitaciones de borde"). Cubre P2, P3.1, P3.2-A/B1/B2/B2c/B3, Gates 1/2a y **P3.4-C2**
> (conservación dinámica a nivel de operador). **C1 (pseudo-unitaridad convergente de S) NO validado:
> RESERVADO** para el solver estabilizado → `STABILIZED_SOLVER_RESERVATION.md` (criterios a priori, sin
> arquitectura). Pregunta física `n_L(λ)` = fase posterior, depende del solver Padé.

Línea aislada sobre el solver FP v1.0 **congelado**. Ningún gate modifica v1.0 ni el benchmark.
Cada entrada: **hipótesis · entrada · resultado esperado · tolerancia · estado**.

Orden de trabajo (ChatGPT, 2026-11-15): migrar tests → registro → prueba A³ → pseudo-herm BdG completo
→ regresión λ=1 (Gate 2b). Solo tras A³ + pseudo-herm se mide física.

---

## PASO 2 — verificaciones analíticas (CERRADAS bajo hipótesis declaradas)
Archivo: `p2_analytic.py`. Documentación: `PADE_SPEC.md §9`.

| # | Gate | Hipótesis | Entrada | Esperado | Tol | Estado |
|---|------|-----------|---------|----------|-----|--------|
| P2.1 | Recuento BC (Jordan) | cadena orden 6, bloque Jordan ±κ tam. 3 | sistema 1er orden 6×6 | subesp. estable dim 3; rank B=3; B·estable=0 | <1e-12 | ✅ |
| P2.2 | Positividad de L | 2G+D Robin(+κ) pos. def. | perfil homog G=1 y tanh real | min autoval >0 (homog 2.003; tanh 0.800) | >0 | ✅ |
| P2.3 | Simetría K_P (inner cancela) | A,A³ autoadj → DAD,D²A³D² simétricos | K_P discreto homog, Dirichlet | asim. rel. ~máquina (salvo cond. D² densa) | <1e-7 rel | ✅ |
| P2.4 | Signo Green §9.3 | J1 = +𝒲(w,p_u)+𝒲(p_w,u) | campos reales sin BC | RHS₊ ≡ LHS; RHS₋ descartado | err FD | ✅ |
| P2.5 | Residuo 𝒲₂ (sympy) | BC (∂±κ)u=0,(∂±κ)u''=0 anulan 𝒲₂ si G'=0 | símbolos borde | 𝒲₂=0 (G'=0); residuo 2G'ūq₃ (G'≠0) | exacto | ✅ |
| P2.6 | Convergencia G' (tanh) | G'(x_b)~e^{−L/w} → residuo→0 | tanh w=2, L∈{10..50} | res_gauss,res_plane ↓ con L | ↓ monótono | ✅ |
| P2.7 | Flujo δ_J | ε_P'(k)≈ε_λ'(k) a orden Padé | λ=2,c²=G=1,k≤k_*=1.6 | max δ_J=0.655%<1% (en k_*) | <1% | ✅ |

**Pendiente teórico (NO cerrado):** maximalidad del dominio (autoadjunción, no solo simetría);
autoadjunción del operador INHOMOGÉNEO (no es multiplicador de Fourier). Ver §9.6.

---

## PASO 3 — gates numéricos (deciden la validez para física)
NO empezar la física hasta pasar A³ + pseudo-herm del BdG completo.

| # | Gate | Hipótesis | Entrada | Esperado | Tol | Estado |
|---|------|-----------|---------|----------|-----|--------|
| P3.1 | **A³ vs Green** (`gate_A3_kernel.py`) | band-solve A³ = recta completa (símbolo 1/(k²+κ²)³) en interior | G=1 homog; gauss + onda plana; BC∈{Dir,Robin}; L∈{10,20,40}; barrido dx | interior O(dx²) (pendiente 2.00, ratio 4.00); borde↓e^{−κL} (loc.); Robin<Dir a L fijo | interior <1e-4 | ✅ |
| P3.2-A | Equivalencia Schur (aumentado) (`gate_P3_2A_schur.py`) | matriz aumentada dispersa (4N) + solve reproduce K_P | homog+tanh, Dir+Robin; barridos N y L; 5 controles ChatGPT | res por bloque p,q1,q2,q3 ~1e−15; aug≡seq≡denso (rel ~1e−12) ∀N,L,BC; cond+residuos A³ | rel <1e−10 | ✅ (código revisado por ChatGPT: componente algebraico OK; barrido L añadido) |
| P3.2-B1 | Simetría por operador/discretización (Control 0) (`gate_P3_2B1_symmetry.py`) | D,A,K_P,C euclídeo+W, abs+rel, por discretización | plana/robin/dir-int; conv N,dx | plana K_P eucl-s.a.; robin A W-s.a. no K_P; C O(dx) (v1.0) | diagnóstico | ✅ (código revisado por ChatGPT) |
| P3.2-B2 | Identidad de Green discreta (`gate_P3_2B2_green.py`) | B2.1 simetría matricial; B2.2 flujo exacto conmutador F=dx·wᵀ[P_R,K_P]u (+ SBP-H trapezoidal observado); B2.3 convergencia al concomitante continuo (+ no-locales) | Dirichlet-fantasma, G=1, lorentzianas asim. | B2.1 K_P†=K_P ~máq ✅; B2.2 F exacto def.+SBP ~8e−15 observado; B2.3 F→BF_disc pendiente=1.07 (O(dx)), dif_naive→4.33e−5≠0 | — | ◐ B2.1 aprobado (ChatGPT); B2.2 OBSERVADO (no identidad SBP formal, alcance limitado por decisión); B2.3 límite continuo CERRADO por B2c |
| P3.2-B2c | **Cierre del límite continuo de B2.3** (`gate_P3_2B2c_continuum.py`) | referencia continua INDEPENDIENTE (convolución Green de línea g₁,g₂,g₃, sin malla FD) → F_matriz y BF_disc convergen al concomitante CONTINUO, errores por separado | G=1, Dirichlet-fantasma, región [−6,6], ±6 nodos exactos | BF_cont=−3.3390424e−3 (7 díg. estable ∀R,h); **BF_disc→BF_cont O(dx²) pend.2.00**; F_matriz→BF_cont O(dx) pend.0.94 (proyector sharp); F_matriz→BF_disc pend.1.01 | pend.>0 def. | ✅ CERRADO (aprobado ChatGPT 2026-11-15; h→O(h²), R→máquina) |
| P3.2-B3 | η-pseudo-herm de M_eff | M_eff con 𝒢=diag(W,−W); Dirichlet-fantasma; pre-req K_Pᵀ=K_P | R=M†𝒢−𝒢M (`gate_P3_2B3_pseudohermiticity.py`); A3-matriz (max/espectral) y A3-suave (QR ortonorm., paquetes baja k); control homog vs waterfall | A1 fuera-diag=0 EXACTO ✓; A2 PKP/Pdia_bulk_rel ~máq, KPsym_rel crece (cond. denso) ✓; **A3-matriz=O(1) (artefacto rejilla); A3-suave=O(dx²) pend 2.00 (QR)**; control homog suave ~máquina; A4a local | A1,A2(G+Qp) ✓; A2(K_P)=limitado por cond. densо PERO filtrado en proyección (proj_KP~máquina); A3-matriz=O(1) caracterizado; A3-suave O(dx²) ✓; A4a ✓ | ◐ EJECUTADO+QR+descomp.proyectada 2026-11-15. **Descomposición proyectada U†RU: proj_C=O(dx²) pend 2.00; proj_KP~máquina (1e−12, la contaminación de K_P NO alcanza el subespacio suave); proj_dia=0; recon_err=0**. El O(dx²) físico procede ENTERAMENTE de C. Defecto matricial O(1) persiste en modos de rejilla (no demostrada pseudo-herm matricial exacta). Pre-producción N grande → realización en BANDA (nueva realización, controles propios). **✅ CERRADO/APROBADO por ChatGPT 2026-11-15 con ALCANCE LIMITADO** (diagnóstico de pseudo-herm PROYECTADA sobre subespacio suave para la realización Dirichlet-fantasma; NO demostración de pseudo-herm matricial exacta ni autoadjunción continua). Siguiente: P3.4. |
| P3.3 | Regresión λ=1 (Gate 2b) | en λ=1 K_P=D/2 local → benchmark | solver scattering λ=1 | n_L(ω), T_upstream ≡ FP v1.0 (β²=0.01008, T/(κ̄/2π)=1.005) | <1% | ⬜ |
| P3.4 | Conservación S | **S†η_out S = η_in** (3 col., métricas in/out SEPARADAS) + balance de flujo región (paquete) + flujo canal a canal | S 3×3; baseline λ=1 (v1.0 bdg_smatrix) → Padé (solver aumentado) | C1 max\|R_S\|; C2 balance paquete **con término de defecto**; C3=diagnóstico de C1 | A1–A4 / R1–R3 (diseño) | ◐ BASELINE λ=1 EJECUTADO+AUDITADO (`gate_P3_4_smatrix.py` + `gate_P3_4_controls.py`). **ESTADOS SEPARADOS (ChatGPT)**: **P3.4-C2 = VALIDADO (nivel operador, alcance documentado)**: balance cierra con defecto medido; barrido nt→O(dt²) pend 2.00 (`gate_P3_4_c2dt.py`); **barrido N HECHO** (`gate_P3_4_closecontrols.py`): ΔQ_R/flujo convergen, ∫defecto→0 O(dx²), |balance| en piso temporal ~1.3e−7 indep. N, recon máquina; estado inicial = estado de PRUEBA (no paquete BdG). **Umbral 1e−8: S y recuento INVARIANTES** (thr 1e−6…1e−12, S−Sref=0.0) → clasificación no frágil. **P3.4-C1 = ejecutado, NO validado por convergencia**: pseudo-unit ~5e−5 (<1e−3 pero floor, dx-convergencia NO demostrada), X0 sin meseta (ventana + blow-up) → transfer INSUFICIENTE para A1/A3. Controles A/B PASAN: S_gate≡B.smatrix (0.0 exacto) + invariante a escala evanescente (~1e−12) → el floor viene de la PROPAGACIÓN por transferencia, no del matching. **P3.4 ABIERTO**, no cerrado. Pendientes (ChatGPT): barrido nt (orden temporal C2), barrido N (ΔQ_R), sensib. umbral 1e−8, paquete BdG asintótico; luego quizá solver scattering estabilizado. NO relajar umbrales. λ≠1 DEPENDE del solver Padé aumentado |
| P3.5 | Convergencia + coste | N, L, resolución ω; coste/freq vs denso | barridos | estable; coste medido | — | ⬜ |

**Nota P3.1 (bug cazado y corregido 2026-11-15):** 1ª versión usaba `np.linspace(endpoint=True)` para la
FFT + `np.interp` → error O(dx) en la REFERENCIA que enmascaraba la convergencia (daba pendiente 0.85).
Corregido con mallas alineadas por `arange` (espaciado exacto dx, sin interp) → pendiente 2.00 limpia.
Lección: verificar que la referencia no introduce su propio error de discretización.
**P3.1 APROBADO por ChatGPT (2026-11-15, revisión línea por línea)** con 3 ajustes documentales
aplicados: (1) Robin es de **1er orden, APROXIMADA** para las colas polinómicas de A³ (NO transparente
exacta; discrepancia rel. O(1/L)); (2) truncamiento de caja **O(L²e^{−κL})** (no O(e^{−κL})); (3) malla
**simétrica** (N impar, ambos extremos ±L). Mejora de rendimiento: LU factorizado una vez y reutilizado
en las 3 aplicaciones de A. Nota para P3.2-B: Robin simétrica respecto al peso trapezoidal
`W=diag(½,1,…,1,½)` (`WM=M†W`), NO al euclídeo ordinario.

**Nota P3.2-A (bug cazado y corregido 2026-11-15):** la referencia densa aplicaba `inv(2G+D)` sin el
proyector de borde, mientras las rutas dispersas anulan el rhs en Dirichlet (`φ_bdy=f_bdy≠0` vs `=0`)
→ `aug~denso=1.0` (100%) con campo no nulo en el borde. Corregido con `A_eff=Adense@Pbc` (misma BC en
todas las rutas). Misma lección que P3.1: las rutas comparadas deben usar BC idénticas. Alcance de
P3.2-A: SOLO equivalencia algebraica (ChatGPT); NO valida scattering ni pseudo-herm.

**Nota P3.2-B2 (veredicto ChatGPT 2026-11-15):** B2.1 APROBADO (K_P†=K_P se deduce de D=Dᵀ, A=Aᵀ).
B2.2 queda como **resultado OBSERVADO, no identidad SBP demostrada**: el flujo de B2.2 usa H trapezoidal
sobre la región y el de B2.3 el euclídeo global ×dx; no son automáticamente la misma formulación. Para
reclamar identidad SBP formal habría que definir `H_R=dx·diag(½,1,…,1,½)` sobre los nodos de la región y
DEMOSTRAR flujo(H_R)=diferencia de concomitantes; mientras sea solo numérico se describe como
comprobación de consistencia.

**DECISIÓN (usuario + ChatGPT 2026-11-15): OPCIÓN B — B2 es un gate NUMÉRICO de alcance limitado, NO
demostración SBP formal.** Condiciones documentales que quedan fijadas:
  1. B2.2 figura como **comprobación numérica de consistencia**, NO como identidad SBP demostrada.
  2. B2.3 figura como **convergencia al concomitante continuo para G=1, Dirichlet-fantasma, región [−6,6]**.
  3. NO se generaliza el resultado a otras BC, perfiles inhomogéneos ni al operador completo M_eff.
  4. **B3 debe traer sus PROPIOS controles de pseudo-hermiticidad** y NO apoyarse en ninguna supuesta
     demostración SBP (que no existe).
Pendiente antes de arrancar B3: ChatGPT revisa el código completo de `gate_P3_2B2c_continuum.py`
(g₁,g₂,g₃, asignación de derivadas a los fuentes, derivadas 2ª/3ª por EDO) → visto bueno definitivo.

**Nota P3.2-B2c (cierre de B2.3, 2026-11-15):** control pedido por ChatGPT. Se construye el concomitante
CONTINUO con una referencia independiente de la malla FD (convolución con las Green de línea
g₁,g₂,g₃ de (2−∂²)⁻ⁿ; derivadas sobre el fuente lorentziano analítico → núcleo continuo, sin aliasing;
2ª/3ª derivadas por recursión de la EDO). Converge a `BF_cont=−3.3390424e−3` (7 díg. estable al refinar
R,h). Resultado: **BF_disc→BF_cont a O(dx²) (pend. 2.00)** y **F_matriz→BF_cont a O(dx) (pend. 0.94)**,
el O(dx) debido al proyector sharp P_R (ya diagnosticado), no a fallo de la correspondencia con Green.
Primer intento (espectral/FFT) DESCARTADO: los símbolos de derivada alta amplificaban el aliasing de las
colas 1/x² y la referencia no convergía (lección: no usar FFT periódica con fuentes de cola polinómica).
Alcance: G=1 homogéneo, Dirichlet-fantasma, región [−6,6]. NO cierra B2.2 (identidad SBP formal).
**APROBADO por ChatGPT (2026-11-15, revisión matemática estática —núcleos g₁₂₃, asignación de fuentes,
recurrencias EDO, definición de los 3 errores y ajuste log-log de 4 puntos—; NO reprodujo la ejecución).**
Ajustes aplicados tras su revisión: κ=√(2G) derivada de G + assert de homogeneidad (no fijada a mano);
malla de cuadratura `linspace(-R,R,n+1)` con paso real; import `lam` eliminado. Controles SEPARADOS
añadidos y pasados: (a) h-convergencia con R=70 fijo → Δ con ratio ~4.0 (cuadratura O(h²) controlada);
(b) R-convergencia con h=6.25e−4 fijo → Δ~1e−18 desde R=30 (truncamiento de dominio a nivel de máquina,
núcleos ~e^{−√2 R}). Autoverificaciones independientes propias: g₂≡g₁⋆g₁, g₃≡g₁⋆g₂ (~cuadratura);
q₂ directa≡anidada (1.7e−9); p″ por EDO ≡ −q₁ (1.9e−9). **B2c CERRADO.**

**Solo tras P3.1–P3.4:** medir n_L(λ), G fijo, λ∈[0.5,2] → ¿persiste el colapso? (pregunta física).

## Precauciones (ChatGPT)
- A³: comprobar dependencia de BC y L, no solo coincidencia puntual. (hecho en P3.1)
- Pseudo-herm: medir sobre el sistema AUMENTADO completo, con todos los auxiliares.
- hermiticidad ≠ corrección física; λ=1 valida infraestructura, no los términos Padé; 0.5% Weyl-Padé
  es diagnóstico, no cota sobre observables. NO tocar v1.0.
- **Robin discreto NO es simétrico en el producto euclídeo ordinario** (ChatGPT): la discretización de
  extremos requiere el peso de borde correspondiente. No afecta a P3.2-A (algebraico), pero SÍ a las
  afirmaciones de autoadjunción en P3.2-B — tenerlo presente al medir `‖M_eff†η−ηM_eff‖`.
