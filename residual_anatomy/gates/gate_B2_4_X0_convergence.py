"""
BLOQUE 2 — Gate B2.4: CONVERGENCIA en X0 (diseno STABILIZED_SOLVER_BLOCK2_4.md).

Eje X0 AISLADO: varia X0 con esquema de caja 2o orden, dx FIJO, representacion fisica (B2.2),
clasificacion/normalizacion de B2.1. NO solve disperso, NO cambia arquitectura.

ESTE ARCHIVO, en __main__, SOLO ejecuta:
  (1) validacion del ESTIMADOR (cond, sigma_min, zona de rango) vs SVD densa en caso moderado,
  (2) prueba de COSTE a N grande,
NO el barrido completo (queda para despues de auditar el codigo).

Reutiliza por import B2.1 (assemble, A_at/A_homog, nL..) y B2.2 (physical_scaling, smatrix_and_solves,
kappa_out, eta_back, rank_zone) SIN modificarlos.
"""
import time
import resource
import numpy as np
import scipy.linalg as la
import gate_B2_1_algebraic as g1
import gate_B2_2_conditioning as g2

RNG = np.random.default_rng(0)


# ------------------- construccion de la representacion fisica ------------------
def build_phys(w, N, X0):
    """Devuelve (Mrep, info, dr, dc). Mrep = D_r M D_c (representacion fisica). Verifica ventana 3x3."""
    M, info = g1.assemble(w, N, X0)
    labL, labR = info["labL"], info["labR"]
    window = (labL.count("in")==1 and labL.count("out")==1 and labL.count("dec")==1 and
              labL.count("grow")==1 and labR.count("in")==2 and labR.count("out")==2 and
              labR.count("dec")==0 and labR.count("grow")==0)
    if not window:
        raise RuntimeError(f"w={w}, X0={X0} fuera de ventana 3x3")
    dr, dc = g2.physical_scaling(info)
    Mrep = dr[:, None] * M * dc[None, :]
    return Mrep, info, dr, dc


# ------------------- metricas de condicionamiento ------------------------------
def svd_metrics(Mrep, ndof):
    """Referencia EXACTA por SVD densa: sigma_min, sigma_max, cond, zona de rango."""
    sv = la.svdvals(Mrep)
    zone, smin, smax, floor = g2.rank_zone(sv, ndof)
    return dict(sigma_min=sv[-1], sigma_max=sv[0], cond=sv[0]/sv[-1], zone=zone, floor=floor)


def _unit(n):
    v = RNG.standard_normal(n) + 1j*RNG.standard_normal(n)
    return v/la.norm(v)


def sigma_max_power(Mrep, iters=300, tol=1e-12):
    """sigma_max por iteracion de potencia sobre M^H M (solo matvecs).
       Devuelve (sigma_max, diag) con diag=dict(iters, converged, crit)."""
    v = _unit(Mrep.shape[0]); lam = 0.0; conv = False; used = 0; crit = np.inf
    for k in range(1, iters+1):
        w = Mrep.conj().T @ (Mrep @ v)
        nw = la.norm(w)
        if not np.isfinite(nw):
            return np.nan, dict(iters=k, converged=False, crit=np.inf, fail="nw no finito")
        if nw == 0:
            # M^H M v = 0: el vector cayo en el nucleo; NO prueba sigma_max(M)=0 -> no concluyente
            return np.nan, dict(iters=k, converged=False, crit=np.inf,
                                fail="M^H M v = 0 (nucleo/degeneracion, no concluyente)")
        v = w/nw; used = k
        crit = abs(nw - lam)/nw
        if crit <= tol:
            conv = True; break
        lam = nw
    Mv = Mrep @ v
    val = np.sqrt(np.vdot(Mv, Mv).real)          # Rayleigh: sigma_max^2 = v^H M^H M v (v unit)
    return val, dict(iters=used, converged=conv, crit=float(crit))


def sigma_min_invpower(Mrep, lu, iters=300, tol=1e-12):
    """sigma_min por iteracion inversa sobre (M^H M)^{-1} usando la LU de M (dos solves/iter).
       Devuelve (sigma_min, diag) con diag=dict(iters, converged, crit)."""
    v = _unit(Mrep.shape[0]); lam = 0.0; conv = False; used = 0; crit = np.inf
    for k in range(1, iters+1):
        t = la.lu_solve(lu, v, trans=2)          # M^{-H} v
        w = la.lu_solve(lu, t, trans=0)          # M^{-1} t = (M^H M)^{-1} v
        nw = la.norm(w); used = k
        if nw == 0 or not np.isfinite(nw):
            return np.nan, dict(iters=k, converged=False, crit=np.inf, fail="nw cero o no finito")
        v = w/nw
        crit = abs(nw - lam)/nw
        if crit <= tol:
            conv = True; break
        lam = nw
    t = la.lu_solve(lu, v, trans=2); w = la.lu_solve(lu, t, trans=0)
    lam = np.vdot(v, w).real                      # v^H (M^H M)^{-1} v (v unit) = lambda_max
    if not np.isfinite(lam) or lam <= 0:
        return np.nan, dict(iters=used, converged=False, crit=float(crit),
                            fail="lambda final no finito o <=0")
    return 1.0/np.sqrt(lam), dict(iters=used, converged=conv, crit=float(crit))


def estim_metrics(Mrep, ndof, lu=None, strict=False):
    """Estimador: sigma_max (potencia), sigma_min (inversa via LU), cond, zona de rango.
       Reusa la LU si se pasa; si no, factoriza. Propaga diagnostico de convergencia.
       Si strict=True y alguna iteracion NO convergio -> RuntimeError (para el barrido)."""
    if lu is None:
        lu = la.lu_factor(Mrep)
    smax, dmax = sigma_max_power(Mrep)
    smin, dmin = sigma_min_invpower(Mrep, lu)
    # validez: ambos estimados finitos y estrictamente positivos (ademas de convergidos)
    valid = (np.isfinite(smax) and np.isfinite(smin) and smax > 0 and smin > 0)
    converged = bool(dmax["converged"] and dmin["converged"] and valid)
    if strict and not converged:
        raise RuntimeError(f"estimador NO fiable: valid={valid}, smax {dmax}, smin {dmin}")
    if not valid:
        # NO convertir silenciosamente una estimacion invalida en NO-RESUELTO
        return dict(sigma_min=smin, sigma_max=smax, cond=np.nan, zone="FALLO", floor=np.nan,
                    converged=False, diag_max=dmax, diag_min=dmin)
    floor = 1.0 * ndof * g2.EPS * smax
    if   smin > 1e3*floor: zone = "RESUELTO"
    elif smin > floor:     zone = "ALERTA"
    else:                  zone = "NO-RESUELTO"
    return dict(sigma_min=smin, sigma_max=smax, cond=smax/smin, zone=zone, floor=floor,
                converged=converged, diag_max=dmax, diag_min=dmin)


# ------------------- error de frontera ----------------------------------------
def E_borde(w, X0):
    """E_borde(+-X0) = ||A_at(+-X0) - A_homog(asint)|| / ||A_homog(asint)||  (2-norma)."""
    AhL = g1.A_homog(g1.nL, g1.vL, g1.VL, w)
    AhR = g1.A_homog(g1.nR, g1.vR, g1.VR, w)
    eL = la.norm(g1.A_at(-X0, w) - AhL, 2)/la.norm(AhL, 2)
    eR = la.norm(g1.A_at(+X0, w) - AhR, 2)/la.norm(AhR, 2)
    return eL, eR


# ------------------- (1) VALIDACION del estimador ------------------------------
def validate_estimator(w=0.10, X0=12.0):
    print("="*80)
    print(f"[VALIDACION estimador vs SVD densa]   w={w}  X0={X0}")
    print("="*80)
    N = int(round(2*X0/0.02))
    Mrep, info, dr, dc = build_phys(w, N, X0)
    ndof = info["ndof"]
    t0 = time.time(); exact = svd_metrics(Mrep, ndof); t_svd = time.time()-t0
    t0 = time.time(); lu = la.lu_factor(Mrep); t_lu = time.time()-t0
    t0 = time.time(); est = estim_metrics(Mrep, ndof, lu=lu); t_est = time.time()-t0
    print(f"  N={N} ndof={ndof}")
    print(f"  {'':>10} {'sigma_min':>13} {'sigma_max':>13} {'cond':>13} {'zona':>12}")
    print(f"  {'SVD densa':>10} {exact['sigma_min']:13.6e} {exact['sigma_max']:13.6e} "
          f"{exact['cond']:13.6e} {exact['zone']:>12}")
    print(f"  {'estimador':>10} {est['sigma_min']:13.6e} {est['sigma_max']:13.6e} "
          f"{est['cond']:13.6e} {est['zone']:>12}")
    rel = lambda a, b: abs(a-b)/abs(b)
    print(f"  err rel: sigma_min={rel(est['sigma_min'],exact['sigma_min']):.2e}  "
          f"sigma_max={rel(est['sigma_max'],exact['sigma_max']):.2e}  "
          f"cond={rel(est['cond'],exact['cond']):.2e}  "
          f"zona_coincide={est['zone']==exact['zone']}")
    print(f"  convergencia estimador: smax={est['diag_max']}  smin={est['diag_min']}")
    print(f"  tiempo: SVD={t_svd:.2f}s | LU={t_lu:.2f}s + estimador={t_est:.2f}s "
          f"-> estimador_total={t_lu+t_est:.2f}s")
    return exact, est


# ------------------- (2) PRUEBA de coste a N grande ---------------------------
def cost_probe(w=0.10, X0=28.0):
    print("\n" + "="*80)
    print(f"[PRUEBA de coste]   w={w}  X0={X0}  (N grande, solo estimador + E_borde)")
    print("="*80)
    N = int(round(2*X0/0.02))
    ndof = 4*N + 12
    mem_one = (ndof**2 * 16) / 1e9
    print(f"  N={N} ndof={ndof}  memoria de UNA matriz densa compleja ~ {mem_one:.2f} GB "
          f"(NO es el pico del proceso)")
    t0 = time.time(); Mrep, info, dr, dc = build_phys(w, N, X0); t_build = time.time()-t0
    t0 = time.time(); lu = la.lu_factor(Mrep); t_lu = time.time()-t0
    t0 = time.time(); est = estim_metrics(Mrep, ndof, lu=lu); t_est = time.time()-t0
    eL, eR = E_borde(w, X0)
    # pico real de memoria del proceso (macOS: ru_maxrss en bytes)
    peak_gb = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    print(f"  tiempos: build={t_build:.1f}s  LU={t_lu:.1f}s  estimador={t_est:.1f}s")
    print(f"  pico RSS del proceso (ru_maxrss, macOS=bytes) ~ {peak_gb:.2f} GB "
          f"(incluye matriz + factores LU + temporales; coexisten durante lu_factor)")
    print(f"  estimador: sigma_min={est['sigma_min']:.3e} sigma_max={est['sigma_max']:.3e} "
          f"cond={est['cond']:.3e} zona={est['zone']}  convergido={est['converged']}")
    print(f"    convergencia: smax={est['diag_max']}  smin={est['diag_min']}")
    print(f"  E_borde(+-X0): izq={eL:.3e} der={eR:.3e}")
    return est


# ------------------- (3) BARRIDO en X0 (estrategia (a): rutas auditadas) -------
# Condiciones: B2.1/2/3 intactos; estimador con strict=True; fallos EXPLICITOS (no NO-RESUELTO);
# dx principal (0.02) y control (0.012) SEPARADOS; cond matricial != convergencia fisica de S.
def analyze_point(w, N, X0, strict=True):
    """Metricas de un punto (w,N,X0) por rutas AUDITADAS (refactorizan; opcion (a)).
       estim_metrics(strict=True): un estimador no fiable -> RuntimeError (fallo explicito)."""
    M, info = g1.assemble(w, N, X0)
    labL, labR = info["labL"], info["labR"]
    window = (labL.count("in")==1 and labL.count("out")==1 and labL.count("dec")==1 and
              labL.count("grow")==1 and labR.count("in")==2 and labR.count("out")==2 and
              labR.count("dec")==0 and labR.count("grow")==0)
    if not window:
        raise RuntimeError(f"w={w}, X0={X0} fuera de ventana 3x3")
    dr, dc = g2.physical_scaling(info)
    Mrep = dr[:, None] * M * dc[None, :]
    ndof = info["ndof"]
    cm = estim_metrics(Mrep, ndof, strict=strict)           # LU propia (refactoriza)
    Mnorm = cm["sigma_max"]
    S, sols, outs, _ = g2.smatrix_and_solves(info, M, dr, dc)  # 3 solves (refactorizan)
    eb = g2.eta_back(Mrep, sols, Mnorm)
    kap = g2.kappa_out(Mrep, dc, outs, sols, Mnorm)          # 3 adjuntos (refactorizan)
    eL, eR = E_borde(w, X0)
    # residuos de matching y de la restriccion 'grow', recuperados de x+info (NO modifica B1/B2):
    # las filas de matching del M original son Y_b - Xi_b a_b; grow es a_L[grow]=0 impuesto.
    iaL, iaR, iyN = info["iaL"], info["iaR"], info["iyN"]
    XiL, XiR = info["XiL"], info["XiR"]
    igrowL = info["labL"].index("grow")
    # Residuos de matching en DOS juegos de variables:
    #  - *_abs/_rel: VARIABLES ORIGINALES (abs y relativa por ||Y_b||).
    #  - *_sc_abs/_sc_rel: variables REESCALADAS del diseno B2.4, T=diag(1,ell,1,ell), ell=xi_u=1/sqrt2.
    T = np.array([1.0, g2.ELL, 1.0, g2.ELL])
    match_abs = match_rel = grow_resid = 0.0
    match_sc_abs = match_sc_rel = 0.0
    for (z, x, brep, _Mrep) in sols:
        Y0 = x[0:4]; YN = x[iyN:iyN+4]; aL = x[iaL:iaL+4]; aR = x[iaR:iaR+4]
        resL = Y0 - XiL@aL; resR = YN - XiR@aR
        rL = la.norm(resL); rR = la.norm(resR)
        match_abs = max(match_abs, rL, rR)
        match_rel = max(match_rel, rL/max(la.norm(Y0), 1e-300), rR/max(la.norm(YN), 1e-300))
        tL = la.norm(T*resL); tR = la.norm(T*resR)                 # ||T(Y_b - Xi_b a_b)||
        match_sc_abs = max(match_sc_abs, tL, tR)
        match_sc_rel = max(match_sc_rel, tL/max(la.norm(T*Y0), 1e-300), tR/max(la.norm(T*YN), 1e-300))
        # RESIDUO de la restriccion impuesta a_L[grow]=0; NO es contaminacion fisica (anulada por
        # construccion) y depende de la normalizacion modal fijada en B2.1.
        grow_resid = max(grow_resid, abs(aL[igrowL]))
    return dict(status="OK", X0=X0, N=N, S=S, cond=cm["cond"], sigma_min=cm["sigma_min"],
                sigma_max=cm["sigma_max"], zone=cm["zone"], converged=cm["converged"],
                eta_back=eb, kappa_b=kap["kb_max"], kappa_M=kap["kM_max"], dadj_abs=kap["dadj_abs"],
                match_err_abs=match_abs, match_err_rel=match_rel,
                match_err_sc_abs=match_sc_abs, match_err_sc_rel=match_sc_rel,
                grow_constraint_resid=grow_resid, EbL=eL, EbR=eR)


def sweep(w, X0_list=(8,10,12,14,16,20,24,28), dx=0.02, ref_X0=12, p_ref=0.97):
    """Barrido en X0 a dx FIJO. R_cond con p de B2.2 (comparacion heuristica, NO ley).
       Fallos se REGISTRAN explicitamente (status='FALLO'), no se reclasifican como NO-RESUELTO."""
    print("="*100)
    print(f"GATE B2.4 — BARRIDO en X0   w={w}  dx={dx}  (esquema caja 2o orden, eje X0)")
    print("="*100)
    pts = []
    for X0 in X0_list:
        N = int(round(2*X0/dx))
        try:
            pt = analyze_point(w, N, X0, strict=True)
        except (RuntimeError, np.linalg.LinAlgError) as e:
            # fallos PREVISTOS (estimador no fiable, ventana, singularidad LU) -> registro explicito,
            # conservando tipo y mensaje. NO se captura Exception generica (ocultaria bugs).
            pt = dict(status="FALLO", X0=X0, N=N, err=f"{type(e).__name__}: {e}")
            print(f"  X0={X0:4.1f} N={N:5d}  *** FALLO [{type(e).__name__}]: {e} ***")
        pts.append(pt)
    ok = [pt for pt in pts if pt["status"] == "OK"]
    # R_cond(X0) = [cond(X0)/cond(ref)]*(N(ref)/N(X0))^p   (p de B2.2; comparacion, no ley)
    ref = next((pt for pt in ok if pt["X0"] == ref_X0), None)
    Nref = int(round(2*ref_X0/dx))
    if ref is None:
        print(f"\n  AVISO: el punto de referencia X0={ref_X0} no esta disponible (fallo o no en la lista)"
              f" -> R_cond NO evaluable en este barrido.")
    for pt in ok:
        pt["R_cond"] = (pt["cond"]/ref["cond"]) * (Nref/pt["N"])**p_ref if ref is not None else float('nan')
        print(f"  X0={pt['X0']:4.1f} N={pt['N']:5d} cond={pt['cond']:.3e} R_cond={pt['R_cond']:.3f} "
              f"smin={pt['sigma_min']:.3e} {pt['zone']} eta={pt['eta_back']:.2e} "
              f"kap_b={pt['kappa_b']:.2e} kap_M={pt['kappa_M']:.2e} dadj={pt['dadj_abs']:.2e} "
              f"Eb={pt['EbL']:.1e}/{pt['EbR']:.1e}")
        print(f"        match_abs={pt['match_err_abs']:.2e} match_rel={pt['match_err_rel']:.2e} "
              f"match_sc_abs={pt['match_err_sc_abs']:.2e} match_sc_rel={pt['match_err_sc_rel']:.2e} "
              f"grow_constraint_resid={pt['grow_constraint_resid']:.2e}")
    # Delta_S: Frobenius Y max elemento a elemento; consecutivos y vs el X0 MAYOR (ref comparativa
    # provisional, NO limite fisico). La comparacion por elemento evita que la norma global
    # oculte un canal sensible.
    def fmt3(D):   # matriz 3x3 real en una linea: filas entre [] separadas por |
        return " | ".join("[" + " ".join(f"{D[i,j]:.1e}" for j in range(3)) + "]" for i in range(3))
    if len(ok) >= 2:
        print("\n  Delta_S(X0): Frob, max y matriz D_ij=|S_ij(a)-S_ij(b)| (abs; adecuada con amplitudes "
              "pequenas).  NO prueba convergencia fisica de S.")
        print("  'next_OK' = siguiente X0 con estado OK: puede SALTAR puntos fallidos (no siempre adyacente).")
        big = ok[-1]["S"]; dbig = max(1.0, la.norm(big))
        for i, pt in enumerate(ok):
            if i+1 < len(ok):
                Dc = np.abs(pt["S"]-ok[i+1]["S"])
                dF = la.norm(pt["S"]-ok[i+1]["S"])/max(1.0, la.norm(ok[i+1]["S"]))
                print(f"    X0={pt['X0']:4.1f} next_OK(vs {ok[i+1]['X0']:.0f}): Frob={dF:.3e} "
                      f"max={Dc.max():.3e}  D={fmt3(Dc)}")
            Dv = np.abs(pt["S"]-big); vF = la.norm(pt["S"]-big)/dbig
            print(f"    X0={pt['X0']:4.1f} vs {ok[-1]['X0']:.0f}:      Frob={vF:.3e} "
                  f"max={Dv.max():.3e}  D={fmt3(Dv)}")
    print("\n  NOTA: cond/R_cond caracterizan el condicionamiento; NO demuestran convergencia fisica de S,")
    print("        pseudo-unitaridad ni n_L. R_cond usa p=0.97 de B2.2 (comparacion heuristica, no ley).")
    return pts


if __name__ == "__main__":
    # Paso 1: validar el estimador en caso MODERADO (SVD factible)
    validate_estimator(w=0.10, X0=12.0)
    validate_estimator(w=0.05, X0=12.0)
    # Paso 2: coste a N grande (NO es el barrido; un punto)
    cost_probe(w=0.10, X0=28.0)
    print("\n  NOTA: barrido completo de X0 NO ejecutado (sweep() NO se llama aqui; requiere autorizacion).")
