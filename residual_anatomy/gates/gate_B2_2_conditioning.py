"""
BLOQUE 2 — Gate B2.2: ESTUDIO DE CONDICIONAMIENTO de M_aug (calibracion).

Implementa el diseno congelado en STABILIZED_SOLVER_BLOCK2_2.md. Objetivo de ESTA ejecucion:
CALIBRAR las constantes provisionales con solves de control; NO interpretar fisica (ni n_L, ni
pseudo-unitaridad como criterio).

Reutiliza el ensamblado YA AUDITADO de B2.1 por import (no lo modifica):
  assemble, build_Xi, rhs_for  de gate_B2_1_algebraic.

Tres representaciones (diseno §2): cruda -> Ruiz -> fisica (ell = xi_u = 1/sqrt(2) fija).
Metricas (2-norma): cond, drop, sigma_min/zona de rango, eta_back, kappa_out,b y kappa_out,M por
salida (adjunto con selector D_c^dag c), acuerdo de S entre representaciones, convergencia Delta_S en N.

Ejes: X0 FIJO, barrido en N (el eje X0 = gate aparte). Ventana 3x3 (guard de B2.1).
"""
import numpy as np
import scipy.linalg as la
import gate_B2_1_algebraic as g1   # ensamblado auditado (import = solo lectura)

ELL = 1.0/np.sqrt(2.0)             # xi_u = 1/sqrt(2 g n_u), g=n_u=1 (congelado)
EPS = np.finfo(float).eps


# ----------------------------- equilibrados -----------------------------------
def ruiz_scaling(M, iters=30):
    """Equilibrado iterativo de Ruiz. Devuelve (dr, dc) con M_ruiz = dr[:,None]*M*dc[None,:]."""
    n = M.shape[0]
    dr = np.ones(n); dc = np.ones(n); A = M.copy()
    for _ in range(iters):
        r = np.sqrt(np.maximum(np.max(np.abs(A), axis=1), 1e-300))
        c = np.sqrt(np.maximum(np.max(np.abs(A), axis=0), 1e-300))
        A = (A / r[:, None]) / c[None, :]
        dr /= r; dc /= c
    return dr, dc


def physical_scaling(info):
    """Escalado FISICO coherente con Ytilde=(u, ell u', v, ell v') aplicado a incognitas Y ecuaciones.
    Columnas (incognitas): campo por nodo diag(1,1/ell,1,1/ell) [Y = D_c z]; amplitudes a_L,a_R -> 1.
    Filas (ecuaciones): bloques interior+matching por diag(1,ell,1,ell); restricciones -> 1.
    Devuelve (dr, dc) diagonales. ell FIJA, independiente de N."""
    ndof = info["ndof"]; N = info["N"]
    dc = np.ones(ndof); dr = np.ones(ndof)
    nfield = 4*(N+1)
    # columnas de campo
    for j in range(N+1):
        dc[4*j:4*j+4] = [1.0, 1.0/ELL, 1.0, 1.0/ELL]
    # a_L, a_R (ultimas 8 columnas) quedan en 1
    # filas: interior (4N) + matching (8) en bloques de 4 ; restricciones (4) en 1
    nrows_blocks = 4*N + 8
    for rb in range(0, nrows_blocks, 4):
        dr[rb:rb+4] = [1.0, ELL, 1.0, ELL]
    # ultimas 4 filas (restricciones) quedan en 1
    return dr, dc


# ----------------------------- metricas ---------------------------------------
def out_indices(info):
    """Indices globales de las 3 amplitudes de SALIDA (a_out) en el vector de incognitas fisico."""
    labL, labR = info["labL"], info["labR"]
    iaL, iaR = info["iaL"], info["iaR"]
    oL = iaL + labL.index("out")                      # a_out^L
    oR = [iaR + i for i, L in enumerate(labR) if L == "out"]  # a_out1^R, a_out2^R
    return [oL, oR[0], oR[1]]


def smatrix_and_solves(info, M, dr, dc):
    """Resuelve los 3 problemas en la representacion (dr,dc). Devuelve:
       S (3x3, convencion fisica), lista de (z, x, b_rep) por problema."""
    outs = out_indices(info)
    Mrep = dr[:, None] * M * dc[None, :]
    sols = []
    S = np.zeros((3, 3), dtype=complex)
    for p, sin in enumerate((np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1]))):
        b = g1.rhs_for(info, sin)
        brep = dr * b
        z = la.solve(Mrep, brep)
        x = dc * z                                    # amplitudes en convencion FISICA
        sols.append((z, x, brep, Mrep))
        for o, idx in enumerate(outs):
            S[o, p] = x[idx]
    return S, sols, outs, Mrep


def kappa_out(Mrep, dc, outs, sols, Mnorm):
    """Sensibilidades por SALIDA y RHS (2-norma). Selector en vars nuevas = D_c^dag c.
    Registra, por (problema, salida): |s|, sens ABS a db (‖y‖), sens ABS a dM (A_out,M=‖y‖‖M‖‖z‖),
    sens REL kappa_out,b y kappa_out,M, y el control adjunto delta_adj=|y^dag b~ - x_salida|.
    Devuelve dict con el detalle completo + agregados (max/min)."""
    adj = {}
    for idx in outs:                                  # un adjunto por salida: M~^dag y = D_c^dag c
        c = np.zeros(Mrep.shape[0], dtype=complex); c[idx] = 1.0
        crep = np.conj(dc) * c
        adj[idx] = la.solve(Mrep.conj().T, crep)
    detail = []
    for pi, (z, x, brep, _Mrep) in enumerate(sols):
        for idx in outs:
            y = adj[idx]
            s = np.vdot(y, brep)                       # debe == x[idx]
            da = abs(s - x[idx])                       # control adjunto (abs)
            absb = la.norm(y)                          # sens ABS a db:  |ds| <= ‖y‖‖db‖
            AM = la.norm(y) * Mnorm * la.norm(z)       # sens ABS a dM:  |ds| <= ‖y‖‖M‖‖z‖·(‖dM‖/‖M‖)
            kb = abs(absb * la.norm(brep) / s) if abs(s) > 0 else np.nan
            kM = abs(AM / s) if abs(s) > 0 else np.nan
            drel = (da/abs(s)) if abs(s) > 0 else np.nan
            detail.append(dict(p=pi, idx=idx, s=abs(s), absb=absb, AM=AM,
                               kb=kb, kM=kM, dadj=da, dadj_rel=drel))
    fin = lambda key: [d[key] for d in detail if np.isfinite(d[key])]
    return dict(detail=detail,
                dadj_abs=max(d["dadj"] for d in detail),
                dadj_rel=(max(fin("dadj_rel")) if fin("dadj_rel") else np.nan),
                kb_max=(max(fin("kb")) if fin("kb") else np.nan),
                kM_max=(max(fin("kM")) if fin("kM") else np.nan),
                absb_max=max(d["absb"] for d in detail),
                AM_max=max(d["AM"] for d in detail),
                s_min=min(d["s"] for d in detail))


def eta_back(Mrep, sols, Mnorm):
    out = []
    for (z, x, brep, _Mrep) in sols:
        out.append(la.norm(Mrep@z - brep) / (Mnorm*la.norm(z) + la.norm(brep)))
    return max(out)


def rank_zone(smax_rep, ndof, c4=1.0):
    smin, smx = smax_rep[-1], smax_rep[0]
    floor = c4 * ndof * EPS * smx
    if smin > 1e3 * floor:   zone = "RESUELTO"
    elif smin > floor:        zone = "ALERTA"
    else:                     zone = "NO-RESUELTO"
    return zone, smin, smx, floor


# ----------------------------- barrido ----------------------------------------
def analyze(w, N, X0=12.0):
    M, info = g1.assemble(w, N, X0)
    ndof = info["ndof"]
    # guard de ventana 3x3 (reutiliza criterio de B2.1)
    labL, labR = info["labL"], info["labR"]
    window = (labL.count("in")==1 and labL.count("out")==1 and labL.count("dec")==1 and
              labL.count("grow")==1 and labR.count("in")==2 and labR.count("out")==2 and
              labR.count("dec")==0 and labR.count("grow")==0)
    if not window:
        return dict(window=False, N=N, w=w)
    reps = {}
    dcr = dict(cruda=(np.ones(ndof), np.ones(ndof)))
    dr_r, dc_r = ruiz_scaling(M);      dcr["ruiz"]   = (dr_r, dc_r)
    dr_f, dc_f = physical_scaling(info); dcr["fisica"] = (dr_f, dc_f)
    S_by_rep = {}
    for name, (dr, dc) in dcr.items():
        S, sols, outs, Mrep = smatrix_and_solves(info, M, dr, dc)
        sv = la.svdvals(Mrep); Mnorm = sv[0]
        cond = sv[0]/sv[-1]
        eb = eta_back(Mrep, sols, Mnorm)
        agg = kappa_out(Mrep, dc, outs, sols, Mnorm)
        zone, smin, smx, floor = rank_zone(sv, ndof)
        reps[name] = dict(cond=cond, sigma_min=sv[-1], sigma_max=sv[0], eta_back=eb,
                          kappa_b_max=agg["kb_max"], kappa_M_max=agg["kM_max"],
                          abs_b_max=agg["absb_max"], absM_max=agg["AM_max"],
                          s_min=agg["s_min"], dadj_abs=agg["dadj_abs"], dadj_rel=agg["dadj_rel"],
                          rank_zone=zone, rank_floor=floor, detail=agg["detail"])
        S_by_rep[name] = S
    cond_cruda = reps["cruda"]["cond"]
    for name in reps:
        reps[name]["drop"] = cond_cruda/reps[name]["cond"]
    # acuerdo de S entre representaciones (max abs elemento a elemento)
    Sc = S_by_rep["cruda"]
    agree = {name: float(np.max(np.abs(S_by_rep[name]-Sc))) for name in ("ruiz","fisica")}
    return dict(window=True, N=N, w=w, reps=reps, S=S_by_rep, agree=agree)


def main(w=0.10, X0=12.0, Ns=(40,60,80,120,160)):
    print("="*94)
    print(f"GATE B2.2 — CONDICIONAMIENTO (calibracion)   w={w}  X0={X0}  ell=1/sqrt2={ELL:.4f}")
    print("="*94)
    results = []
    for N in Ns:
        r = analyze(w, N, X0)
        if not r["window"]:
            print(f"  N={N}: FUERA DE VENTANA 3x3 -> saltado"); continue
        results.append(r)
        print(f"\n--- N={N}  (ndof={4*N+12}) ---")
        print(f"  {'repr':>7} {'cond2':>11} {'drop':>8} {'sigma_min':>11} {'eta_back':>10} "
              f"{'kap_out,b':>11} {'kap_out,M':>11} {'rango':>11}")
        for name in ("cruda","ruiz","fisica"):
            d = r["reps"][name]
            print(f"  {name:>7} {d['cond']:11.3e} {d['drop']:8.2f} {d['sigma_min']:11.3e} "
                  f"{d['eta_back']:10.2e} {d['kappa_b_max']:11.3e} {d['kappa_M_max']:11.3e} "
                  f"{d['rank_zone']:>11}")
        for name in ("cruda","ruiz","fisica"):
            d = r["reps"][name]
            print(f"      {name:>6}: delta_adj(abs/rel)={d['dadj_abs']:.2e}/{d['dadj_rel']:.2e}  "
                  f"abs_sens_b(‖y‖)={d['abs_b_max']:.2e}  abs_sens_M={d['absM_max']:.2e}")
        print(f"  acuerdo S (max|S_rep-S_cruda|):  ruiz={r['agree']['ruiz']:.3e}  "
              f"fisica={r['agree']['fisica']:.3e}")
        print(f"  |S_cruda| Frobenius = {la.norm(r['S']['cruda']):.4f}   s_min(out) = "
              f"{min(r['reps'][n]['s_min'] for n in r['reps']):.2e}")
    # ley de crecimiento cond ~ C N^p (sobre cada representacion) y Delta_S convergencia
    if len(results) >= 3:
        Ns_ok = np.array([r["N"] for r in results], float)
        print("\n" + "-"*94)
        print("  AJUSTE cond(N) ~ C * N^p   (log-log, minimos cuadrados):")
        for name in ("cruda","ruiz","fisica"):
            conds = np.array([r["reps"][name]["cond"] for r in results])
            p, logC = np.polyfit(np.log(Ns_ok), np.log(conds), 1)
            print(f"    {name:>7}:  p = {p:+.3f}   C = {np.exp(logC):.3e}   "
                  f"(cond {conds[0]:.2e} -> {conds[-1]:.2e})")
        print("\n  CONVERGENCIA Delta_S (repr=fisica, vs malla mas fina N="
              f"{results[-1]['N']}):")
        Sfin = results[-1]["S"]["fisica"]
        den = max(1.0, la.norm(Sfin))
        for r in results[:-1]:
            dS = la.norm(r["S"]["fisica"] - Sfin)/den
            print(f"    N={r['N']:4d} vs {results[-1]['N']}:  Delta_S = {dS:.3e}")
        print("  Delta_S entre mallas consecutivas:")
        for a in range(len(results)-1):
            dS = la.norm(results[a]["S"]["fisica"]-results[a+1]["S"]["fisica"])/max(1.0, la.norm(results[a+1]["S"]["fisica"]))
            print(f"    N={results[a]['N']:4d} vs {results[a+1]['N']:4d}:  Delta_S = {dS:.3e}")
    print("\n" + "="*94)
    print("  NOTA: calibracion de constantes provisionales. NO valida fisica, S, ni pseudo-unitaridad.")
    print("="*94)
    return results


if __name__ == "__main__":
    main(w=0.10, X0=12.0, Ns=(40, 60, 80, 120, 160))
    print("\n\n########## control: segunda frecuencia en ventana ##########")
    main(w=0.05, X0=12.0, Ns=(40, 60, 80, 120, 160))
