"""
RESIDUAL ANATOMY — PASO 0 (gate de viabilidad numerica).

Objetivo: decidir si el residuo R(w,lam) = T_upstream/(kappa_local/2pi) - 1
          (w=2) tiene dependencia con lam distinguible del error numerico.

NO modifica la entrega congelada v1.0. Importa el solver raiz sin tocarlo:
  - correlations.scatter_field  (matching modal con knob lam, VALIDADO para beta/correlaciones)
  - mode_matching.solve_inject  (ruta del benchmark, lam=1) -> solo para control de regresion
Dependencias (SHA-256) registradas en residual_anatomy/README.md.

Observable = n_L (canal u upstream saliente de norma +1, Ec.42 FP), NUNCA n_tot.
Referencia = kappa_local geometrica (geom(), copiada verbatim de FP_benchmark_table.py),
             independiente de lam por construccion.
Estimador  = T_omega = omega/ln(1+1/n_L), T_0 = lim omega->0 por ajuste cuadratico
             (identico al benchmark; lin y cub solo para el spread por orden).

GATE: no fuerza PASS. Emite DOS veredictos separados:
  (i)  residuo total R vs presupuesto     -> ya detectado en el benchmark (-4.28% vs 0.06%)
  (ii) dependencia dispersiva dR(lam)=R(lam)-R(1) vs presupuesto (PAREADO por orden)
"""
import sys, os, csv, contextlib
sys.path.insert(0, os.path.join(os.path.dirname(os.path.dirname(os.path.abspath(__file__))), 'solver'))
import numpy as np
with contextlib.redirect_stdout(open(os.devnull, 'w')):     # silencia auto-tests de import
    import correlations as C
    import mode_matching as MM

HERE = os.path.dirname(os.path.abspath(__file__))
REF_TUP_OVER_KLOC = 0.9572413965031795        # out/fp_benchmark.csv, w=2 (congelado)
TOL_REGRESION = 1e-4

# ---- geom() COPIADO VERBATIM de FP_benchmark_table.py (geometrico, lam-indep) ----
g = C.g; n_u, mach_u, md = 1.0, 0.5, 2.0
vLb = mach_u*np.sqrt(g*n_u); jc = n_u*vLb; n_d = (jc/(md*np.sqrt(g)))**(2/3)
def box(w): L = max(40., 8*w+20.); return L, int(round(L*240))
def geom(w):
    L, _ = box(w); x = np.linspace(-L, L, 80000); s = 0.5*(1+np.tanh(x/w)); n = n_u+(n_d-n_u)*s
    c = np.sqrt(n); v = jc/n; f = c-v; i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    t = f[i]/(f[i]-f[i+1]); xH = x[i]+t*(x[i+1]-x[i]); cH = np.interp(xH, x, c)
    kK = abs(np.interp(xH, x, np.gradient(f, x)))
    Lam = 2*cH**2; dxi = 1.58*cH/(kK*Lam**2)**(1/3)
    kbar = abs(np.interp(xH+0.5*dxi, x, f)-np.interp(xH-0.5*dxi, x, f))/dxi
    return kK, kbar, cH, xH

# ---- extraccion de n_L por RUTA NUEVA (scatter_field, con knob lam) ----
def occ_new(w, om, L, N, lam):
    """n_L, n_tot, conservacion-de-norma-de-la-columna-inyectada (NO pseudo-u completa)."""
    op = 'validated' if lam == 1.0 else 'lambda'
    r = C.scatter_field(L, N, om, 'd2', w=w, lam=lam, op=op)
    if r['Rin_neg'] is None or r['nchan'] != 2:
        return None
    oa = r['out_amp']
    nL = abs(oa.get(('L', 1), 0.0))**2
    nRp = abs(oa.get(('R', 1), 0.0))**2
    cons = abs(r['flux_out'] - r['eta_in'])        # sum_out eta|a|^2 == eta_in  (columna d2)
    return nL, nL+nRp, cons

# ---- extraccion de n_L por RUTA DEL BENCHMARK (solve_inject, lam=1) ----
def occ_bench(w, om, L, N):
    MM.w = w; MM.omega = om; x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
    Rin = [m for m in MM.modes_full(n[-1], v[-1], V[-1]) if m['prop'] and m['vg'] < 0]
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
    out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=incR)
    nL = sum(a for (q, s, a) in out if q > 0 and s == 'L')
    nRp = sum(a for (q, s, a) in out if q > 0 and s == 'R')
    pu = abs(sum(q*a for (q, s, a) in out)+1.0)
    return nL, nL+nRp, pu

# ---- extrapolacion T0 identica al benchmark ----
XR = np.array([0.05, 0.065, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.27])
def T0_orders(w, L, N, lam, bench=False):
    """Devuelve (T0_quad, T0_lin, T0_cub, cons_max, nvalid, oms, nL_arr)."""
    kK, kbar, cH, xH = geom(w); T0g = kK/(2*np.pi); oms = XR*T0g
    nL, cons = [], []
    for om in oms:
        r = occ_bench(w, om, L, N) if bench else occ_new(w, om, L, N, lam)
        if r is None:
            nL.append(np.nan); cons.append(np.nan)
        else:
            nL.append(r[0]); cons.append(r[2])
    nL = np.array(nL); cons = np.array(cons)
    m = ~np.isnan(nL)
    if m.sum() < 4:
        return (np.nan,)*3 + (np.nanmax(cons), int(m.sum()), oms, nL)
    Tw = oms[m]/np.log(1+1/nL[m])
    Tq = np.polyfit(oms[m], Tw, 2)[-1]
    Tl = np.polyfit(oms[m], Tw, 1)[-1]
    Tc = np.polyfit(oms[m], Tw, 3)[-1]
    return Tq, Tl, Tc, np.nanmax(cons), int(m.sum()), oms, nL

# =====================================================================
# FASE A  — PUERTA DURA (si falla algo, se aborta sin tocar la Fase B)
# =====================================================================
print("="*70); print("RESIDUAL ANATOMY — PASO 0  (gate de viabilidad)"); print("="*70)
w = 2.0; L, N = box(w); kK, kbar, cH, xH = geom(w); T0g = kK/(2*np.pi)
print(f"\nw={w}  L={L}  N={N}  kappa_local={kK:.10f}  kloc/2pi={T0g:.10f}")
print(f"kappa_bar/2pi={kbar/(2*np.pi):.10f}  cH={cH:.6f}  xH={xH:.6f}")

A_ok = True

# --- A1: ruta nueva (scatter_field,lam=1) vs ruta benchmark (solve_inject) punto a punto ---
print("\n[A1] Ruta NUEVA (scatter_field, lam=1, op=validated) vs ruta BENCHMARK (solve_inject)")
print("     mismo operador -> deben coincidir a precision de maquina")
oms = XR*T0g; maxrel = 0.0
print(f"     {'omega':>9} {'nL_new':>14} {'nL_bench':>14} {'|rel.dif|':>10}")
for om in oms:
    rn = occ_new(w, om, L, N, 1.0); rb = occ_bench(w, om, L, N)
    rel = abs(rn[0]-rb[0])/abs(rb[0]); maxrel = max(maxrel, rel)
    print(f"     {om:9.5f} {rn[0]:14.8e} {rb[0]:14.8e} {rel:10.2e}")
print(f"     -> max|rel.dif| = {maxrel:.2e}   {'OK' if maxrel<1e-6 else 'FALLO (>1e-6)'}")
A_ok &= maxrel < 1e-6

# --- A2: regresion de T_upstream/(kloc/2pi) en lam=1 contra el valor congelado ---
print("\n[A2] Regresion: Tup/(kloc/2pi) [ruta nueva, lam=1] vs benchmark congelado")
Tq, Tl, Tc, cons, nv, _, _ = T0_orders(w, L, N, 1.0)
ratio_new = Tq/T0g
print(f"     Tup/kloc2pi (ruta nueva) = {ratio_new:.10f}")
print(f"     referencia congelada     = {REF_TUP_OVER_KLOC:.10f}")
print(f"     |dif| = {abs(ratio_new-REF_TUP_OVER_KLOC):.2e}   tol={TOL_REGRESION:.0e}   "
      f"{'OK' if abs(ratio_new-REF_TUP_OVER_KLOC)<TOL_REGRESION else 'FALLO'}")
A_ok &= abs(ratio_new-REF_TUP_OVER_KLOC) < TOL_REGRESION

# --- A3: fondo y operador invariantes / consistentes al variar lam ---
print("\n[A3] Invarianza del fondo y consistencia del operador lam")
x = np.linspace(-L, L, N)
n0, v0, V0 = C.bg(x, w, 'horizon')
# el fondo no depende de lam (bg no recibe lam); verificamos operador lam=1 == build_MBdG
with contextlib.redirect_stdout(open(os.devnull, 'w')):
    from fd_bdg_operator import build_MBdG
M1 = C.build_op_lambda(x, n0, v0, V0, 1.0)
Mb = build_MBdG(x, n0, v0, V0, C.mu)
opdif = abs((M1-Mb)).max()
# kappa_local no depende de lam (geom no usa lam): control explicito
kK_b = geom(2.0)[0]
print(f"     max|build_op_lambda(1) - build_MBdG| = {opdif:.2e}   "
      f"{'OK' if opdif<1e-10 else 'FALLO'}")
print(f"     kappa_local es geometrico (geom no usa lam) -> invariante por construccion: {kK_b:.10f}")
A_ok &= opdif < 1e-10

print("\n" + "-"*70)
if not A_ok:
    print("FASE A: FALLO -> se ABORTA antes del barrido. Revisar la ruta nueva.")
    sys.exit(1)
print("FASE A: OK -> la ruta nueva reproduce el benchmark. Continuo a Fase B.")

# =====================================================================
# FASE B — barrido dispersivo en lam (w=2)
# =====================================================================
print("\n" + "="*70); print("FASE B — barrido dispersivo  w=2,  lam in {0.5, 1, 2, 4}"); print("="*70)
LAMS = [0.5, 1.0, 2.0, 4.0]
res = {}
print(f"\n{'lam':>5} {'Tup_quad':>11} {'R_quad':>10} {'R_lin':>10} {'R_cub':>10} "
      f"{'cons_max':>9} {'nvalid':>6}")
for lam in LAMS:
    Tq, Tl, Tc, cons, nv, oms, nLa = T0_orders(w, L, N, lam)
    Rq = Tq/T0g-1; Rl = Tl/T0g-1; Rc = Tc/T0g-1
    res[lam] = dict(Rq=Rq, Rl=Rl, Rc=Rc, cons=cons, nv=nv, Tq=Tq,
                    wmax_ok=(nv == len(XR)))
    print(f"{lam:5.2f} {Tq:11.7f} {Rq:10.4%} {Rl:10.4%} {Rc:10.4%} {cons:9.1e} {nv:6d}")

# --- dR PAREADO por orden (control #2 de ChatGPT) ---
print("\n[dR PAREADO]  dR(lam)=R(lam)-R(1) por orden de extrapolacion (misma banda/puntos)")
R1 = res[1.0]
print(f"{'lam':>5} {'dR_quad':>11} {'dR_lin':>11} {'dR_cub':>11} {'spread_dR':>11}")
dR_table = {}
for lam in LAMS:
    if lam == 1.0: continue
    dRq = res[lam]['Rq']-R1['Rq']; dRl = res[lam]['Rl']-R1['Rl']; dRc = res[lam]['Rc']-R1['Rc']
    spread = max(abs(dRq-dRl), abs(dRq-dRc))      # spread de dR entre ordenes = presupuesto pareado
    dR_table[lam] = dict(dRq=dRq, dRl=dRl, dRc=dRc, spread=spread)
    print(f"{lam:5.2f} {dRq:11.4%} {dRl:11.4%} {dRc:11.4%} {spread:11.2e}")

# --- convergencia pareada en N (control #2): repetir dR con N mas fino ---
print("\n[CONVERGENCIA PAREADA]  dR_quad con N y 1.5N (mismo L); debe mantenerse")
N2 = int(round(N*1.5))
print(f"{'lam':>5} {'dR_quad(N)':>12} {'dR_quad(1.5N)':>14} {'|cambio|':>10}")
for lam in LAMS:
    if lam == 1.0: continue
    Tq_a = T0_orders(w, L, N, lam)[0]; Tq1_a = T0_orders(w, L, N, 1.0)[0]
    Tq_b = T0_orders(w, L, N2, lam)[0]; Tq1_b = T0_orders(w, L, N2, 1.0)[0]
    dR_a = (Tq_a-Tq1_a)/T0g; dR_b = (Tq_b-Tq1_b)/T0g
    print(f"{lam:5.2f} {dR_a:12.4%} {dR_b:14.4%} {abs(dR_a-dR_b):10.2e}")

# =====================================================================
# VEREDICTO
# =====================================================================
print("\n" + "="*70); print("VEREDICTO"); print("="*70)
# presupuesto dispersivo = max spread pareado entre ordenes sobre los lam!=1
budget = max(d['spread'] for d in dR_table.values())
signal = max(abs(d['dRq']) for d in dR_table.values())
print(f"(i)  Residuo TOTAL: R(1)={R1['Rq']:.4%}  (congelado -4.28%, presupuesto ~0.06%)")
print(f"     -> ya detectado por el benchmark (factor ~70 sobre el ruido).")
print(f"(ii) Dependencia DISPERSIVA: max|dR(lam)|={signal:.4%}  vs presupuesto pareado={budget:.2e}")
ratio = signal/budget if budget > 0 else np.inf
if ratio > 3:
    print(f"     senal/presupuesto = {ratio:.1f}  (>3) -> DEPENDENCIA DISPERSIVA DETECTADA")
    print(f"     RECOMENDACION: continuar al Paso 1 (barrido lam fino + eje w).")
else:
    print(f"     senal/presupuesto = {ratio:.1f}  (<=3) -> dependencia NO distinguible")
    print(f"     del ruido en este rango de lam. (NO implica 'dispersion sin efecto').")
# ventana
allwin = all(res[l]['wmax_ok'] for l in LAMS)
print(f"\nVentana espectral: los 9 puntos dentro de la ventana propagante en todo lam: "
      f"{'SI' if allwin else 'NO (revisar)'}")

# --- CSV ---
with open(os.path.join(HERE, 'out', 'residual_gate.csv'), 'w', newline='') as fh:
    wr = csv.writer(fh)
    wr.writerow(['lam', 'Tup_quad', 'R_quad', 'R_lin', 'R_cub', 'cons_max', 'nvalid',
                 'dR_quad', 'dR_lin', 'dR_cub', 'spread_dR_pareado'])
    for lam in LAMS:
        d = res[lam]; dr = dR_table.get(lam, dict(dRq=0, dRl=0, dRc=0, spread=0))
        wr.writerow([lam, d['Tq'], d['Rq'], d['Rl'], d['Rc'], d['cons'], d['nv'],
                     dr['dRq'], dr['dRl'], dr['dRc'], dr['spread']])
print(f"\nCSV: residual_anatomy/out/residual_gate.csv")
