"""
DIAGNOSTICO T2 (no gate) — ¿desaparece el colapso de n_L(lam) con el fondo consistente?

Mismas condiciones que el PASO 0 (residual_gate.py): w=2, (L,N)=box(2), banda omega = XR*kappa_local/2pi
(9 puntos), observable n_L (canal u upstream, norma +, inyectando d2), estimador T_omega=omega/ln(1+1/n_L),
T0 por extrapolacion polinomica (cuad = referencia; lin y cub para el spread), ruta correlations.scatter_field.

Dos variantes PAREADAS (mismo solver, misma malla, misma banda):
  QP_FIJO : build_op_lambda tal cual (Qp lam-indep)  -> reproduce el Paso 0 (colapso)
  QP_LAM  : build_op_lambda + diag((sqrt(lam)-1)Qp) en L+ y -(...) en -L-  -> BdG consistente con
            masa m'=1/sqrt(lam), g'n=c^2/sqrt(lam) (NOVELTY_AUDIT Apendice A; Goldstone limpio, diag_goldstone).
NO modifica correlations.py: sustituye C.build_op_lambda en tiempo de ejecucion por un envoltorio.
Integridad previa: (I1) el envoltorio en lam=1 es identico al original; (I2) regresion Tup/(kloc/2pi)
en lam=1 contra el benchmark congelado 0.9572413965 (tol 1e-4, como el Paso 0).
Salidas: ../out/diag_T2.csv, ../out/diag_T2_nL.csv, ../out/diag_T2_manifest.json
"""
import os, sys, io, csv, json, time, hashlib, platform, contextlib
import numpy as np, scipy, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C

REF_TUP_OVER_KLOC = 0.9572413965031795; TOL_REG = 1e-4
LAMS = [0.25, 0.5, 1.0, 2.0, 4.0]          # Paso 0 = {0.5,1,2,4}; 0.25 = ampliacion (rango de diag_goldstone)
XR = np.array([0.05, 0.065, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.27])

# ---- geom()/box() copiados verbatim de residual_gate.py (a su vez de FP_benchmark_table.py) ----
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

ORIG = C.build_op_lambda
def build_qp_lam(x, n, v, V, lam=1.0):
    M = ORIG(x, n, v, V, lam)
    Qp = V - C.mu + C.g*n + 0.5*v**2
    d = (np.sqrt(lam)-1.0)*Qp
    return (M + sp.bmat([[sp.diags(d), None], [None, sp.diags(-d)]], format='csr')).tocsr()
VARIANTS = {'QP_FIJO': ORIG, 'QP_LAM': build_qp_lam}

def occ(w, om, L, N, lam, variant):
    C.build_op_lambda = VARIANTS[variant]
    try:
        op = 'validated' if lam == 1.0 else 'lambda'
        r = C.scatter_field(L, N, om, 'd2', w=w, lam=lam, op=op)
    finally:
        C.build_op_lambda = ORIG
    if r['Rin_neg'] is None or r['nchan'] != 2:
        return None
    oa = r['out_amp']; nL = abs(oa.get(('L', 1), 0.0))**2; nRp = abs(oa.get(('R', 1), 0.0))**2
    return nL, nL+nRp, abs(r['flux_out'] - r['eta_in'])

def T0_orders(w, L, N, lam, variant, oms):
    nL, cons = [], []
    for om in oms:
        r = occ(w, om, L, N, lam, variant)
        nL.append(np.nan if r is None else r[0]); cons.append(np.nan if r is None else r[2])
    nL = np.array(nL); cons = np.array(cons); m = ~np.isnan(nL)
    if m.sum() < 4:
        return dict(Tq=np.nan, Tl=np.nan, Tc=np.nan, cons=np.nanmax(cons) if m.any() else np.nan, nv=int(m.sum()), nL=nL)
    Tw = oms[m]/np.log(1+1/nL[m])
    return dict(Tq=np.polyfit(oms[m], Tw, 2)[-1], Tl=np.polyfit(oms[m], Tw, 1)[-1],
                Tc=np.polyfit(oms[m], Tw, 3)[-1], cons=float(np.nanmax(cons)), nv=int(m.sum()), nL=nL)

t0 = time.time(); w = 2.0; L, N = box(w); kK, kbar, cH, xH = geom(w); T0g = kK/(2*np.pi); oms = XR*T0g
print(f"T2  w={w} L={L} N={N} kloc/2pi={T0g:.10f}")
# ---- integridad ----
x = np.linspace(-L, L, N)
with contextlib.redirect_stdout(io.StringIO()):
    n0, v0, V0 = C.bg(x, w, 'horizon')
I1 = float(abs(build_qp_lam(x, n0, v0, V0, 1.0) - ORIG(x, n0, v0, V0, 1.0)).max())
reg = T0_orders(w, L, N, 1.0, 'QP_LAM', oms)
I2 = abs(reg['Tq']/T0g - REF_TUP_OVER_KLOC)
ok = (I1 == 0.0) and (I2 < TOL_REG)
print(f"I1 max|envoltorio(1)-original(1)| = {I1:.1e}   I2 |Tup/kloc2pi - ref| = {I2:.2e}  -> {'OK' if ok else 'STOP'}")
rows, rows_nL = [], []
if ok:
    for Nn in (N, int(round(1.5*N))):
        for variant in ('QP_FIJO', 'QP_LAM'):
            for lam in LAMS:
                r = T0_orders(w, L, Nn, lam, variant, oms)
                rows.append(dict(variante=variant, lam=lam, N=Nn, L=L, Tq=r['Tq'], Tl=r['Tl'], Tc=r['Tc'],
                                 Teff_over_TH=r['Tq']/T0g, spread_ordenes=(max(abs(r['Tq']-r['Tl']), abs(r['Tq']-r['Tc']))/T0g
                                                                           if np.isfinite(r['Tq']) else np.nan),
                                 cons_max=r['cons'], n_validas=r['nv']))
                for om, nl in zip(oms, r['nL']):
                    rows_nL.append(dict(variante=variant, lam=lam, N=Nn, omega=om, omega_over_T0g=om/T0g, n_L=nl))
                print(f"  N={Nn:5d} {variant:8s} lam={lam:4.2f}  Teff/TH={r['Tq']/T0g:9.5f}  "
                      f"(lin {r['Tl']/T0g:8.5f}, cub {r['Tc']/T0g:8.5f})  cons={r['cons']:.1e}  nval={r['nv']}")
secs = time.time()-t0
def wcsv(p, rr):
    with open(p, 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rr[0].keys())); wr.writeheader(); wr.writerows(rr)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
man = dict(tipo='DIAGNOSTICO T2 (no gate)', script='residual_anatomy/diag_T2_qp_consistent.py',
           timestamp=time.strftime('%Y-%m-%d %H:%M:%S %Z'), segundos=secs, integridad='PASS' if ok else 'STOP',
           I1_max_abs=I1, I2_regresion=I2, ref=REF_TUP_OVER_KLOC, w=w, L=L, N=[N, int(round(1.5*N))],
           lams=LAMS, XR=XR.tolist(), kloc_over_2pi=T0g,
           dependencias_sha256={p: sha(os.path.join(ROOT, p)) for p in
                                ('solver/correlations.py', 'solver/fd_bdg_operator.py', 'residual_anatomy/residual_gate.py',
                                 'residual_anatomy/diag_T2_qp_consistent.py')},
           entorno=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                        platform=platform.platform()))
if rows:
    wcsv(os.path.join(OUT, 'diag_T2.csv'), rows); wcsv(os.path.join(OUT, 'diag_T2_nL.csv'), rows_nL)
    man['csv'] = {p: sha(os.path.join(OUT, p)) for p in ('diag_T2.csv', 'diag_T2_nL.csv')}
json.dump(man, open(os.path.join(OUT, 'diag_T2_manifest.json'), 'w'), indent=1)
print(f"fin ({secs:.0f} s)  integridad={'PASS' if ok else 'STOP'}")
sys.exit(0 if ok else 2)
