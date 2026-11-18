"""
DIAGNOSTICO (no gate) — separar discretizacion, redondeo y extrapolacion en omega.
Rama diag-condicionamiento (2026-11-17). No modifica solver ni gates.

(E1) coeficiente a de Larre et al. 2012 (m_u = 0.5). Fuentes de error que se separan:
  - redondeo del solve: 'base' (spsolve) frente a 'eqref' (equilibrado + refinamiento, diag_conditioning)
  - discretizacion en dx: dx = 0.02, 0.01, 0.005, 0.0025 con Richardson (p observado y p = 2)
  - extrapolacion omega -> 0 (ajuste polinomico de omega*n_L): ventanas cada vez mas bajas y ordenes 2-4
(FP) T_up/(kappa_local/2pi) del benchmark para w = 2 y 8: N x (1, 1.5, 2.25) con Richardson, base frente a eqref.
Salida: out/diag_convergence_e1.csv, out/diag_convergence_fp.csv
"""
import os, sys, io, csv, time, types, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver')); sys.path.insert(0, HERE); sys.path.insert(0, os.path.join(HERE, 'gates'))
with contextlib.redirect_stdout(io.StringIO()):
    import diag_conditioning as D
    import gate_E1_lrcp as E
    import mode_matching as MM
t0 = time.time()


def e1_nL(m, dx, om, kind, L=30.0):
    x, n, v, V, mu = E.bg_lrcp(E.grid(L, dx), L, m)
    A, rhs, nL_of = D.assemble(x, n, v, V, mu, om)
    return nL_of(D.solve(A, rhs, kind))


# ---------------- E1 ----------------
m = 0.5; a_th, b_th, T57 = E.lrcp_theory(m)
DX4 = (0.02, 0.01, 0.005, 0.0025)
WINS = {'W1 [0.1,0.6]': (0.1, 0.6), 'W2 [0.05,0.3]': (0.05, 0.3), 'W3 [0.025,0.15]': (0.025, 0.15), 'W4 [0.0125,0.075]': (0.0125, 0.075)}
rows_e1 = []; nLs = {}
for wname, (lo, hi) in WINS.items():
    oms = T57*np.geomspace(lo, hi, 12)
    for dx in DX4:
        for kind in (('base', 'eqref') if wname.startswith('W1') else ('eqref',)):
            nLs[(wname, dx, kind)] = np.array([e1_nL(m, dx, om, kind) for om in oms])
    print(f"  E1 {wname} hecho [{time.time()-t0:.0f} s]", flush=True)
    for dx in DX4:
        for kind in (('base', 'eqref') if wname.startswith('W1') else ('eqref',)):
            for order in (2, 3, 4):
                a, b = E.fit_ab(oms, nLs[(wname, dx, kind)], order)
                rows_e1.append(dict(ventana=wname, dx=dx, solver=kind, orden=order, a=a, a_rel=(a-a_th)/a_th))

def a_of(w, dx, kind, order):
    return next(r['a'] for r in rows_e1 if r['ventana'] == w and r['dx'] == dx and r['solver'] == kind and r['orden'] == order)

print(f"\nE1 a_th = {a_th:.10f}")
print("(1) redondeo: max |a(base) - a(eqref)|/a_th en W1, por dx y orden:")
for dx in DX4:
    print(f"   dx={dx}: " + '  '.join(f"o{o}: {abs(a_of('W1 [0.1,0.6]', dx, 'base', o)-a_of('W1 [0.1,0.6]', dx, 'eqref', o))/a_th:.1e}" for o in (2, 3, 4)))
print("(2) discretizacion: a(dx) con eqref, orden 3, y Richardson con los 3 dx mas finos (p observado):")
rich = {}
for w in WINS:
    seq = [a_of(w, dx, 'eqref', 3) for dx in DX4]
    p = np.log2(abs(seq[1]-seq[2])/abs(seq[2]-seq[3]))
    r_ = seq[3] + (seq[3]-seq[2])/(2**p-1); rich[w] = r_
    print(f"   {w:18s}: " + ', '.join(f'{(s-a_th)/a_th:+.2e}' for s in seq) + f"  (err. rel. por dx)  p={p:.2f}  Richardson err {(r_-a_th)/a_th:+.2e}")
print("(3) extrapolacion en omega: Richardson-dx por ventana y orden (err. rel. frente a a_th):")
for w in WINS:
    out = []
    for o in (2, 3, 4):
        seq = [a_of(w, dx, 'eqref', o) for dx in DX4]
        r_ = seq[3] + (seq[3]-seq[2])/3.0
        out.append(f"o{o}: {(r_-a_th)/a_th:+.2e}")
    print(f"   {w:18s}: " + '  '.join(out))

# ---------------- FP ----------------
n_u, md, jc = 1.0, 2.0, 0.5; n_d = (jc/md)**(2/3)
def kloc(w):
    L = max(40., 8*w+20.); x = np.linspace(-L, L, 80000); n = n_u+(n_d-n_u)*0.5*(1+np.tanh(x/w))
    f = np.sqrt(n)-jc/n; i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    xH = x[i]+f[i]/(f[i]-f[i+1])*(x[i+1]-x[i]); return abs(np.interp(xH, x, np.gradient(f, x)))
SPLA_BASE = MM.spla
SPLA_REF = types.SimpleNamespace(spsolve=lambda A, b: D.solve(A.tocsc(), b, 'eqref'))
def T_ratio(w, N, kind):
    L = max(40., 8*w+20.); T0 = kloc(w)/(2*np.pi)
    oms = np.array([0.05, 0.065, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.27])*T0; nL = []
    MM.spla = SPLA_REF if kind == 'eqref' else SPLA_BASE
    try:
        for om in oms:
            MM.w = w; MM.omega = om; x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
            Rin = [q for q in MM.modes_full(n[-1], v[-1], V[-1]) if q['prop'] and q['vg'] < 0]
            out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=[1.0 if q['norm'] < 0 else 0.0 for q in Rin])
            nL.append(sum(a for (q, s, a) in out if q > 0 and s == 'L'))
    finally:
        MM.spla = SPLA_BASE
    Tw = oms/np.log(1+1/np.array(nL)); return np.polyfit(oms, Tw, 2)[-1]/T0
rows_fp = []
print("\nFP T_up/(kappa_local/2pi): N0 = L*240 (benchmark), factores 1, 1.5, 2.25")
for w in (2.0, 8.0):
    L = max(40., 8*w+20.); N0 = int(round(L*240))
    res = {}
    for f in (1.0, 1.5, 2.25):
        N = int(round(N0*f))
        for kind in ('base', 'eqref'):
            res[(f, kind)] = T_ratio(w, N, kind)
            rows_fp.append(dict(w=w, N=N, factor=f, solver=kind, T_over_kloc2pi=res[(f, kind)]))
    seq = [res[(f, 'eqref')] for f in (1.0, 1.5, 2.25)]
    p = np.log(abs(seq[0]-seq[1])/abs(seq[1]-seq[2]))/np.log(1.5)
    r_ = seq[2] + (seq[2]-seq[1])/(1.5**p-1)
    print(f"   w={w}: eqref " + ', '.join(f'{s:.8f}' for s in seq) + f"  p={p:.2f}  Richardson {r_:.8f}"
          f"  | redondeo |base-eqref| max {max(abs(res[(f,'base')]-res[(f,'eqref')]) for f in (1.0,1.5,2.25)):.1e}"
          f"  | discretizacion en N0 {abs(seq[0]-r_):.1e}   [{time.time()-t0:.0f} s]", flush=True)

for name, rows in (('diag_convergence_e1.csv', rows_e1), ('diag_convergence_fp.csv', rows_fp)):
    with open(os.path.join(OUT, name), 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
print(f"fin ({time.time()-t0:.0f} s)")
