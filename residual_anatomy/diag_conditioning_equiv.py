"""
DIAGNOSTICO (no gate) — equivalencia lam <-> w/sqrt(lam) con AMBOS lados resueltos con refinamiento iterativo
en precision mixta (residuo en long double), para separar la identidad del ruido de redondeo del solve.
Rama diag-condicionamiento (2026-11-17). Sin modificar correlations.py: se sustituye en tiempo de ejecucion
C.spla por un envoltorio cuyo spsolve es 'eqref' (equilibrado de Ruiz + refinamiento) de diag_conditioning.

Lado A: operador lam con fondo consistente (Qp_lam = sqrt(lam) Qp; envoltorio identico al de diag_T2), w=2, L=40.
Lado B: v1.0 (op='validated', lam=1) en w'=2/sqrt(lam), L'=40/sqrt(lam), mismo N, omega'=sqrt(lam) omega.
Se recalculan ambos con 'base' (spsolve) y con 'eqref' en los 90 puntos de out/diag_T2_nL.csv (QP_LAM).
Salida: out/diag_conditioning_equiv.csv
"""
import os, sys, io, csv, time, types, contextlib
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver')); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C
    import diag_conditioning as D

ORIG = C.build_op_lambda
def build_qp_lam(x, n, v, V, lam=1.0):                     # = diag_T2_qp_consistent.build_qp_lam
    M = ORIG(x, n, v, V, lam)
    Qp = V - C.mu + C.g*n + 0.5*v**2
    d = (np.sqrt(lam)-1.0)*Qp
    return (M + sp.bmat([[sp.diags(d), None], [None, sp.diags(-d)]], format='csr')).tocsr()

SPLA_BASE = C.spla
SPLA_REF = types.SimpleNamespace(spsolve=lambda A, b: D.solve(A.tocsc(), b, 'eqref'))


def nL_side(side, lam, N, om, kind):
    C.spla = SPLA_REF if kind == 'eqref' else SPLA_BASE
    try:
        if side == 'A':
            C.build_op_lambda = build_qp_lam
            op = 'validated' if lam == 1.0 else 'lambda'
            r = C.scatter_field(40.0, N, om, 'd2', w=2.0, lam=lam, op=op)
        else:
            s = np.sqrt(lam)
            r = C.scatter_field(40.0/s, N, s*om, 'd2', w=2.0/s, lam=1.0, op='validated')
    finally:
        C.build_op_lambda = ORIG; C.spla = SPLA_BASE
    return abs(r['out_amp'].get(('L', 1), 0.0))**2


t0 = time.time()
ref = [q for q in csv.DictReader(open(os.path.join(OUT, 'diag_T2_nL.csv'))) if q['variante'] == 'QP_LAM']
rows = []
for q in ref:
    lam, N, om = float(q['lam']), int(q['N']), float(q['omega'])
    row = dict(lam=lam, N=N, omega=om, nL_T2_csv=float(q['n_L']))
    for kind in ('base', 'eqref'):
        row[f'A_{kind}'] = nL_side('A', lam, N, om, kind); row[f'B_{kind}'] = nL_side('B', lam, N, om, kind)
    row['rel_AB_base'] = abs(row['A_base']-row['B_base'])/row['B_base']
    row['rel_AB_eqref'] = abs(row['A_eqref']-row['B_eqref'])/row['B_eqref']
    row['rel_A_base_vs_csv'] = abs(row['A_base']-row['nL_T2_csv'])/row['nL_T2_csv']
    rows.append(row)
by = {}
for r in rows:
    by.setdefault((r['lam'], r['N']), []).append(r)
print("max dif. relativa n_L, lado A (lam, fondo consistente) vs lado B (v1.0 reescalado):")
for (lam, N), rr in sorted(by.items()):
    print(f"  lam={lam:4.2f} N={N:5d}:  base {max(r['rel_AB_base'] for r in rr):.1e}   eqref {max(r['rel_AB_eqref'] for r in rr):.1e}"
          f"   (A_base vs CSV de T2: {max(r['rel_A_base_vs_csv'] for r in rr):.1e})")
print(f"GLOBAL: base {max(r['rel_AB_base'] for r in rows):.1e}   eqref {max(r['rel_AB_eqref'] for r in rows):.1e}   [{time.time()-t0:.0f} s]")
with open(os.path.join(OUT, 'diag_conditioning_equiv.csv'), 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
