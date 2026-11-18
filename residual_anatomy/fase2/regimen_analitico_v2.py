"""
FASE 2 — sensibilidad de la estimacion analitica del regimen (addendum 3). Igual que regimen_analitico.py, pero con el
termino -p^5/(128 c^3) de la expansion de Omega que faltaba (lo senalo el revisor externo):
  Omega ~ c p + p^3/(8c) + a p^5,   a = gamma c/2 - 1/(128 c^3)
  fase  Phi(p) = p^3/(24 c kappa) + a p^5/(5 kappa);  R = 4.8 a c p*^2;  d ~ [p*^2/(8c) + a p*^4]/kappa
Se toma la PRIMERA raiz de Phi(p) = umbral al crecer p (con a < 0, Phi tiene maximo; si no llega al umbral: 'sin p*').
"""
import csv, os, numpy as np
from scipy.optimize import brentq
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rows = list(csv.DictReader(open(os.path.join(ROOT, 'out', 'fp_benchmark.csv'))))
K = np.array([float(r['kappa_local']) for r in rows]); cH = float(rows[0]['cH'])
def first_root(f, hi=50.0):
    ps = np.linspace(1e-4, hi, 200000); v = f(ps); i = np.where(np.sign(v[:-1]) != np.sign(v[1:]))[0]
    return brentq(f, ps[i[0]], ps[i[0]+1]) if len(i) else None
for variante, extra in (('sin el termino (addendum 2)', 0.0), ('con -p^5/(128c^3)', -1/(128*cH**3))):
    print(f"=== {variante} ===")
    for thr in (0.5, 1.0, 2.0):
        for g in (0.0, 1.0, 10.0):
            a = g*cH/2 + extra
            ps = [first_root(lambda p, k=k: p**3/(24*cH*k) + a*p**5/(5*k) - thr) for k in K]
            if any(p is None for p in ps):
                print(f"  umbral {thr} gamma={g:4}: sin p* en alguna anchura"); continue
            Rs = [4.8*a*cH*p**2 for p in ps]; ds = [(p**2/(8*cH) + a*p**4)/k for p, k in zip(ps, K)]
            peff = -np.polyfit(np.log(K), np.log(ds), 1)[0]
            reg = ['k6' if R > 10 else ('k4' if abs(R) < 0.1 else 'transicion') for R in Rs]
            print(f"  umbral {thr} gamma={g:4}: R = " + ', '.join(f'{R:+.2f}' for R in Rs) + f" | {reg} | p_eff = {peff:.3f}")
