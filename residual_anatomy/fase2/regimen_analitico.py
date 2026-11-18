"""
FASE 2 — estimacion ANALITICA del regimen (sin resolver BdG; no es la medicion). Segundo addendum al pre-registro.
Argumento de la fase (docs/FASE2_DISENO.md §3) con la dispersion MIXTA: Omega ~ c p + p^3/(8c) + gamma c p^5/2, asi que la
fase dispersiva es Phi(p) = p^3/(24 c kappa) + gamma c p^5/(10 kappa).
  p*: Phi(p*) = umbral (de orden 1; se prueban 0.5, 1 y 2)
  R = (termino k^6)/(termino k^4) en p* = 2.4 gamma c^2 p*^2  ->  k4 si R < 0.1, k6 si R > 10, transicion en medio
  d ~ [p*^2/(8c) + gamma c p*^4/2]/kappa  ->  p_eff = -pendiente de log d frente a log kappa en las cinco anchuras
c = c_H; kappa = kappa_local de out/fp_benchmark.csv (w = 2, 3, 4, 6, 8).
"""
import csv, os, numpy as np
from scipy.optimize import brentq
ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
rows = list(csv.DictReader(open(os.path.join(ROOT, 'out', 'fp_benchmark.csv'))))
W = [float(r['w']) for r in rows]; K = np.array([float(r['kappa_local']) for r in rows]); cH = float(rows[0]['cH'])
print(f"c_H = {cH:.4f}; kappa_local: " + ', '.join(f"w={w:g}: {k:.4f}" for w, k in zip(W, K)))
for thr in (0.5, 1.0, 2.0):
    print(f"--- umbral de fase {thr} ---")
    for g in (0.0, 1.0, 10.0):
        ps = [brentq(lambda p: p**3/(24*cH*k) + g*cH*p**5/(10*k) - thr, 1e-6, 1e3) for k in K]
        Rs = [2.4*g*cH**2*p**2 for p in ps]
        ds = [(p**2/(8*cH) + g*cH*p**4/2)/k for p, k in zip(ps, K)]
        peff = -np.polyfit(np.log(K), np.log(ds), 1)[0]
        reg = ['k6' if R > 10 else ('k4' if R < 0.1 else 'transicion') for R in Rs]
        print(f"  gamma={g:4}: R = " + ', '.join(f'{R:.2f}' for R in Rs) + f" | {reg} | p_eff = {peff:.3f}")
