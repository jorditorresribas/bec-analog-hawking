"""
Espectro de Hawking para la web (panel «espectro cuántico»): n_L(ω) del solver v1.0 del repositorio bec-analog-hawking
(tag v2.0), SIN modificarlo: mode_matching.solve_inject, perfil tanh w = 2 (Mach 0,5 → 2), L = 40, N = 9600 (la malla del
benchmark). n_L = ocupación del canal Hawking aguas arriba (Finazzi-Parentani Ec. 42), el mismo observable del benchmark.
Frecuencias: 40 valores log-espaciados entre 0,05·T_H y 0,97·ω_max. T_H = κ_local/2π.
Salida: hawking_spectrum.json (lo lee la web) y hawking_spectrum.csv. Uso: python3 hawking_spectrum.py /ruta/al/repo
"""
import sys, os, io, json, csv, hashlib, subprocess, contextlib
import numpy as np
REPO = os.path.abspath(sys.argv[1] if len(sys.argv) > 1 else '/Users/mba15/bec-analog-hawking')
sys.path.insert(0, os.path.join(REPO, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import mode_matching as MM
W, L, N = 2.0, 40.0, 9600
x = np.linspace(-L, L, 80000); MM.w = W
n, v, V = MM.bg(x, 'horizon'); c = np.sqrt(n); f = c - v
i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]; xH = x[i] + f[i]/(f[i]-f[i+1])*(x[i+1]-x[i])
kappa = abs(np.interp(xH, x, np.gradient(f, x))); TH = kappa/(2*np.pi)
def nL(om):
    MM.omega = om; xx = np.linspace(-L, L, N); nn, vv, VV = MM.bg(xx, 'horizon')
    Rin = [q for q in MM.modes_full(nn[-1], vv[-1], VV[-1]) if q['prop'] and q['vg'] < 0]
    if not any(q['norm'] < 0 for q in Rin): return None
    with contextlib.redirect_stdout(io.StringIO()):
        out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=[1.0 if q['norm'] < 0 else 0.0 for q in Rin])
    return sum(a for (q, s, a) in out if q > 0 and s == 'L')
lo, hi = 0.1, 1.0                                    # omega_max: ultima omega con canal de norma negativa
for _ in range(40):
    m = 0.5*(lo+hi); lo, hi = (m, hi) if nL(m) is not None else (lo, m)
omax = lo
oms = np.geomspace(0.05*TH, 0.97*omax, 40)
rows = [dict(omega=float(o), omega_over_TH=float(o/TH), n_L=float(nL(o))) for o in oms]
for r in rows:
    r['n_planck_TH'] = 1/np.expm1(r['omega_over_TH']); r['T_omega_over_TH'] = r['omega_over_TH']/np.log1p(1/r['n_L'])
sha = {p: hashlib.sha256(open(os.path.join(REPO, p), 'rb').read()).hexdigest() for p in ('solver/mode_matching.py', 'solver/fd_bdg_operator.py')}
head = subprocess.run(['git', '-C', REPO, 'rev-parse', 'HEAD'], capture_output=True, text=True).stdout.strip()
tag = subprocess.run(['git', '-C', REPO, 'describe', '--tags', '--exact-match'], capture_output=True, text=True).stdout.strip()
meta = dict(w=W, L=L, N=N, kappa_local=float(kappa), T_H=float(TH), omega_max=float(omax), repo_commit=head, repo_tag=tag, sha256=sha)
here = os.path.dirname(os.path.abspath(__file__))
json.dump(dict(meta=meta, data=rows), open(os.path.join(here, 'hawking_spectrum.json'), 'w'), indent=1)
with open(os.path.join(here, 'hawking_spectrum.csv'), 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=list(rows[0])); wr.writeheader(); wr.writerows(rows)
print(json.dumps(meta, indent=1)); print("omega/T_H, n_L, Planck(T_H), T_omega/T_H:")
for r in rows[::4]: print(f"  {r['omega_over_TH']:7.3f}  {r['n_L']:.4e}  {r['n_planck_TH']:.4e}  {r['T_omega_over_TH']:.4f}")
