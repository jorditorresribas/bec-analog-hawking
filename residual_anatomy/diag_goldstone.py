"""
DIAGNOSTICO (no gate) — test del modo de Goldstone sobre los operadores del proyecto.

Pregunta: el vector de fase del condensado Phi_G = (psi0, -psi0*) = (sqrt(n), -sqrt(n))
(convencion del proyecto: fase del condensado ya extraida, L± = -½(d ± iv)² + W)
¿es modo de frecuencia cero de cada operador?  Si el fondo es solucion estacionaria de la
GPE cuyo linealizado es M, entonces M·Phi_G = 0 (salvo discretizacion).

Operadores (NINGUNO se modifica; solo se importan/copian):
  A. v1.0 build_MBdG (fd_bdg_operator.py)                         lam = 1
  B. build_op_lambda (correlations.py), Qp lam-indep (como se uso)  lam en LAMS
  B'. idem con Qp_lam = sqrt(lam)·Qp (BdG consistente con masa m'=1/sqrt(lam))
  C. no local Weyl (residual_anatomy/gate0_pseudoherm.build_nonlocal_M), Qp lam-indep
  D. Pade 2o orden (Kpade2 copiado verbatim de gate_pade_audit2.py) + conveccion skew de gate0
Prediccion (paper):
  B : residuo = ½(1-sqrt(lam))·(sqrt n)''  (en ambas componentes)
  B': residuo = O(dx²) (independiente de lam)
  C/D: residuo = (K_lam - K_1)·sqrt(n)  (≠0 para lam≠1)
Salidas: ../out/diag_goldstone.csv + ../out/diag_goldstone_manifest.json
"""
import os, sys, io, json, csv, hashlib, time, platform, contextlib
import numpy as np, scipy, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__))
ROOT = os.path.dirname(HERE)
OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver')); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):          # silencia autotests de import
    import fd_bdg_operator as F
    import correlations as C
    import gate0_pseudoherm as G0M

LAMS = (0.25, 0.5, 1.0, 2.0, 4.0)
W = 2.0

def bg(x):
    return G0M.bg(x, W, 'horizon')        # n, v, V con Qp analitico (identico a mode_matching.bg)

def sqrtn_dd(x):
    """(sqrt n)'' analitico para el perfil tanh."""
    n_u = G0M.n_u; md = 2.0; n_d = (G0M.jcur/(md*np.sqrt(G0M.g)))**(2/3)
    t = np.tanh(x/W); Dn = n_d-n_u; n = n_u+Dn*0.5*(1+t)
    n1 = Dn*(1-t**2)/(2*W); n2 = -Dn*t*(1-t**2)/W**2
    return n2/(2*np.sqrt(n)) - n1**2/(4*n**1.5)

MARGIN = 10.0   # el residuo se mide SOLO en |x| <= L - MARGIN (las filas de borde de la FD truncan
                # la malla y el matching modal las sustituye; no son fisica). Fijado antes de mirar resultados v2.
def res(M, phi, mask):
    r = M @ phi
    m2 = np.concatenate([mask, mask])
    return float(np.linalg.norm(r[m2])/np.linalg.norm(phi[m2])), r

def Dmat(x):                                  # verbatim gate_pade_audit2.Dmat (D = -d²)
    N = len(x); dx = x[1]-x[0]; e = np.ones(N-1)
    return -(np.diag(e, 1) - 2*np.diag(np.ones(N)) + np.diag(e, -1))/dx**2

def Kpade2(x, lam):                           # verbatim gate_pade_audit2.Kpade2
    n, v, V = G0M.bg(x, 2.0, 'horizon'); G = G0M.g*n
    D = Dmat(x); B = np.diag(2*G) + D; Binv = np.linalg.inv(B)
    term1 = D @ Binv @ D
    term2 = (D@D) @ np.linalg.matrix_power(Binv, 3) @ (D@D)
    return 0.5*D + (lam-1)/4.0*term1 - (lam-1)**2/16.0*term2

def M_from_K(x, K):
    """BdG con cinetico K (denso), conveccion skew y diagonal G+Qp de gate0 (Qp lam-indep)."""
    N = len(x); dx = x[1]-x[0]; n, v, V = bg(x); e = np.ones(N-1)
    D1 = (np.diag(e, 1) - np.diag(e, -1))/(2*dx)
    conv = 0.5*(np.diag(v)@D1 + D1@np.diag(v))
    G = G0M.g*n; Qp = V - G0M.mu + G0M.g*n + 0.5*v**2
    Lp = K - 1j*conv + np.diag(G+Qp); Lm = K + 1j*conv + np.diag(G+Qp)
    return np.block([[Lp, np.diag(G)], [-np.diag(G), -Lm]])

rows = []
def add(op, lam, N, L, r, pred=None, note=''):
    rows.append(dict(operador=op, lam=lam, N=N, L=L, residuo_rel=r,
                     prediccion=(np.nan if pred is None else pred),
                     cociente=(np.nan if (pred is None or pred == 0) else r/pred), nota=note))
    print(f"  {op:28s} lam={lam:4.2f} N={N:5d}  ||M Phi||/||Phi|| = {r:.3e}"
          + ("" if pred is None else f"   pred = {pred:.3e}  ratio = {r/pred if pred else float('nan'):.4f}")
          + (f"  [{note}]" if note else ""))

t0 = time.time()
print("DIAGNOSTICO GOLDSTONE  (Phi_G = (sqrt n, -sqrt n))")
# ---------- A, B, B' : operadores locales dispersos, varias N (escalado dx²) ----------
L = 40.0
for N in (2000, 4000, 8000):
    x = np.linspace(-L, L, N); n, v, V = bg(x); s = np.sqrt(n)
    phi = np.concatenate([s, -s]).astype(complex)
    mk = np.abs(x) <= L - MARGIN
    sdd = sqrtn_dd(x); base = np.linalg.norm(sdd[mk])/np.linalg.norm(s[mk])
    rA, _ = res(F.build_MBdG(x, n, v, V, G0M.mu), phi, mk)
    add('A v1.0 build_MBdG', 1.0, N, L, rA)
    Qp = V - C.mu + C.g*n + 0.5*v**2
    for lam in LAMS:
        MB = C.build_op_lambda(x, n, v, V, lam)
        pred = abs(1-np.sqrt(lam))/2*base
        rB, _ = res(MB, phi, mk); add('B build_op_lambda (Qp fijo)', lam, N, L, rB, pred if lam != 1 else None)
        d = (np.sqrt(lam)-1)*Qp
        MBc = MB + sp.bmat([[sp.diags(d), None], [None, sp.diags(-d)]], format='csr')
        rBc, _ = res(MBc, phi, mk); add("B' idem, Qp_lam=sqrt(lam)Qp", lam, N, L, rBc)
    # separacion de terminos (flujo) en lam=2, N actual: conveccion sobre sqrt(n)
    e = np.ones(N); dx = x[1]-x[0]
    D1 = sp.diags([-e[:-1], e[:-1]], [-1, 1], (N, N))/(2*dx)
    conv = sp.diags(v)@D1 + 0.5*sp.diags(np.gradient(v, x))
    rows.append(dict(operador='termino conveccion (v d + v\'/2) sqrt n', lam=np.nan, N=N, L=L,
                     residuo_rel=float(np.linalg.norm((conv@s)[mk])/np.linalg.norm(s[mk])),
                     prediccion=0.0, cociente=np.nan, nota='n v = j constante -> 0 salvo dx^2'))
    print(f"  termino de conveccion sobre sqrt n: {rows[-1]['residuo_rel']:.3e}")

# ---------- C, D : operadores densos (Weyl no local, Pade) ----------
for (Ld, Nd) in ((30.0, 400), (30.0, 800)):
    x = np.linspace(-Ld, Ld, Nd); n, v, V = bg(x); s = np.sqrt(n)
    phi = np.concatenate([s, -s]).astype(complex); mk = np.abs(x) <= Ld - MARGIN
    K1w = G0M.weyl_kinetic(x, 1.0, W, 'horizon'); K1p = Kpade2(x, 1.0)
    for lam in LAMS:
        Mw = G0M.build_nonlocal_M(x, lam, W, 'horizon')
        Kw = G0M.weyl_kinetic(x, lam, W, 'horizon')
        predw = np.linalg.norm(((Kw-K1w)@s)[mk])/np.linalg.norm(s[mk])
        rC, _ = res(Mw, phi, mk); add('C Weyl no local (Qp fijo)', lam, Nd, Ld, rC, predw if lam != 1 else None)
        Kp = Kpade2(x, lam); Mp = M_from_K(x, Kp)
        predp = np.linalg.norm(((Kp-K1p)@s)[mk])/np.linalg.norm(s[mk])
        rD, _ = res(Mp, phi, mk); add('D Pade 2o orden (Qp fijo)', lam, Nd, Ld, rD, predp if lam != 1 else None)

secs = time.time()-t0
os.makedirs(OUT, exist_ok=True)
pcsv = os.path.join(OUT, 'diag_goldstone.csv')
with open(pcsv, 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
deps = {p: sha(os.path.join(ROOT, p)) for p in
        ('solver/fd_bdg_operator.py', 'solver/correlations.py', 'residual_anatomy/gate0_pseudoherm.py',
         'residual_anatomy/gate_pade_audit2.py', 'residual_anatomy/diag_goldstone.py')}
man = dict(tipo='DIAGNOSTICO (no gate)', version='v2 (residuo solo en interior |x|<=L-MARGIN; v1 dominada por filas de borde)', margin=MARGIN, script='residual_anatomy/diag_goldstone.py',
           timestamp=time.strftime('%Y-%m-%d %H:%M:%S %Z'), segundos=secs, lams=LAMS, w=W,
           csv=dict(path='out/diag_goldstone.csv', sha256=sha(pcsv)), dependencias_sha256=deps,
           entorno=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                        platform=platform.platform()))
with open(os.path.join(OUT, 'diag_goldstone_manifest.json'), 'w') as f:
    json.dump(man, f, indent=1)
print(f"\nCSV -> {pcsv}  ({secs:.0f} s)")
