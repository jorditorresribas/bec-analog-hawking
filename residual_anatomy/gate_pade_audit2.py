"""
AUDITORIA PADE fase 2 (secuencia de ChatGPT antes del solver):
  (2) Omega_P^2 vs c^2k^2+lam k^4/4 DIRECTAMENTE, malla densa, tol 1%.
  (3) Validar implementacion de D^2 (2G+D)^-3 D^2: directa (densa) vs 3 band-solves vs adjunto.
  (4) Padé(2o orden) vs Weyl, separando error de MALLA (lam=1) de error ESTRUCTURAL (lam=2).
Weyl (gate0) como referencia. Operadores densos para auditoria.
"""
import numpy as np, scipy.linalg as sla
src = open('residual_anatomy/gate0_pseudoherm.py').read().split('# ===========')[0]
ns = {}; exec(src, ns)
weyl_kinetic = ns['weyl_kinetic']; bg = ns['bg']; g = ns['g']
G0 = g*1.0

def eps_exact(k, lam): return -G0 + np.sqrt(G0**2 + G0*k**2 + lam*k**4/4)
def eps_P2(k, lam):
    s = k**2
    return s/2 + (lam-1)*s**2/(4*(2*G0+s)) - (lam-1)**2*s**4/(16*(2*G0+s)**3)

print("="*70); print("AUDITORIA PADE fase 2 (pruebas 2,3,4 de ChatGPT)"); print("="*70)

# --- (2) Omega_P^2 vs c^2k^2+lam k^4/4 DIRECTO ---
print("\n[2] Omega_P^2 = eps_P2(eps_P2+2G)  vs  Omega^2 = c^2 k^2 + lam k^4/4 (DIRECTO), tol 1%")
kb = np.linspace(0.02, 2.5, 500)
print(f"    {'lam':>5} {'k@1% (Om_P^2)':>14} {'max err(k<=1.6)':>16}")
for lam in (0.5, 1.0, 2.0):
    ep = eps_P2(kb, lam); Om2p = ep*(ep+2*G0)
    Om2_obj = G0*kb**2 + lam*kb**4/4             # c^2=G0; objetivo DIRECTO
    err = np.abs(Om2p-Om2_obj)/Om2_obj
    over = kb[err > 0.01]; kcut = f"{over[0]:.2f}" if len(over) else ">2.5"
    m = kb <= 1.6
    print(f"    {lam:5.1f} {kcut:>14} {err[m].max():16.2e}")

# --- (3) Implementacion de D^2 (2G+D)^-3 D^2 ---
def Dmat(x):
    N = len(x); dx = x[1]-x[0]; e = np.ones(N-1)
    return -(np.diag(e, 1) - 2*np.diag(np.ones(N)) + np.diag(e, -1))/dx**2   # D=-d^2
print("\n[3] Implementacion D^2(2G+D)^-3 D^2: directa(densa) vs 3 band-solves vs adjunto")
print(f"    {'N':>5} {'||direct-bandsolve||/||.||':>26} {'||T-T^H||/||T|| (adj)':>22}")
for N in (200, 400, 800):
    x = np.linspace(-40, 40, N); n, v, V = bg(x, 2.0, 'horizon'); G = g*n
    D = Dmat(x); B = np.diag(2*G) + D
    D2 = D @ D
    # directa (densa)
    T_direct = D2 @ np.linalg.matrix_power(np.linalg.inv(B), 3) @ D2
    # por 3 band-solves sobre un conjunto de vectores de prueba (aleatorios + localizados)
    rng = np.random.default_rng(0)
    Xtest = np.column_stack([rng.standard_normal(N), rng.standard_normal(N),
                             np.exp(-(x/3.0)**2)*np.cos(x), np.exp(-(x/15.0)**2)])
    Y = D2 @ Xtest
    Z = sla.solve(B, Y); Z = sla.solve(B, Z); Z = sla.solve(B, Z)   # (2G+D)^-3
    T_bs = D2 @ Z
    rel = np.linalg.norm(T_direct @ Xtest - T_bs)/np.linalg.norm(T_direct @ Xtest)
    adj = np.abs(T_direct - T_direct.conj().T).max()/np.abs(T_direct).max()
    print(f"    {N:5d} {rel:26.2e} {adj:22.2e}")

# --- (4) Padé(2o orden) vs Weyl, separar malla (lam=1) de estructural (lam=2) ---
def Kpade2(x, lam):
    n, v, V = bg(x, 2.0, 'horizon'); G = g*n
    D = Dmat(x); B = np.diag(2*G) + D; Binv = np.linalg.inv(B)
    term1 = D @ Binv @ D                       # D(2G+D)^-1 D  = s^2/(2G+s)
    term2 = (D@D) @ np.linalg.matrix_power(Binv, 3) @ (D@D)  # D^2(2G+D)^-3 D^2
    return 0.5*D + (lam-1)/4.0*term1 - (lam-1)**2/16.0*term2
print("\n[4] Kpade2 vs Kweyl en fondo horizonte: lam=1 = error de MALLA; lam=2 = + estructural")
x = np.linspace(-60, 60, 600)
ihor = np.abs(x) < 6.0; iasy = np.abs(x) > 35.0
fs = np.exp(-(x/15.0)**2)*np.cos(0.5*x); fl = np.exp(-(x/3.0)**2)*np.cos(1.0*x)
base = {}
for lam in (1.0, 2.0):
    Kp = Kpade2(x, lam); Kw = weyl_kinetic(x, lam, mode='horizon')
    for name, f in (("suave", fs), ("loc.hor", fl)):
        d = (Kp@f)-(Kw@f)
        rh = np.linalg.norm(d[ihor])/max(np.linalg.norm((Kw@f)[ihor]), 1e-30)
        print(f"    lam={lam} {name:>8}: |Kp2-Kw| rel horizonte = {rh:.2e}")
        base[(lam, name)] = rh
print("    -> error ESTRUCTURAL ~ (error lam=2) menos (error de malla lam=1):")
for name in ("suave", "loc.hor"):
    est = max(base[(2.0, name)] - base[(1.0, name)], 0.0)
    print(f"       {name:>8}: estructural ~ {est:.2e}  (vs total lam=2 {base[(2.0,name)]:.2e})")
