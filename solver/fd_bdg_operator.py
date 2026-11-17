"""
SOLVER INDEPENDIENTE (paso 1) — operador BdG discretizado en diferencias finitas.

Metodo COMPLETAMENTE DISTINTO del transfer: se ensambla el operador BdG como
matriz dispersa 2N x 2N (sobre u,v en una malla) para un solve de scattering
disperso (estable, NO exponencia evanescentes). Diseñado para poder REFUTAR el
resultado del transfer si corresponde.

Operador BdG (hbar=m=1), eigenproblema  M_BdG [u;v] = omega [u;v]:
  fila u:  [ -1/2 (d_x+iv)^2 + W ] u + g n v
  fila v:  -g n u + [ -(-1/2 (d_x-iv)^2 + W) ] v,   W = 2 g n + V - mu
es decir  M = [[ L+ ,  g n ],[ -g n , -L- ]],  L± = -1/2 (d_x±iv)^2 + W.

VALIDACION (paso 1): en region homogenea, una onda plana e^{ikx}(U,V) con k raiz
del cuartico debe satisfacer M_BdG Psi = omega Psi (residuo ~ error de malla).
Esto valida el operador INDEPENDIENTEMENTE del transfer. (La extraccion de beta y
la comparacion con el transfer son el paso 2, en sesion dedicada.)
"""
import numpy as np
import scipy.sparse as sp

g = 1.0


def build_MBdG(x, n, v, V, mu):
    """Operador BdG disperso 2N x 2N con FD de 2o orden. u = [0:N], v = [N:2N]."""
    N = len(x); dx = x[1]-x[0]
    W = 2*g*n + V - mu
    dv = np.gradient(v, x)
    # operadores FD: D1 (central), D2 (laplaciano)
    e = np.ones(N)
    D1 = sp.diags([-e[:-1], e[:-1]], [-1, 1], (N, N)).tolil()/(2*dx)
    D2 = sp.diags([e[:-1], -2*e, e[:-1]], [-1, 0, 1], (N, N)).tolil()/dx**2
    I = sp.identity(N)
    # L+ = -1/2 (d_x+iv)^2 + W = -1/2 d_x^2 - i v d_x + (1/2 v^2 - 1/2 i v' + W)
    diagP = sp.diags(0.5*v**2 - 0.5j*dv + W)
    diagM = sp.diags(0.5*v**2 + 0.5j*dv + W)
    Lp = -0.5*D2 - 1j*sp.diags(v)@D1 + diagP   # -1/2(d+iv)^2+W
    Lm = -0.5*D2 + 1j*sp.diags(v)@D1 + diagM   # -1/2(d-iv)^2+W
    gnI = sp.diags(g*n)
    M = sp.bmat([[Lp, gnI], [-gnI, -Lm]], format='csr')
    return M


def plane_wave_mode(k, n, v, V, mu, omega):
    """Vector (U,V) de la onda plana e^{ikx} para la raiz k (null de M_sym-omega)."""
    W = 2*g*n + V - mu
    Lp = 0.5*(k+v)**2 + W          # simbolo de L+ en e^{ikx}
    Lm = 0.5*(k-v)**2 + W          # simbolo de L-
    Msym = np.array([[Lp, g*n], [-g*n, -Lm]])
    vals, vecs = np.linalg.eig(Msym)
    i = np.argmin(np.abs(vals - omega))   # el que da omega
    return vecs[:, i], vals[i]


def quartic_roots(omega, v, c):
    return np.roots([0.25, 0.0, c**2 - v**2, 2*omega*v, -omega**2])


# ---------------------- VALIDACION paso 1 (homogeneo) ------------------------
print("="*68)
print("SOLVER INDEPENDIENTE paso 1 — operador BdG en diferencias finitas")
print("="*68)
n_u, vL, mu = 1.0, 0.5, 1.125; VL = 0.0; cL = np.sqrt(g*n_u)
L, N = 40.0, 4000
x = np.linspace(-L/2, L/2, N)
n = np.full(N, n_u); v = np.full(N, vL); V = np.full(N, VL)
M = build_MBdG(x, n, v, V, mu)

omega = 0.1
print(f"\n  Homogeneo: n={n_u} v={vL} c={cL}  omega={omega}")
print(f"  Validando cada raiz propagante del cuartico contra M_BdG Psi = omega Psi:")
maxres = 0.0
for k in quartic_roots(omega, vL, cL):
    if abs(k.imag) > 1e-8:
        continue                          # solo propagantes
    kr = k.real
    (U, Vv), lam = plane_wave_mode(kr, n_u, vL, VL, mu, omega)
    Psi = np.concatenate([U*np.exp(1j*kr*x), Vv*np.exp(1j*kr*x)])
    MPsi = M @ Psi
    core = np.r_[np.abs(x) < 0.8*L/2, np.abs(x) < 0.8*L/2]   # interior (evita bordes)
    res = np.max(np.abs((MPsi - omega*Psi)[core])) / np.max(np.abs(omega*Psi[core]))
    maxres = max(maxres, res)
    print(f"    k={kr:+.4f}:  |M Psi - omega Psi|/|omega Psi| = {res:.2e}")
print(f"\n  residuo maximo = {maxres:.2e}")
print(f"  -> {'OK: operador FD independiente reproduce la dispersion (converge con dx)' if maxres<1e-3 else 'REVISAR'}")

# Convergencia con la malla (confirma que el residuo es error de discretizacion)
print("\n  Convergencia con resolucion (debe decrecer ~dx^2):")
for Nt in (2000, 4000, 8000):
    xt = np.linspace(-L/2, L/2, Nt)
    Mt = build_MBdG(xt, np.full(Nt, n_u), np.full(Nt, vL), np.full(Nt, VL), mu)
    k = [r.real for r in quartic_roots(omega, vL, cL) if abs(r.imag) < 1e-8][0]
    (U, Vv), _ = plane_wave_mode(k, n_u, vL, VL, mu, omega)
    Psi = np.concatenate([U*np.exp(1j*k*xt), Vv*np.exp(1j*k*xt)])
    core = np.r_[np.abs(xt) < 0.8*L/2, np.abs(xt) < 0.8*L/2]
    res = np.max(np.abs((Mt@Psi - omega*Psi)[core]))/np.max(np.abs(omega*Psi[core]))
    print(f"    N={Nt:5d} (dx={L/Nt:.4f}): residuo = {res:.2e}")
print("\n  (Paso 2, sesion dedicada: solve de scattering con fuente + BC salientes,")
print("   extraer beta INDEPENDIENTE y comparar con el transfer.)")
