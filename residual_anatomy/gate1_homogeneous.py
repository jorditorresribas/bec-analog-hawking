"""
GATE 1 — validacion COMPLETA del operador no local en fondo HOMOGENEO (antes del horizonte).
Fondo upstream CONSISTENTE: n=1, v=vLb=0.5, V=0 -> Qp=0, G=c^2=1 (ChatGPT).
No toca el benchmark congelado. Reutiliza build_nonlocal_M/weyl_kinetic de gate0_pseudoherm.py.

Criterios (ChatGPT):
  1. dispersion homogenea Omega^2 = c^2 k^2 + lam k^4/4
  2. autovalores del simbolo BdG = v k ± Omega(k)
  3. modos: |U|^2-|V|^2 = ±1, signo de norma correcto
  4. corriente J_inf = eps'(k)(|U|^2+|V|^2)+v(|U|^2-|V|^2) == v_g (|U|^2-|V|^2)
  5. lam=1 recupera Bogoliubov estandar
  6. convergencia con N (parte discreta, k~=sin(k dx)/dx -> k)
"""
import numpy as np
src = open('residual_anatomy/gate0_pseudoherm.py').read().split('# ===========')[0]
ns = {}; exec(src, ns)
weyl_kinetic = ns['weyl_kinetic']; build_nonlocal_M = ns['build_nonlocal_M']
g = ns['g']; vLb = ns['vLb']

# --- fondo homogeneo consistente ---
n0, v0 = 1.0, vLb; G0 = g*n0; c2 = g*n0; c0 = np.sqrt(c2)

def eps_lam(k, lam): return -G0 + np.sqrt(G0**2 + c2*k**2 + lam*k**4/4.0)
def Omega(k, lam):   e = eps_lam(k, lam); return np.sqrt(e*(e + 2*G0))
def epsp(k, lam):    return (c2*k + lam*k**3/2)/(eps_lam(k, lam) + G0)   # eps'(k)
def Msym(k, lam):
    e = eps_lam(k, lam)
    return np.array([[e+G0+v0*k, G0], [-G0, -e-G0+v0*k]], dtype=complex)

print("="*70); print("GATE 1 — fondo HOMOGENEO consistente (n=1,v=0.5,V=0,Qp=0,G=c^2=1)"); print("="*70)

KS = [-2.0, -1.0, -0.5, 0.5, 1.0, 2.0]        # k±, lejos de Nyquist
LAMS = [0.5, 1.0, 2.0, 4.0]

# ---- 1+2: dispersion y autovalores del simbolo ----
print("\n[1+2] Dispersion Omega^2=c^2k^2+lam k^4/4  y  autovalores de M(k) = v k ± Omega")
print(f"    {'lam':>4} {'k':>6} {'Om^2_obj':>10} {'eps(eps+2G)':>12} {'eig1-(vk+Om)':>13} {'eig2-(vk-Om)':>13}")
emax = 0.0
for lam in LAMS:
    for k in KS:
        Om2_obj = c2*k**2 + lam*k**4/4
        e = eps_lam(k, lam); Om2_fac = e*(e+2*G0)
        w = np.sort(np.linalg.eigvals(Msym(k, lam)).real)
        tgt = np.sort([v0*k+Omega(k, lam), v0*k-Omega(k, lam)])
        d1, d2 = w[0]-tgt[0], w[1]-tgt[1]
        emax = max(emax, abs(Om2_obj-Om2_fac), abs(d1), abs(d2))
        if lam in (1.0, 2.0) and k in (0.5, 1.0, 2.0):
            print(f"    {lam:4.1f} {k:6.2f} {Om2_obj:10.5f} {Om2_fac:12.5f} {d1:13.2e} {d2:13.2e}")
print(f"    -> max|error| sobre todos (lam,k) = {emax:.2e}")

# ---- 3+4: norma de los modos y corriente ----
print("\n[3+4] Norma |U|^2-|V|^2=±1 y corriente  J_inf = v_g (|U|^2-|V|^2)")
print(f"    {'lam':>4} {'k':>6} {'rama':>5} {'norm=|U|^2-|V|^2':>16} {'J_inf':>11} {'vg*norm':>11} {'dif':>9}")
for lam in (1.0, 2.0):
    for k in (0.5, 1.0):
        w, V = np.linalg.eig(Msym(k, lam))
        for idx in range(2):
            U, Vv = V[:, idx]
            norm = abs(U)**2 - abs(Vv)**2
            s = 1.0/np.sqrt(abs(norm)); U, Vv = U*s, Vv*s      # normaliza a |norm|=1
            norm = abs(U)**2 - abs(Vv)**2
            om = w[idx].real
            branch = +1 if om > v0*k else -1
            vg = v0 + branch*(c2*k + lam*k**3/2)/Omega(k, lam)
            Jinf = epsp(k, lam)*(abs(U)**2+abs(Vv)**2) + v0*(abs(U)**2-abs(Vv)**2)
            print(f"    {lam:4.1f} {k:6.2f} {branch:+5d} {norm:16.6f} {Jinf:11.6f} {vg*norm:11.6f} {abs(Jinf-vg*norm):9.1e}")

# ---- 5: lam=1 recupera Bogoliubov estandar (eps_1=k^2/2) ----
print("\n[5] lam=1: eps_1(k)=k^2/2 exacto  ->  Bogoliubov estandar")
for k in (0.5, 1.0, 2.0):
    print(f"    k={k}: eps_1={eps_lam(k,1.0):.8f}  k^2/2={k**2/2:.8f}  dif={abs(eps_lam(k,1.0)-k**2/2):.1e}")

# ---- 6: parte DISCRETA — simbolo del operador discreto, convergencia con N ----
# Extraer M_disc(k) 2x2 aplicando M a (e^{ikx},0) y (0,e^{ikx}) y leyendo el interior.
# conveccion centrada -> k~ = sin(k dx)/dx (no exactamente k); debe converger a k.
print("\n[6] DISCRETO: autovalores del operador vs v*k~ ± Omega_disc, k~=sin(k dx)/dx; convergencia N")
def disc_symbol(x, lam, k):
    N = len(x); M = build_nonlocal_M(x, lam, mode='homog')
    e = np.exp(1j*k*x); z = np.zeros(N, complex)
    i0, i1 = N//2-30, N//2+30
    def blk(col):  # col=(u_part, v_part) plane wave; devuelve [[Mu_u,Mu_v],[Mv_u,Mv_v]] en interior
        psi = np.concatenate(col)
        r = M @ psi
        return r[:N], r[N:]
    ru_u, rv_u = blk((e, z))     # inyecta en u
    ru_v, rv_v = blk((z, e))     # inyecta en v
    sl = slice(i0, i1)
    Md = np.array([[np.mean((ru_u/e)[sl]), np.mean((ru_v/e)[sl])],
                   [np.mean((rv_u/e)[sl]), np.mean((rv_v/e)[sl])]])
    return Md
print(f"    {'N':>5} {'k':>6} {'eig_disc(+)':>12} {'vk+Omega(cont)':>14} {'|dif|':>9}")
for N in (400, 800, 1600):
    x = np.linspace(-40, 40, N); dx = x[1]-x[0]
    for k in (1.0,):
        Md = disc_symbol(x, 2.0, k)
        wd = np.sort(np.linalg.eigvals(Md).real)
        ktil = np.sin(k*dx)/dx
        cont = v0*k + Omega(k, 2.0)
        disc_pred = v0*ktil + Omega(k, 2.0)     # prediccion con k~ (convergente)
        print(f"    {N:5d} {k:6.2f} {wd[1]:12.6f} {cont:14.6f} {abs(wd[1]-cont):9.2e}"
              f"   (k~={ktil:.5f}, pred_disc={disc_pred:.6f}, |eig-pred|={abs(wd[1]-disc_pred):.1e})")
print("\n(nota: la discrepancia con el continuo viene de k~=sin(k dx)/dx -> k; converge O(dx^2).)")
