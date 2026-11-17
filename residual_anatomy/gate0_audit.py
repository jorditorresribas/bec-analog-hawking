"""
AUDITORIA del GATE 0 (cuestiones A, B, C de ChatGPT). Reutiliza las funciones de
gate0_pseudoherm.py sin re-ejecutar su bloque main.
"""
import numpy as np
src = open('residual_anatomy/gate0_pseudoherm.py').read().split('# ===========')[0]
ns = {}; exec(src, ns)
weyl_kinetic = ns['weyl_kinetic']; build_nonlocal_M = ns['build_nonlocal_M']
pseudoherm_residual = ns['pseudoherm_residual']; g = ns['g']; nval = ns['nval']

print("="*68); print("AUDITORIA GATE 0 (A: operador inhomogeneo, B: residuo, C: bordes)"); print("="*68)

# [A] En lam=1, eps_1(x,k)=k^2/2 NO depende de x -> K_weyl(lam=1) debe ser
#     independiente del fondo y reproducir -1/2 d^2/dx^2.
print("\n[A] Operador inhomogeneo de Weyl, control lam=1 (eps_1=k^2/2 indep de x):")
x = np.linspace(-40, 40, 600); dx = x[1]-x[0]
Kh = weyl_kinetic(x, 1.0, mode='horizon')
Kg = weyl_kinetic(x, 1.0, mode='homog')
print(f"   ||K(lam1,horizonte) - K(lam1,homog)||/||K|| = {np.abs(Kh-Kg).max()/np.abs(Kg).max():.2e}"
      f"   (debe ~0: el simbolo no usa n)")
f = np.exp(-(x/8.0)**2)*np.cos(0.5*x)
d2 = np.gradient(np.gradient(f, x), x); ref = -0.5*d2
sl = slice(120, 480)
print(f"   K(lam1) f  vs  -1/2 f''  (interior):  err rel = "
      f"{np.abs((Kg@f)[sl]-ref[sl]).max()/np.abs(ref[sl]).max():.2e}")

# [A''] consistencia LOCAL en fondo con HORIZONTE, region asintotica (n~n_u)
print("\n[A''] Consistencia local en fondo con HORIZONTE (lam=2), asintotica upstream x~-50:")
x = np.linspace(-80, 80, 1200)
K = weyl_kinetic(x, 2.0, mode='horizon')
i0 = int(np.argmin(np.abs(x-(-50.0))))
for kphys in (0.5, 1.0):
    u = np.exp(1j*kphys*x); Ku = K @ u
    nloc = nval(x[i0]); G = g*nloc; c2 = g*nloc
    es = -G + np.sqrt(G**2 + c2*kphys**2 + 2.0*kphys**4/4)
    ratio = np.real(np.mean((Ku/u)[i0-40:i0+40]))
    print(f"   k={kphys}: K.u/u={ratio:.6f}  eps_lam(k;n_loc)={es:.6f}  err={abs(ratio-es)/es:.2e}")

# [C] pseudo-herm con L creciente a dx CONSTANTE (N escala con L): borde vs bulk
print("\n[C] Pseudo-herm con L a dx~CONSTANTE (N escala con L), lam=2:")
print(f"   {'L':>5} {'N':>6} {'dx':>7} {'residuo rel':>13}")
for L in (20, 40, 80):
    N = int(round(2*L/0.2))
    xx = np.linspace(-L, L, N); M = build_nonlocal_M(xx, 2.0)
    print(f"   {L:5d} {N:6d} {xx[1]-xx[0]:7.4f} {pseudoherm_residual(M):13.2e}")
print("   (si BAJA con L a dx fijo -> residuo de BORDE diluido; si se mantiene -> bulk)")

print("\n[B] Residuo = max|M_dag eta - eta M| / max|M|, eta=diag(I_N,-I_N), sobre M COMPLETO.")
