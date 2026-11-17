"""
DERIVACION GEOMETRICA de la gravedad superficial kappa desde la METRICA acustica
(no una formula recordada). ChatGPT: dejar explicitos flujo, tiempo, normalizacion
del vector de Killing. NO introducir/quitar factores por intuicion.

Metrica acustica 1+1D (Unruh/Visser):
    ds^2 = -(c_s^2 - v^2) dt^2 - 2 v dt dx + dx^2
Vector de Killing xi = d/dt (estacionario). Horizonte: g_tt=0 -> v^2 = c_s^2.
Gravedad superficial invariante:  kappa^2 = -1/2 (nabla_mu xi_nu)(nabla^mu xi^nu)|_H.

Se deriva kappa(c,v,c',v') simbolicamente, se compara con la candidata
|d(c_s-v)/dx|_H, y se evalua numericamente para cada anchura w comparando
kappa/2pi con T_fit (sin re-ajustar el espectro).
"""
import sympy as sp

t, x = sp.symbols('t x', real=True)
c = sp.Function('c')(x)      # c_s(x)
v = sp.Function('v')(x)      # v(x) velocidad del fluido

# metrica y su inversa
g = sp.Matrix([[-(c**2 - v**2), -v], [-v, sp.Integer(1)]])
ginv = g.inv()
coords = [t, x]

# Christoffel
def christoffel(g, ginv, coords):
    n = len(coords)
    Gamma = [[[0]*n for _ in range(n)] for _ in range(n)]
    for l in range(n):
        for m in range(n):
            for k in range(n):
                s = 0
                for r in range(n):
                    s += ginv[l, r]*(sp.diff(g[r, m], coords[k])
                                     + sp.diff(g[r, k], coords[m])
                                     - sp.diff(g[m, k], coords[r]))
                Gamma[l][m][k] = sp.simplify(s/2)
    return Gamma

Gamma = christoffel(g, ginv, coords)

# xi = d/dt  (componentes arriba), xi_low = g . xi
xi_up = sp.Matrix([1, 0])
xi_low = g*xi_up

# nabla_mu xi_nu = d_mu xi_nu - Gamma^l_{mu nu} xi_l
nab = sp.zeros(2, 2)
for mu in range(2):
    for nu in range(2):
        term = sp.diff(xi_low[nu], coords[mu])
        for l in range(2):
            term -= Gamma[l][mu][nu]*xi_low[l]
        nab[mu, nu] = sp.simplify(term)

# kappa^2 = -1/2 (nabla_mu xi_nu)(nabla^mu xi^nu)
k2 = 0
for mu in range(2):
    for nu in range(2):
        for a in range(2):
            for b in range(2):
                k2 += ginv[mu, a]*ginv[nu, b]*nab[mu, nu]*nab[a, b]
k2 = sp.simplify(-k2/2)

# en el horizonte v = c (flujo = sonido). Sustituir v -> c y v' -> (dv/dx)
cp, vp = sp.symbols("cp vp", real=True)
subs_deriv = {sp.diff(c, x): cp, sp.diff(v, x): vp}
k2_h = k2.subs(subs_deriv)
k2_h = k2_h.subs(v, c)                       # horizonte v=c
k2_h = sp.simplify(k2_h)
kappa = sp.sqrt(sp.simplify(k2_h))
print("kappa^2 en el horizonte (v=c):")
sp.pprint(sp.simplify(k2_h))
print("\nkappa geometrica =")
sp.pprint(sp.simplify(kappa))
print("\ncandidata |d(c-v)/dx| en horizonte = |cp - vp|")
print("ratio kappa_geom / |cp-vp| :")
sp.pprint(sp.simplify(kappa/sp.Abs(cp - vp)))

# ------- evaluacion numerica para cada anchura w -------
import numpy as np
g_, n_u, mu_u, md, = 1.0, 1.0, 0.5, 2.0
vLb = mu_u*np.sqrt(g_*n_u); jc = n_u*vLb; mu_ = g_*n_u + 0.5*vLb**2
nd = (jc/(md*np.sqrt(g_)))**(2/3)
# kappa geometrica como funcion lambdificada de (c, cp, vp)
cc, ccp, vvp = sp.symbols('cc ccp vvp', real=True)
kappa_expr = sp.simplify(kappa).subs({c: cc, cp: ccp, vp: vvp})
kf = sp.lambdify((cc, ccp, vvp), kappa_expr, 'numpy')

print("\n  w   c_H     cp_H     vp_H    |cp-vp|   kappa_geom  T_g(cand)  T_g(geom)")
Tfit = {2:0.0457, 3:0.0310, 4:0.0234, 6:0.0158, 8:0.0119}
for w in (2.0,3.0,4.0,6.0,8.0):
    xs=np.linspace(-8*w-12,8*w+12,20000)
    tt=np.tanh(xs/w); s=0.5*(1+tt); n=n_u+(nd-n_u)*s; vv=jc/n; cs=np.sqrt(g_*n)
    f=cs-vv; iH=np.argmin(np.abs(f))
    cH=cs[iH]; cpH=np.gradient(cs,xs)[iH]; vpH=np.gradient(vv,xs)[iH]
    cand=abs(cpH-vpH); kg=float(kf(cH,cpH,vpH))
    print(f" {w:3.0f} {cH:7.4f} {cpH:+8.4f} {vpH:+8.4f} {cand:9.4f} {kg:11.4f} "
          f"{cand/(2*np.pi):9.4f} {kg/(2*np.pi):9.4f}   T_fit={Tfit[int(w)]:.4f} "
          f"ratio_geom={Tfit[int(w)]/(kg/(2*np.pi)):.3f}")
