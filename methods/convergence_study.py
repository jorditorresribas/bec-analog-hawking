"""
CONVERGENCIA HIDRODINAMICA — T_H_fit/T_H_grav vs anchura del horizonte w.

Si al suavizar el horizonte (w grande) el cociente -> 1 y el espectro se vuelve
termico, confirma que la desviacion a w=2 era DISPERSIVA (fisica) y valida el
pipeline. Fondo ANALITICO (sin interpolacion) parametrizado por w; para cada w se
reconstruye el fondo, se recalcula gravedad superficial y se rehace la matriz S.

Aviso: la matriz de transferencia tiene limite de conditioning (evanescentes);
para w grande el dominio es largo y la pseudo-unitaridad se degrada -> solo se
confia en los w donde S^dag eta S = eta se mantiene.
"""
import numpy as np
import scipy.linalg as la
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g = 1.0
n_u, mach_u, mach_d = 1.0, 0.5, 2.0
vL = mach_u*np.sqrt(g*n_u); jcur = n_u*vL; mu = g*n_u + 0.5*vL**2
n_d = (jcur/(mach_d*np.sqrt(g)))**(2/3)
nL, cL = n_u, np.sqrt(g*n_u)
nR, vR, cR = n_d, jcur/n_d, np.sqrt(g*n_d)
VL, VR = 0.0, mu - g*n_d - 0.5*(jcur/n_d)**2


def bg(x, w):
    """Fondo analitico waterfall: n,v,V,v' en x (metodo inverso)."""
    t = np.tanh(x/w); s = 0.5*(1+t)
    Dn = n_d - n_u
    n = n_u + Dn*s
    n1 = Dn*(1-t**2)/(2*w)
    n2 = -Dn*t*(1-t**2)/w**2
    v = jcur/n; dv = -jcur*n1/n**2
    Qp = 0.5*(n2/(2*n) - n1**2/(4*n**2))
    V = mu - g*n - 0.5*v**2 + Qp
    return n, v, V, dv


def A_at(x, omega, w):
    n, v, V, dv = bg(x, w)
    a_uu = -1j*dv + v**2 + 4*g*n + 2*V - 2*mu - 2*omega
    b_vv = 1j*dv + v**2 + 4*g*n + 2*V - 2*mu + 2*omega
    return np.array([[0,1,0,0],[a_uu,-2j*v,2*g*n,0],
                     [0,0,0,1],[2*g*n,0,b_vv,2j*v]], dtype=complex)


def A_homog(n, v, V, omega):
    a_uu = v**2 + 4*g*n + 2*V - 2*mu - 2*omega
    b_vv = v**2 + 4*g*n + 2*V - 2*mu + 2*omega
    return np.array([[0,1,0,0],[a_uu,-2j*v,2*g*n,0],
                     [0,0,0,1],[2*g*n,0,b_vv,2j*v]], dtype=complex)


def ncur(Y, Vf):
    u, up, v, vp = Y
    return np.imag(np.conj(u)*up)+Vf*abs(u)**2+np.imag(np.conj(v)*vp)-Vf*abs(v)**2


def mode_set(n, v, V, c, omega, Vf):
    lam, vec = la.eig(A_homog(n, v, V, omega))
    out = []
    for i in range(4):
        k = -1j*lam[i]; Y = vec[:, i]
        if abs(k.imag) < 1e-8:
            Om = omega - v*k.real; vg = v + (c**2*k.real + k.real**3/2)/Om
            Y = Y/np.sqrt(abs(ncur(Y, Vf)))
            out.append(dict(prop=True, vg=vg, norm=int(np.sign(Om)), Y=Y))
        else:
            out.append(dict(prop=False, decay=(k.imag < 0), Y=Y))
    return out


def smatrix(omega, w, X0, dx):
    M = np.eye(4, dtype=complex); x = -X0; ns = int(2*X0/dx)
    for _ in range(ns):
        k1 = A_at(x, omega, w)@M; k2 = A_at(x+dx/2, omega, w)@(M+dx/2*k1)
        k3 = A_at(x+dx/2, omega, w)@(M+dx/2*k2); k4 = A_at(x+dx, omega, w)@(M+dx*k3)
        M = M + dx/6*(k1+2*k2+2*k3+k4); x += dx
    mL = mode_set(nL, vL, VL, cL, omega, vL)
    mR = mode_set(nR, vR, VR, cR, omega, vR)
    Lin = next(m for m in mL if m['prop'] and m['vg'] > 0)
    Lout = next(m for m in mL if m['prop'] and m['vg'] < 0)
    Lev = [m for m in mL if not m['prop']]
    Lev_d = next(m for m in Lev if m['decay']); Lev_g = next(m for m in Lev if not m['decay'])
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    Rout = [m for m in mR if m['prop'] and m['vg'] > 0]
    PhiL = np.column_stack([Lin['Y'], Lout['Y'], Lev_d['Y'], Lev_g['Y']])
    PhiR = np.column_stack([Rin[0]['Y'], Rin[1]['Y'], Rout[0]['Y'], Rout[1]['Y']])
    T = la.solve(PhiR, M@PhiL)
    S = np.zeros((3, 3), dtype=complex)
    for col, src in enumerate(['Lin', 'Rin1', 'Rin2']):
        aLin = 1.0*(src == 'Lin'); aR1 = 1.0*(src == 'Rin1'); aR2 = 1.0*(src == 'Rin2')
        Mm = np.array([[T[0,1], T[0,2]], [T[1,1], T[1,2]]])
        aLout, aLev = la.solve(Mm, [aR1 - T[0,0]*aLin, aR2 - T[1,0]*aLin])
        aL = np.array([aLin, aLout, aLev, 0.0])
        S[:, col] = [aLout, T[2,:]@aL, T[3,:]@aL]
    eta_in = np.array([Lin['norm'], Rin[0]['norm'], Rin[1]['norm']])
    eta_out = np.array([Lout['norm'], Rout[0]['norm'], Rout[1]['norm']])
    return S, eta_in, eta_out


def T_grav(w):
    xs = np.linspace(-5*w-5, 5*w+5, 4000)
    n, v, V, dv = bg(xs, w); c = np.sqrt(g*n)
    f = c - v; iH = np.argmin(np.abs(f))
    ksg = 0.5*abs(np.gradient(f, xs)[iH])
    return ksg/(2*np.pi), ksg


print("="*70)
print("CONVERGENCIA HIDRODINAMICA: T_H_fit/T_H_grav vs anchura w")
print("="*70)
widths = [2.0, 3.0, 4.0, 5.0]
res = []
for w in widths:
    X0 = 12.0; dx = 0.004           # X0 fijo: conditioning estable; mas w = menos homogeneo
    THg, ksg = T_grav(w)
    oms = np.linspace(0.02, 0.30, 14)
    Ntot, pus = [], []
    for om in oms:
        try:
            S, ei, eo = smatrix(om, w, X0, dx)
        except Exception:
            Ntot.append(np.nan); pus.append(np.nan); continue
        pu = np.max(np.abs(S.conj().T@np.diag(eo.astype(float))@S - np.diag(ei.astype(float))))
        jneg = np.where(ei < 0)[0]; ipos = np.where(eo > 0)[0]
        N = sum(sum(abs(S[i, jj])**2 for jj in jneg) for i in ipos)
        Ntot.append(N); pus.append(pu)
    oms, Ntot, pus = np.array(oms), np.array(Ntot), np.array(pus)
    ok = (pus < 1e-2) & np.isfinite(Ntot) & (Ntot > 0)
    if ok.sum() >= 4:
        lo = ok & (oms < 0.18)
        A = np.polyfit(oms[lo], np.log(1+1/Ntot[lo]), 1)
        THf = 1/A[0]; ratio = THf/THg
    else:
        THf, ratio = np.nan, np.nan
    res.append((w, THg, THf, ratio, np.nanmedian(pus[ok]) if ok.any() else np.nan, ok.sum()))
    print(f"  w={w:4.1f}: T_grav={THg:.4f}  T_fit={THf:.4f}  ratio={ratio:.2f}  "
          f"pseudo-unit~{np.nanmedian(pus[ok]) if ok.any() else np.nan:.1e}  (pts:{ok.sum()})")

res = np.array([(r[0], r[1], r[2], r[3], r[4]) for r in res])
print("\n" + "="*70)
print("  Si ratio -> 1 al crecer w: desviacion era DISPERSIVA, pipeline validado.")
print("="*70)

fig, ax = plt.subplots(1, 2, figsize=(12, 5))
ax[0].plot(res[:,0], res[:,3], 'ko-', ms=7)
ax[0].axhline(1.0, color='g', ls='--', label='hidrodinamico (ratio=1)')
ax[0].set_xlabel('anchura horizonte w (/ξ)'); ax[0].set_ylabel(r'$T_H^{fit}/T_H^{grav}$')
ax[0].set_title('Convergencia al limite hidrodinamico'); ax[0].legend(); ax[0].grid(alpha=0.3)
ax[1].plot(res[:,0], res[:,1], 'bs-', label=r'$T_H^{grav}=\kappa/2\pi$')
ax[1].plot(res[:,0], res[:,2], 'ro-', label=r'$T_H^{fit}$ (del espectro)')
ax[1].set_xlabel('anchura horizonte w (/ξ)'); ax[1].set_ylabel(r'$T_H$')
ax[1].set_title('Temperaturas vs anchura'); ax[1].legend(); ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.savefig("out/convergence_study.png", dpi=120); plt.close()
print("\nFigura: out/convergence_study.png")
