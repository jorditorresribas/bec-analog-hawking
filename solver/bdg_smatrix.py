"""
FASE 3B — Matriz S de scattering de Bogoliubov sobre el horizonte waterfall.
Criterio de aceptacion: S^dag eta S = eta (pseudo-unitaridad). SIN CAP, SIN termal.

Corriente de norma conservada (derivada y verificada analiticamente):
   J = Im(u* u') + V_f|u|^2 + Im(v* v') - V_f|v|^2
Se normaliza cada canal a |J|=1. eta = signo de la norma de Bogoliubov = sign(Omega).

Metodo: matriz de transferencia sobre el perfil real n(x),v(x),V(x), integrando
la region de transicion (fuera de +-X0 el fondo es homogeneo). Se excluye el modo
evanescente que CRECE en -inf (condicion de contorno fisica).
"""
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d

g = 1.0
d = np.load("out/waterfall_state.npz")
xg, n0g, v0g, Vg = d["x"], d["n"], d["v"], d["V"]
mu = float(d["mu"]); jcur = float(d["j"])
dvg = np.gradient(v0g, xg)
fn = interp1d(xg, n0g, kind='cubic'); fv = interp1d(xg, v0g, kind='cubic')
fV = interp1d(xg, Vg, kind='cubic'); fdv = interp1d(xg, dvg, kind='cubic')

nL, vL, VL = n0g[:50].mean(), v0g[:50].mean(), Vg[:50].mean()
nR, vR, VR = n0g[-50:].mean(), v0g[-50:].mean(), Vg[-50:].mean()
cL, cR = np.sqrt(g*nL), np.sqrt(g*nR)


def A_at(x, omega):
    n, v, V, dv = float(fn(x)), float(fv(x)), float(fV(x)), float(fdv(x))
    a_uu = -1j*dv + v**2 + 4*g*n + 2*V - 2*mu - 2*omega
    b_vv = 1j*dv + v**2 + 4*g*n + 2*V - 2*mu + 2*omega
    return np.array([[0,1,0,0],
                     [a_uu, -2j*v, 2*g*n, 0],
                     [0,0,0,1],
                     [2*g*n, 0, b_vv, 2j*v]], dtype=complex)


def A_homog(n, v, V, omega):
    a_uu = v**2 + 4*g*n + 2*V - 2*mu - 2*omega
    b_vv = v**2 + 4*g*n + 2*V - 2*mu + 2*omega
    return np.array([[0,1,0,0],[a_uu,-2j*v,2*g*n,0],
                     [0,0,0,1],[2*g*n,0,b_vv,2j*v]], dtype=complex)


def norm_current(Y, Vf):
    u, up, v, vp = Y
    return (np.imag(np.conj(u)*up) + Vf*abs(u)**2
            + np.imag(np.conj(v)*vp) - Vf*abs(v)**2)


def modes(n, v, V, c, omega, Vf):
    """Devuelve los 4 modos clasificados y normalizados a |J|=1 (los propagantes)."""
    lam, vec = la.eig(A_homog(n, v, V, omega))
    out = []
    for i in range(4):
        k = -1j*lam[i]
        Y = vec[:, i]
        prop = abs(k.imag) < 1e-8
        if prop:
            kr = k.real
            Om = omega - v*kr
            vg = v + (c**2*kr + kr**3/2)/Om
            J = norm_current(Y, Vf)
            Y = Y/np.sqrt(abs(J))                 # normalizar a |J|=1
            out.append(dict(k=kr, prop=True, Om=Om, vg=vg, J=np.sign(J),
                            norm=int(np.sign(Om)), Y=Y))
        else:
            out.append(dict(k=k, prop=False, decay=(k.imag < 0), Y=Y))
    return out


def transfer(omega, X0=12.0, dx=0.002):
    """Matriz de transferencia Y(+X0) = Prop Y(-X0) sobre el perfil real."""
    M = np.eye(4, dtype=complex)
    x = -X0; n = int(2*X0/dx)
    for _ in range(n):
        k1 = A_at(x, omega)@M
        k2 = A_at(x+dx/2, omega)@(M+dx/2*k1)
        k3 = A_at(x+dx/2, omega)@(M+dx/2*k2)
        k4 = A_at(x+dx, omega)@(M+dx*k3)
        M = M + dx/6*(k1+2*k2+2*k3+k4)
        x += dx
    return M


def smatrix(omega, X0=12.0):
    Prop = transfer(omega, X0)
    mL = modes(nL, vL, VL, cL, omega, vL)
    mR = modes(nR, vR, VR, cR, omega, vR)
    # Ordenar LEFT: in(vg>0), out(vg<0), evan_decay(-inf), evan_grow(-inf)
    Lin  = next(m for m in mL if m['prop'] and m['vg'] > 0)
    Lout = next(m for m in mL if m['prop'] and m['vg'] < 0)
    Lev  = [m for m in mL if not m['prop']]
    Lev_d = next(m for m in Lev if m['decay']); Lev_g = next(m for m in Lev if not m['decay'])
    # RIGHT: 4 propagantes -> in(vg<0), out(vg>0)
    Rins  = [m for m in mR if m['prop'] and m['vg'] < 0]
    Routs = [m for m in mR if m['prop'] and m['vg'] > 0]
    PhiL = np.column_stack([Lin['Y'], Lout['Y'], Lev_d['Y'], Lev_g['Y']])
    PhiR = np.column_stack([Rins[0]['Y'], Rins[1]['Y'], Routs[0]['Y'], Routs[1]['Y']])
    T = la.solve(PhiR, Prop@PhiL)          # a_R = T a_L
    # a_L=(aLin,aLout,aLev_d,aLev_g=0); a_R=(aRin1,aRin2,aRout1,aRout2)
    # incognitas de matching (aLout,aLev_d) de las filas Rin; luego Rout.
    IN = ['Lin', 'Rin1', 'Rin2']; OUT = ['Lout', 'Rout1', 'Rout2']
    S = np.zeros((3, 3), dtype=complex)
    for col, src in enumerate(IN):
        aLin = 1.0 if src == 'Lin' else 0.0
        aRin1 = 1.0 if src == 'Rin1' else 0.0
        aRin2 = 1.0 if src == 'Rin2' else 0.0
        # filas 0,1 de T dan aRin1,aRin2 = T[r,0]aLin+T[r,1]aLout+T[r,2]aLev_d
        Mmat = np.array([[T[0,1], T[0,2]], [T[1,1], T[1,2]]])
        rhs = np.array([aRin1 - T[0,0]*aLin, aRin2 - T[1,0]*aLin])
        aLout, aLev_d = la.solve(Mmat, rhs)
        aL = np.array([aLin, aLout, aLev_d, 0.0])
        aRout1 = T[2,:]@aL; aRout2 = T[3,:]@aL
        S[:, col] = [aLout, aRout1, aRout2]
    eta_out = np.array([Lout['norm'], Routs[0]['norm'], Routs[1]['norm']])
    eta_in = np.array([Lin['norm'], Rins[0]['norm'], Rins[1]['norm']])
    return S, eta_in, eta_out, (Lout, Routs)


if __name__ == "__main__":
 print("="*68)
 print("FASE 3B — MATRIZ S DE SCATTERING DE BOGOLIUBOV")
 print("="*68)
 omega = 0.1
 S, eta_in, eta_out, info = smatrix(omega)
 print(f"\n  omega = {omega}")
 print(f"  canales IN  (norma): {list(eta_in)}")
 print(f"  canales OUT (norma): {list(eta_out)}")
 print(f"\n  |S| =\n{np.abs(S).round(4)}")
 Eout = np.diag(eta_out.astype(float)); Ein = np.diag(eta_in.astype(float))
 pu = np.max(np.abs(S.conj().T @ Eout @ S - Ein))
 print(f"\n  TEST S^dag eta S = eta:  max|residuo| = {pu:.2e}")
 print(f"  -> {'PASA: matriz S fisica' if pu < 1e-3 else 'FALLA'}")
