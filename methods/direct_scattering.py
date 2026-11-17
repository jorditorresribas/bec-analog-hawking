"""
PASO 2 — extraccion INDEPENDIENTE de beta por solve de scattering directo.
Metodo distinto del transfer: solve disperso sobre el operador FD validado
(fd_bdg_operator), fuente inyectada + condiciones salientes (CAP), beta por
ajuste del campo a los modos asintoticos. MISMA normalizacion de flujo que el
transfer, implementacion distinta. Debe poder REFUTAR el transfer.

beta_directo (|beta|^2) se compara con |S[neg-out, pos-in]|^2 del transfer, con
convergencia en N y dominio. NO se ajusta uno para reproducir al otro.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from fd_bdg_operator import build_MBdG, plane_wave_mode, quartic_roots, g

# ---- parametros del fondo (identicos al transfer / waterfall, w=2) ----------
n_u, mach_u, mach_d, w = 1.0, 0.5, 2.0, 2.0
vLb = mach_u*np.sqrt(g*n_u); jcur = n_u*vLb; mu = g*n_u + 0.5*vLb**2
n_d = (jcur/(mach_d*np.sqrt(g)))**(2/3)


def bg(x):
    t = np.tanh(x/w); s = 0.5*(1+t); Dn = n_d-n_u
    n = n_u+Dn*s; n1 = Dn*(1-t**2)/(2*w); n2 = -Dn*t*(1-t**2)/w**2
    v = jcur/n; dv = -jcur*n1/n**2; Qp = 0.5*(n2/(2*n)-n1**2/(4*n**2))
    return n, v, mu-g*n-0.5*v**2+Qp


def flux_mode(k, n, v, V, omega):
    """(U,V) normalizado a |flujo de norma|=1, + grupo + norma."""
    (U, Vv), lam = plane_wave_mode(k, n, v, V, mu, omega)
    J = (k+v)*abs(U)**2 + (k-v)*abs(V)**2     # corriente de norma de la onda plana
    U, Vv = U/np.sqrt(abs(J)), Vv/np.sqrt(abs(J))
    Om = omega - v*k
    vg = v + (g*n*k + k**3/2)/Om          # c_s^2 = g n
    return U, Vv, vg, int(np.sign(Om))


def asymptotic_modes(n, v, V, omega):
    c = np.sqrt(g*n); out = []
    for k in quartic_roots(omega, v, c):
        if abs(k.imag) < 1e-8:
            U, Vv, vg, nrm = flux_mode(k.real, n, v, V, omega)
            out.append(dict(k=k.real, prop=True, U=U, V=Vv, vg=vg, norm=nrm))
        else:
            out.append(dict(k=k, prop=False, decay=(k.imag < 0)))
    return out


def direct_beta(L=20.0, N=4000, cap0=3.0, verbose=False):
    x = np.linspace(-L, L, N); dx = x[1]-x[0]
    n, v, V = bg(x)
    nL, vL_, VL = n[0], v[0], V[0]; nR, vR, VR = n[-1], v[-1], V[-1]
    M = build_MBdG(x, n, v, V, mu)

    mL = asymptotic_modes(nL, vL_, VL, omega)
    mR = asymptotic_modes(nR, vR, VR, omega)
    Lin = next(m for m in mL if m['prop'] and m['vg'] > 0)   # entrante upstream (+)

    # onda incidente con taper (1 en la izquierda homogenea -> 0 antes del horizonte)
    taper = 0.5*(1 - np.tanh((x+8.0)/1.2))
    phase = np.exp(1j*Lin['k']*x)
    psi_inc = np.concatenate([Lin['U']*phase*taper, Lin['V']*phase*taper])

    # CAP (condiciones salientes): franja de anchura FIJA en los bordes
    xcap = L - 6.0
    Gam = np.where(np.abs(x) > xcap, cap0*((np.abs(x)-xcap)/(L-xcap))**2, 0.0)
    CAP = sp.diags(np.concatenate([Gam, Gam]))
    Op = (M - omega*sp.identity(2*N) - 1j*CAP).tocsc()

    src = -(M - omega*sp.identity(2*N)) @ psi_inc
    psi_scat = spla.spsolve(Op, src)
    psi = psi_inc + psi_scat
    u = psi[:N]; vv = psi[N:]

    # ajustar el campo TRANSMITIDO (region derecha homogenea, antes del CAP) a los
    # 4 modos downstream -> amplitudes c_j (modos normalizados a flujo 1)
    win = (x > 9) & (x < 13)
    props = [m for m in mR if m['prop']]
    Phi_u = np.array([m['U']*np.exp(1j*m['k']*x[win]) for m in props]).T
    Phi_v = np.array([m['V']*np.exp(1j*m['k']*x[win]) for m in props]).T
    Amat = np.vstack([Phi_u, Phi_v]); b = np.concatenate([u[win], vv[win]])
    c, *_ = np.linalg.lstsq(Amat, b, rcond=None)
    fit_res = np.linalg.norm(Amat@c - b)/np.linalg.norm(b)

    # beta = amplitud del canal SALIENTE de norma NEGATIVA (vg>0, norm<0)
    beta2 = 0.0; trans = {}
    for ci, m in zip(c, props):
        trans[(m['norm'], np.sign(m['vg']))] = abs(ci)**2
        if m['vg'] > 0 and m['norm'] < 0:
            beta2 += abs(ci)**2
    if verbose:
        print(f"    residuo del ajuste a modos = {fit_res:.2e}")
        print(f"    amplitudes^2 por (norma,dir): {[(k,round(val,4)) for k,val in trans.items()]}")
    return beta2, fit_res


omega = 0.1
print("="*66)
print("PASO 2 — beta INDEPENDIENTE por scattering directo (FD, estable)")
print("="*66)

# beta del TRANSFER (metodo anterior), misma definicion
import bdg_smatrix as B
S, ei, eo = B.smatrix(omega)[:3]
ineg_out = np.where(eo < 0)[0][0]; ipos_in = np.where(ei > 0)[0][0]
beta2_transfer = abs(S[ineg_out, ipos_in])**2
print(f"\n  beta^2 TRANSFER  (|S[neg-out, pos-in]|^2) = {beta2_transfer:.5f}")

print(f"\n  beta^2 DIRECTO, convergencia en N (L=20, mismo omega={omega}):")
for N in (2000, 4000, 8000):
    b2, r = direct_beta(L=20.0, N=N, verbose=(N == 4000))
    print(f"    N={N:5d}: beta^2_directo = {b2:.5f}   (res.ajuste {r:.1e})")

print(f"\n  convergencia en tamano de dominio (N=4000):")
for L in (16.0, 20.0, 26.0):
    b2, r = direct_beta(L=L, N=int(L*200))
    print(f"    L={L:4.0f}: beta^2_directo = {b2:.5f}")

b2f, _ = direct_beta(L=20.0, N=8000)
print("\n" + "="*66)
print(f"  COMPARACION:  transfer={beta2_transfer:.5f}  directo={b2f:.5f}  "
      f"ratio={b2f/beta2_transfer:.3f}")
ok = abs(b2f-beta2_transfer)/beta2_transfer < 0.1
print(f"  -> {'ESCENARIO A: dos metodos coinciden -> beta validado' if ok else 'DISCREPAN -> revisar (escenario B/C)'}")
print("="*66)
