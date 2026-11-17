"""
JUEZ INDEPENDIENTE v2 — scattering por MATCHING MODAL EXACTO (no CAP).

Formulacion: operador FD en el interior + condiciones de contorno que imponen,
en los extremos, la forma modal asintotica exacta. Las amplitudes de los modos
salientes y evanescentes son incognitas; se resuelve un unico sistema lineal
disperso. Estable (no propaga exponenciales) e INDEPENDIENTE del transfer.

Orden de trabajo (ChatGPT):
  (3) verificar el matching en HOMOGENEO (transmision=1, reflexion=0, beta=0)
  (4) control SIN horizonte (subsonico): beta hacia norma negativa ~ 0
  (5) perfil con HORIZONTE: extraer beta, convergencia en N y L
  (6) comparar con el transfer

Incognitas: u_0..u_{N-1}, v_0..v_{N-1}, a_L (salientes izq.), a_R (salientes der.)
Ecuaciones: FD interior (filas 1..N-2) + matching modal en los 2 puntos de cada
extremo. A_in=1 (canal entrante upstream, norma +).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
from fd_bdg_operator import build_MBdG, plane_wave_mode, quartic_roots, g

n_u, mach_u, w = 1.0, 0.5, 2.0
vLb = mach_u*np.sqrt(g*n_u); jcur = n_u*vLb; mu = g*n_u + 0.5*vLb**2
omega = 0.1


def bg(x, mode):
    """mode: 'horizon' (Mach_d=2), 'subsonic' (Mach_d=0.8), 'homog' (constante)."""
    if mode == 'homog':
        return np.full_like(x, n_u), np.full_like(x, vLb), np.zeros_like(x)
    md = 2.0 if mode == 'horizon' else 0.8
    n_d = (jcur/(md*np.sqrt(g)))**(2/3)
    t = np.tanh(x/w); s = 0.5*(1+t); Dn = n_d-n_u
    n = n_u+Dn*s; n1 = Dn*(1-t**2)/(2*w); n2 = -Dn*t*(1-t**2)/w**2
    v = jcur/n; dv = -jcur*n1/n**2; Qp = 0.5*(n2/(2*n)-n1**2/(4*n**2))
    return n, v, mu-g*n-0.5*v**2+Qp


def modes_full(n, v, V):
    """4 modos (propagantes + evanescentes), normalizados a |flujo|=1 los propag."""
    c = np.sqrt(g*n); out = []
    for k in quartic_roots(omega, v, c):
        (U, Vv), lam = plane_wave_mode(k, n, v, V, mu, omega)
        prop = abs(k.imag) < 1e-7
        if prop:
            kr = k.real; Om = omega - v*kr; vg = v + (c**2*kr + kr**3/2)/Om
            J = (kr+v)*abs(U)**2 + (kr-v)*abs(Vv)**2
            U, Vv = U/np.sqrt(abs(J)), Vv/np.sqrt(abs(J))
            out.append(dict(k=k, U=U, V=Vv, prop=True, vg=vg,
                            norm=int(np.sign(Om)), Jsign=int(np.sign(J))))
        else:
            out.append(dict(k=k, U=U, V=Vv, prop=False))
    return out


def solve(L, N, mode):
    x = np.linspace(-L, L, N); dx = x[1]-x[0]
    n, v, V = bg(x, mode)
    Mop = (build_MBdG(x, n, v, V, mu) - omega*sp.identity(2*N)).tolil()
    mL = modes_full(n[0], v[0], V[0]); mR = modes_full(n[-1], v[-1], V[-1])

    Lin = next(m for m in mL if m['prop'] and m['vg'] > 0)           # entrante (+)
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0)            # saliente izq
             or (not m['prop'] and m['k'].imag < 0)]                # evan decae en -inf
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0)           # saliente der
             or (not m['prop'] and m['k'].imag > 0)]                # evan decae en +inf
    mLn, mRn = len(Laway), len(Raway)
    nT = 2*N + mLn + mRn
    iu = lambda i: i; iv = lambda i: N+i
    iaL = lambda j: 2*N+j; iaR = lambda j: 2*N+mLn+j

    A = sp.lil_matrix((nT, nT), dtype=complex); rhs = np.zeros(nT, dtype=complex)
    # interior FD (filas 1..N-2 para u y para v)
    for i in range(1, N-1):
        A.rows[iu(i)] = list(Mop.rows[iu(i)]); A.data[iu(i)] = list(Mop.data[iu(i)])
        A.rows[iv(i)] = list(Mop.rows[iv(i)]); A.data[iv(i)] = list(Mop.data[iv(i)])

    def match_left(prow, p):      # fila prow: campo en punto p = suma modal izq
        A[prow, prow] = 0.0
        for comp, idx in ((Lin['U'], iu(p)), (Lin['V'], iv(p))):
            pass
        # u_p
        A[iu(p), :] = 0; A[iu(p), iu(p)] = 1.0
        for j, m in enumerate(Laway): A[iu(p), iaL(j)] = -m['U']*np.exp(1j*m['k']*x[p])
        rhs[iu(p)] = Lin['U']*np.exp(1j*Lin['k']*x[p])
        # v_p
        A[iv(p), :] = 0; A[iv(p), iv(p)] = 1.0
        for j, m in enumerate(Laway): A[iv(p), iaL(j)] = -m['V']*np.exp(1j*m['k']*x[p])
        rhs[iv(p)] = Lin['V']*np.exp(1j*Lin['k']*x[p])

    def match_right(p):
        A[iu(p), :] = 0; A[iu(p), iu(p)] = 1.0
        for j, m in enumerate(Raway): A[iu(p), iaR(j)] = -m['U']*np.exp(1j*m['k']*x[p])
        A[iv(p), :] = 0; A[iv(p), iv(p)] = 1.0
        for j, m in enumerate(Raway): A[iv(p), iaR(j)] = -m['V']*np.exp(1j*m['k']*x[p])

    # matching en punto 0 y N-1 (reemplazan las filas de borde del FD)
    match_left(iu(0), 0); match_right(N-1)
    # 4 ecuaciones extra (segundo punto de cada extremo) -> filas 2N..2N+3
    extra = 2*N
    for (p, side) in ((1, 'L'), (N-2, 'R')):
        for comp in ('U', 'V'):
            A[extra, :] = 0
            fld = iu(p) if comp == 'U' else iv(p)
            A[extra, fld] = 1.0
            if side == 'L':
                for j, m in enumerate(Laway): A[extra, iaL(j)] = -m[comp]*np.exp(1j*m['k']*x[p])
                rhs[extra] = Lin[comp]*np.exp(1j*Lin['k']*x[p])
            else:
                for j, m in enumerate(Raway): A[extra, iaR(j)] = -m[comp]*np.exp(1j*m['k']*x[p])
            extra += 1

    z = spla.spsolve(A.tocsc(), rhs)
    aL = z[2*N:2*N+mLn]; aR = z[2*N+mLn:2*N+mLn+mRn]

    # amplitudes de canales propagantes
    refl = {m['norm']: abs(aL[j])**2 for j, m in enumerate(Laway) if m['prop']}
    tran = {m['norm']: abs(aR[j])**2 for j, m in enumerate(Raway) if m['prop']}
    beta2 = sum(abs(aR[j])**2 for j, m in enumerate(Raway) if m['prop'] and m['norm'] < 0)
    # pseudo-unitaridad: sum_out eta_NORMA |a|^2 = eta_in  (NORMA, no corriente/vg)
    flux_out = (sum(abs(aL[j])**2*m['norm'] for j, m in enumerate(Laway) if m['prop'])
                + sum(abs(aR[j])**2*m['norm'] for j, m in enumerate(Raway) if m['prop']))
    return dict(beta2=beta2, refl=refl, tran=tran, flux_out=flux_out,
                flux_in=Lin['norm'], mLn=mLn, mRn=mRn,
                psi=z[:2*N], x=x, n=n, v=v, N=N,
                aL=aL, aR=aR, Laway=Laway, Raway=Raway, Lin=Lin)


def solve_inject(L, N, mode, incL=1.0, incR=None):
    """Generalizacion: inyecta canal(es) entrante(s) con amplitudes dadas.
       incL = amplitud del entrante upstream (Lin, vg>0).
       incR = lista de amplitudes de los entrantes downstream (Rin, vg<0).
       Devuelve amplitudes salientes clasificadas por (norma, lado)."""
    x = np.linspace(-L, L, N); dx = x[1]-x[0]
    n, v, V = bg(x, mode)
    Mop = (build_MBdG(x, n, v, V, mu) - omega*sp.identity(2*N)).tolil()
    mL = modes_full(n[0], v[0], V[0]); mR = modes_full(n[-1], v[-1], V[-1])
    Lin = next(m for m in mL if m['prop'] and m['vg'] > 0)
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]          # entrantes downstream
    if incR is None: incR = [0.0]*len(Rin)
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0) or (not m['prop'] and m['k'].imag < 0)]
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0) or (not m['prop'] and m['k'].imag > 0)]
    mLn, mRn = len(Laway), len(Raway); nT = 2*N + mLn + mRn
    iu = lambda i: i; iv = lambda i: N+i; iaL = lambda j: 2*N+j; iaR = lambda j: 2*N+mLn+j
    A = sp.lil_matrix((nT, nT), dtype=complex); rhs = np.zeros(nT, dtype=complex)
    for i in range(1, N-1):
        A.rows[iu(i)] = list(Mop.rows[iu(i)]); A.data[iu(i)] = list(Mop.data[iu(i)])
        A.rows[iv(i)] = list(Mop.rows[iv(i)]); A.data[iv(i)] = list(Mop.data[iv(i)])
    def set_row(row, fld_idx, away, iax, inc_terms, p):
        A[row, :] = 0; A[row, fld_idx] = 1.0
        for j, m in enumerate(away):
            comp = m['U'] if fld_idx < N else m['V']
            A[row, iax(j)] = -comp*np.exp(1j*m['k']*x[p])
        s = 0j
        for amp, m in inc_terms:
            comp = m['U'] if fld_idx < N else m['V']
            s += amp*comp*np.exp(1j*m['k']*x[p])
        rhs[row] = s
    Lterms = [(incL, Lin)]; Rterms = list(zip(incR, Rin))
    for p in (0, 1):
        set_row(iu(p) if p == 0 else 2*N, iu(p), Laway, iaL, Lterms, p)
        set_row(iv(p) if p == 0 else 2*N+1, iv(p), Laway, iaL, Lterms, p)
    for k, p in enumerate((N-1, N-2)):
        set_row(iu(p) if p == N-1 else 2*N+2, iu(p), Raway, iaR, Rterms, p)
        set_row(iv(p) if p == N-1 else 2*N+3, iv(p), Raway, iaR, Rterms, p)
    z = spla.spsolve(A.tocsc(), rhs)
    aL = z[2*N:2*N+mLn]; aR = z[2*N+mLn:2*N+mLn+mRn]
    out = []
    for j, m in enumerate(Laway):
        if m['prop']: out.append((m['norm'], 'L', abs(aL[j])**2))
    for j, m in enumerate(Raway):
        if m['prop']: out.append((m['norm'], 'R', abs(aR[j])**2))
    return out, Rin


print("="*66); print("JUEZ INDEPENDIENTE v2 — matching modal exacto"); print("="*66)

print("\n(3) VERIFICACION HOMOGENEA (sin transicion): transmision=1, refl=0, beta=0")
r = solve(L=20, N=3000, mode='homog')
print(f"    transmision (norma+)={r['tran'].get(1,0):.4f}  reflexion={sum(r['refl'].values()):.2e}")
print(f"    beta^2={r['beta2']:.2e}  flujo: in={r['flux_in']:+d} out={r['flux_out']:+.4f}")
ok3 = abs(r['tran'].get(1,0)-1) < 1e-2 and sum(r['refl'].values()) < 1e-2
print(f"    -> {'OK: matching correcto en homogeneo' if ok3 else 'REVISAR matching'}")

print("\n(4) CONTROL SIN HORIZONTE (subsonico Mach_d=0.8): sin canal norma negativa")
r = solve(L=20, N=3000, mode='subsonic')
print(f"    canales norma- en Raway: {r['mRn']} modos; beta^2={r['beta2']:.2e}")
print(f"    conservacion flujo: in={r['flux_in']:+d} out={r['flux_out']:+.4f}  "
      f"(dif {abs(r['flux_in']-r['flux_out']):.1e})")

print("\n(5)+(6) PERFIL CON HORIZONTE: beta y convergencia")
import bdg_smatrix as B
S, ei, eo = B.smatrix(omega)[:3]
b2_tr = abs(S[np.where(eo<0)[0][0], np.where(ei>0)[0][0]])**2
print(f"    beta^2 TRANSFER = {b2_tr:.5f}")
print(f"    beta^2 MATCHING, convergencia en N (L=20):")
for N in (2000, 4000, 6000):
    r = solve(L=20, N=N, mode='horizon')
    print(f"      N={N}: beta^2={r['beta2']:.5f}  flujo_out={r['flux_out']:+.4f}")
print(f"    convergencia en L (N=4000):")
for L in (16, 20, 26, 32):
    r = solve(L=L, N=int(L*200), mode='horizon')
    print(f"      L={L}: beta^2={r['beta2']:.5f}  flujo_out={r['flux_out']:+.4f}")
