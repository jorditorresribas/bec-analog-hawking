"""
FASE 2 — modelo con dispersion de dos escalas (docs/FASE2_DISENO.md).
Interaccion de alcance finito U(k) = g(1 + gamma k^4)  =>  (omega - v k)^2 = c^2 k^2 + k^4/4 + gamma c^2 k^6.

NO modifica v1.0: reutiliza fd_bdg_operator.build_MBdG (operador v1.0), correlations.bg (perfil tanh de v1.0) y
diag_conditioning.solve (equilibrado + refinamiento). Todo lo nuevo esta aqui.

Operador (fondo consistente, V_ext = V_gamma = V0 - gamma g n''''):
    M_gamma = build_MBdG(x, n, v, V0, mu) + gamma g [[P, P], [-P, -P]],   P = diag(sqrt n) D4 diag(sqrt n)
Con fondo NO consistente (V_ext = V0) la diagonal lleva ademas + gamma g n'''' (control negativo de F0).
"""
import os, sys, io, contextlib
import numpy as np, scipy.sparse as sp

HERE = os.path.dirname(os.path.abspath(__file__)); RA = os.path.dirname(HERE); ROOT = os.path.dirname(RA)
for p in (os.path.join(ROOT, 'solver'), RA):
    if p not in sys.path: sys.path.insert(0, p)
with contextlib.redirect_stdout(io.StringIO()):
    import fd_bdg_operator as F
    import correlations as C
    import diag_conditioning as D
g = F.g; MU = C.mu; JCUR = C.jcur; N_U = C.n_u


# ---------------- fondo ----------------
def background(x, w):
    """Perfil tanh de v1.0 (Mach 0.5 -> 2, anchura w): n, v, V0 (potencial de v1.0) y n'''' analitica."""
    n, v, V0 = C.bg(x, w, 'horizon')
    n_d = (JCUR/(2.0*np.sqrt(g)))**(2/3); Dn = n_d - N_U
    t = np.tanh(x/w)
    n4 = -4.0*Dn*t*(3*t**2 - 2)*(1 - t**2)/w**4
    return n, v, V0, n4


def D4(N, dx):
    e = np.ones(N)
    return sp.diags([e[:-2], -4*e[:-1], 6*e, -4*e[:-1], e[:-2]], [-2, -1, 0, 1, 2], (N, N))/dx**4


def operator(x, n, v, V0, n4, gamma, consistent=True):
    """M_gamma (2N x 2N, csr). consistent=False: V_ext = V0 (fondo NO estacionario con la interaccion nueva)."""
    N = len(x); dx = x[1] - x[0]
    Vdiag = V0 if consistent else V0 + gamma*g*n4      # build_MBdG usa W = 2gn + V - mu; Hartree aporta +gamma g n''''
    M = F.build_MBdG(x, n, v, Vdiag, MU)
    if gamma == 0:
        return M
    s = sp.diags(np.sqrt(n)); P = (s @ D4(N, dx) @ s)*(gamma*g)
    return (M + sp.bmat([[P, P], [-P, -P]], format='csr')).tocsr()


# ---------------- modos homogeneos ----------------
def roots(omega, v, c, gamma):
    if gamma == 0:
        return np.roots([-0.25, 0.0, v**2 - c**2, -2*omega*v, omega**2])
    return np.roots([-gamma*c**2, 0.0, -0.25, 0.0, v**2 - c**2, -2*omega*v, omega**2])


def symbol(k, n, v, V, gamma):
    W = 2*g*n + V - MU; Gk = g*n*(1 + gamma*k**4)
    Lp = 0.5*(k + v)**2 + W + g*n*gamma*k**4
    Lm = 0.5*(k - v)**2 + W + g*n*gamma*k**4
    return np.array([[Lp, Gk], [-Gk, -Lm]], dtype=complex)


def modes(omega, n, v, V, gamma, prop_tol=1e-7):
    c = np.sqrt(g*n); out = []
    for k in roots(omega, v, c, gamma):
        vals, vecs = np.linalg.eig(symbol(k, n, v, V, gamma))
        i = int(np.argmin(np.abs(vals - omega))); U, Vv = vecs[:, i]
        if abs(k.imag) < prop_tol:
            kr = k.real; Om = omega - v*kr
            vg = v + (c**2*kr + kr**3/2 + 3*gamma*c**2*kr**5)/Om
            # simbolo en kr real (mismo eigvec que v1.0 con gamma=0)
            vals, vecs = np.linalg.eig(symbol(kr, n, v, V, gamma)); i = int(np.argmin(np.abs(vals - omega))); U, Vv = vecs[:, i]
            J = vg*(abs(U)**2 - abs(Vv)**2); s = 1/np.sqrt(abs(J))
            out.append(dict(k=kr, U=U*s, V=Vv*s, prop=True, vg=vg, norm=int(np.sign(Om)), sym_res=abs(vals[i]-omega)))
        else:
            out.append(dict(k=k, U=U, V=Vv, prop=False, sym_res=abs(vals[i]-omega)))
    return out


def classify(omega, nL, vL, VL, nR, vR, VR, gamma):
    mL = modes(omega, nL, vL, VL, gamma); mR = modes(omega, nR, vR, VR, gamma)
    Lin = [m for m in mL if m['prop'] and m['vg'] > 0]
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0) or (not m['prop'] and m['k'].imag < 0)]
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0) or (not m['prop'] and m['k'].imag > 0)]
    return Lin, Rin, Laway, Raway


# ---------------- scattering ----------------
def scatter(x, n, v, V0, n4, omega, gamma, inject, kind='eqref', consistent=True):
    """inject: ('L', j) o ('R', 'pos'|'neg'). Devuelve dict con amplitudes salientes, n_L, defecto y canales."""
    N = len(x); deg = 4 if gamma == 0 else 6; npt = deg//2
    Mop = (operator(x, n, v, V0, n4, gamma, consistent) - omega*sp.identity(2*N)).tolil()
    Lin, Rin, Laway, Raway = classify(omega, n[0], v[0], V0[0], n[-1], v[-1], V0[-1], gamma)
    chan = dict(Lin=len(Lin), Rin=len(Rin), Rin_neg=sum(m['norm'] < 0 for m in Rin), nLa=len(Laway), nRa=len(Raway),
                Lout=sum(m['prop'] for m in Laway), Rout=sum(m['prop'] for m in Raway))
    if len(Laway) != npt or len(Raway) != npt:
        return dict(ok=False, chan=chan)
    incL = [0.0]*len(Lin); incR = [0.0]*len(Rin)
    if inject[0] == 'L':
        incL[inject[1]] = 1.0; eta_in = Lin[inject[1]]['norm']
    else:
        want = -1 if inject[1] == 'neg' else +1
        cand = [j for j, m in enumerate(Rin) if m['norm'] == want]
        if not cand: return dict(ok=False, chan=chan)
        incR[cand[0]] = 1.0; eta_in = want
    mLn, mRn = len(Laway), len(Raway); nT = 2*N + mLn + mRn
    iu = lambda i: i; iv = lambda i: N + i
    A = sp.lil_matrix((nT, nT), dtype=complex); rhs = np.zeros(nT, dtype=complex)
    def copy_row(dst, src):
        A.rows[dst] = list(Mop.rows[src]); A.data[dst] = list(Mop.data[src])
    lo, hi = (1, N-2) if deg == 4 else (2, N-3)
    for i in range(lo, hi+1):
        copy_row(iu(i), iu(i)); copy_row(iv(i), iv(i))
    if deg == 6:      # fila suma (u + v): el bloque P se cancela, alcanza solo +-1
        for i in (1, N-2):
            r = (Mop[iu(i), :] + Mop[iv(i), :]).tolil()
            A.rows[iu(i)] = list(r.rows[0]); A.data[iu(i)] = list(r.data[0])
    xs = x
    def set_row(row, fld, p, away, off, inc):
        """fld: 'u', 'v', 's' (= u + v) o 'd' (= u - v). Fase de cada modo referida a SU borde (x_ref = x[0] o
        x[-1]); para modos propagantes solo cambia la fase de la amplitud, no |a|^2."""
        cu, cv = {'u': (1, 0), 'v': (0, 1), 's': (1, 1), 'd': (1, -1)}[fld]
        xref = xs[0] if p < N//2 else xs[-1]
        A.rows[row] = []; A.data[row] = []
        if cu: A[row, iu(p)] = cu
        if cv: A[row, iv(p)] = cv
        for j, m in enumerate(away):
            A[row, off + j] = -(cu*m['U'] + cv*m['V'])*np.exp(1j*m['k']*(xs[p] - xref))
        rhs[row] = sum(a*(cu*m['U'] + cv*m['V'])*np.exp(1j*m['k']*(xs[p] - xref)) for a, m in inc)
    offL, offR = 2*N, 2*N + mLn
    Linc = list(zip(incL, Lin)); Rinc = list(zip(incR, Rin))
    if deg == 4:      # disposicion identica a v1.0 / gate_E1_lrcp.solve_d2
        slots = {('L', 0, 'u'): iu(0), ('L', 0, 'v'): iv(0), ('L', 1, 'u'): 2*N, ('L', 1, 'v'): 2*N+1,
                 ('R', N-1, 'u'): iu(N-1), ('R', N-1, 'v'): iv(N-1), ('R', N-2, 'u'): 2*N+2, ('R', N-2, 'v'): 2*N+3}
    else:
        # El estado del sistema discreto de orden 6 es s = u+v en 4 nodos y d = u-v en 2 (la fila u alcanza +-2 en s y
        # +-1 en d; la fila suma, +-1 en ambos). Matching de s en p = 0..3 y de d en p = 0, 1 por lado: asi ninguna fila
        # usa solo nodos fijados por el matching (lo que dejaria el sistema casi singular).
        slots = {('L', 0, 's'): iu(0), ('L', 0, 'd'): iv(0), ('L', 1, 's'): iv(1), ('L', 1, 'd'): 2*N,
                 ('L', 2, 's'): 2*N+1, ('L', 3, 's'): 2*N+2,
                 ('R', N-1, 's'): iu(N-1), ('R', N-1, 'd'): iv(N-1), ('R', N-2, 's'): iv(N-2), ('R', N-2, 'd'): 2*N+3,
                 ('R', N-3, 's'): 2*N+4, ('R', N-4, 's'): 2*N+5}
    for (side, p, fld), row in slots.items():
        if side == 'L': set_row(row, fld, p, Laway, offL, Linc)
        else: set_row(row, fld, p, Raway, offR, Rinc)
    if kind == 'debug':
        return dict(A=A.tocsc(), rhs=rhs, Laway=Laway, Raway=Raway, Lin=Lin, Rin=Rin, offL=offL, offR=offR)
    z = D.solve(A.tocsc(), rhs, kind)
    aL = z[offL:offL+mLn]; aR = z[offR:offR+mRn]
    outs = ([(m['norm'], 'L', abs(aL[j])**2) for j, m in enumerate(Laway) if m['prop']] +
            [(m['norm'], 'R', abs(aR[j])**2) for j, m in enumerate(Raway) if m['prop']])
    nL = sum(a for (q, s_, a) in outs if q > 0 and s_ == 'L')
    defect = abs(sum(q*a for (q, s_, a) in outs) - eta_in)
    return dict(ok=True, chan=chan, outs=outs, n_L=nL, defect=defect, z=z, N=N)


def grid(w, f=1.0):
    L = max(40.0, 8*w + 20.0); N = int(round(L*240*f))
    return np.linspace(-L, L, N), L, N


def nL_horizon(w, omega, gamma, f=1.0, kind='eqref'):
    x, L, N = grid(w, f); n, v, V0, n4 = background(x, w)
    r = scatter(x, n, v, V0, n4, omega, gamma, ('R', 'neg'), kind)
    return r
