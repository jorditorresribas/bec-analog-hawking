"""
FASE 2 — formulacion revisada (docs/FASE2_DISENO.md §6). Modelo fenomenologico de dispersion U(k) = g(1 + gamma k^4):
(omega - v k)^2 = c^2 k^2 + k^4/4 + gamma c^2 k^6.

Ensamblado en (s, d) = (u + v, u - v):   fila 1:  -iB s + A d = omega s ;   fila 2:  (A + 2K) s - iB d = omega d
  A = -1/2 D2 + a0,  a0 = 1/2 v^2 + g n + V0 - mu  (= Q_p con el fondo consistente V_gamma = V0 - gamma g n'''')
  B = v D1 + 1/2 v'  ,  K = g sqrt(n) (1 + gamma D4) sqrt(n)
Diferencias finitas de 4o orden. Matching con los modos EXACTOS de la red (raices del polinomio caracteristico discreto).
No modifica v1.0 ni sextic.py (que se conserva como registro de F0-F3).
"""
import os, sys, io, contextlib
import numpy as np, scipy.sparse as sp, mpmath

HERE = os.path.dirname(os.path.abspath(__file__)); RA = os.path.dirname(HERE); ROOT = os.path.dirname(RA)
for p in (os.path.join(ROOT, 'solver'), RA, HERE):
    if p not in sys.path: sys.path.insert(0, p)
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C
    import diag_conditioning as Dg
g = C.g; MU = C.mu; JCUR = C.jcur; N_U = C.n_u

# esténciles de 4o orden: {desplazamiento: coeficiente} (sin el factor de h)
S1 = {-2: 1/12, -1: -8/12, 1: 8/12, 2: -1/12}
S2 = {-2: -1/12, -1: 16/12, 0: -30/12, 1: 16/12, 2: -1/12}
S4 = {-3: -1/6, -2: 12/6, -1: -39/6, 0: 56/6, 1: -39/6, 2: 12/6, 3: -1/6}


def band(st, N, h, p):
    return sp.diags([np.full(N - abs(o), c) for o, c in st.items()], list(st.keys()), (N, N), format='csr')/h**p


# ---------------- fondo (perfil tanh de v1.0, derivadas analiticas) ----------------
def background(x, w):
    n, v, V0 = C.bg(x, w, 'horizon')
    n_d = (JCUR/(2.0*np.sqrt(g)))**(2/3); Dn = n_d - N_U
    t = np.tanh(x/w); n1 = Dn*(1 - t**2)/(2*w)
    vp = -JCUR*n1/n**2
    n4 = -4.0*Dn*t*(3*t**2 - 2)*(1 - t**2)/w**4
    return dict(n=n, v=v, vp=vp, V0=V0, n4=n4)


def homog(x, n0, v0):
    N = len(x); V = MU - g*n0 - 0.5*v0**2
    return dict(n=np.full(N, n0), v=np.full(N, v0), vp=np.zeros(N), V0=np.full(N, V), n4=np.zeros(N))


def blocks(x, bg, gamma, consistent=True):
    N = len(x); h = x[1] - x[0]
    a0 = 0.5*bg['v']**2 + g*bg['n'] + bg['V0'] - MU + (0 if consistent else gamma*g*bg['n4'])
    A = -0.5*band(S2, N, h, 2) + sp.diags(a0)
    B = sp.diags(bg['v']) @ band(S1, N, h, 1) + sp.diags(0.5*bg['vp'])
    sq = sp.diags(np.sqrt(bg['n']))
    K = g*(sq @ (sp.identity(N) + (gamma*band(S4, N, h, 4) if gamma else 0*sp.identity(N))) @ sq)
    return A.tocsr(), B.tocsr(), K.tocsr(), a0


def operator(x, bg, gamma, consistent=True):
    """M en (s, d): [[-iB, A], [A + 2K, -iB]] (2N x 2N)."""
    A, B, K, _ = blocks(x, bg, gamma, consistent)
    return sp.bmat([[-1j*B, A], [A + 2*K, -1j*B]], format='csr')


# ---------------- modos exactos de la red ----------------
def lpoly(st, scale=1):
    return {o: c*scale for o, c in st.items()}


def lmul(a, b):
    out = {}
    for i, x in a.items():
        for j, y in b.items():
            out[i+j] = out.get(i+j, 0) + x*y
    return out


def ladd(*ps):
    out = {}
    for p in ps:
        for i, x in p.items(): out[i] = out.get(i, 0) + x
    return out


def leval(p, z): return sum(c*z**o for o, c in p.items())
def lderiv(p): return {o-1: c*o for o, c in p.items() if o != 0}


MP = mpmath.mp
# Coeficientes de los esténciles como RACIONALES EXACTOS (numerador, denominador); se convierten a mpf dentro de
# workdps. (Ejecucion 1 de F1': se crearon como mpf al importar, con 15 digitos, y el redondeo, multiplicado por
# gamma g/h^4 ~ 1e7, desplazaba las raices de la red como h^-4.)
RS1 = {-2: (1, 12), -1: (-8, 12), 1: (8, 12), 2: (-1, 12)}
RS2 = {-2: (-1, 12), -1: (16, 12), 0: (-30, 12), 1: (16, 12), 2: (-1, 12)}
RS4 = {-3: (-1, 6), -2: (12, 6), -1: (-39, 6), 0: (56, 6), 1: (-39, 6), 2: (12, 6), 3: (-1, 6)}


def mpst(rs):
    return {o: MP.mpf(a)/b for o, (a, b) in rs.items()}


DPS = 60


def symbols(n0, v0, a0, gamma, h):
    """Simbolos de la red en alta precision (los datos double se toman como exactos)."""
    n0, v0, a0, gamma, h = (MP.mpf(float(q)) for q in (n0, v0, a0, gamma, h))
    At = ladd(lpoly(mpst(RS2), -MP.mpf(1)/2/h**2), {0: a0})
    Bt = lpoly(mpst(RS1), -1j*v0/h)                                 # -iB(z), con v' = 0 en la region asintotica
    Kt = ladd({0: g*n0}, lpoly(mpst(RS4), gamma*g*n0/h**4) if gamma else {})
    return At, Bt, Kt


def lattice_modes(omega, n0, v0, a0, gamma, h, side, prop_tol=1e-20):
    """Raices z de F(z) = (B~ - w)^2 - A~(A~ + 2K~) = 0 y vectores (S, D), todo con 60 digitos; despues a double."""
    with mpmath.workdps(DPS):
        At, Bt, Kt = symbols(n0, v0, a0, gamma, h); om = MP.mpf(float(omega)); hh = MP.mpf(float(h)); vv = MP.mpf(float(v0))
        Bw = ladd(Bt, {0: -om})
        F = ladd(lmul(Bw, Bw), lpoly(lmul(At, ladd(At, lpoly(Kt, 2))), -1))
        lo, hi = min(F), max(F)
        coeffs = [MP.mpc(F.get(o, 0)) for o in range(hi, lo-1, -1)]
        zs = mpmath.polyroots(coeffs, maxsteps=500, extraprec=400)
        Fz = lderiv(F); out = []
        for z in zs:
            Mz = mpmath.matrix([[leval(Bw, z), leval(At, z)], [leval(At, z) + 2*leval(Kt, z), leval(Bw, z)]])
            # vector nulo: de la fila con mayor norma
            r0 = [Mz[0, 0], Mz[0, 1]] if abs(Mz[0, 0]) + abs(Mz[0, 1]) >= abs(Mz[1, 0]) + abs(Mz[1, 1]) else [Mz[1, 0], Mz[1, 1]]
            Sv, Dv = r0[1], -r0[0]; nr = mpmath.sqrt(abs(Sv)**2 + abs(Dv)**2); Sv, Dv = Sv/nr, Dv/nr
            res = max(abs(Mz[0, 0]*Sv + Mz[0, 1]*Dv), abs(Mz[1, 0]*Sv + Mz[1, 1]*Dv))
            k = -1j*mpmath.log(z)/hh
            m = dict(z=complex(z), k=complex(k), S=complex(Sv), D=complex(Dv), null_res=float(res))
            if abs(abs(z) - 1) < prop_tol:
                kr = mpmath.re(k)
                Fw = -2*leval(Bw, z); Fk = leval(Fz, z)*1j*hh*z
                vg = mpmath.re(-Fk/Fw); dens = mpmath.re(mpmath.conj(Sv)*Dv); J = vg*dens
                sc = 1/mpmath.sqrt(abs(J))
                m.update(k=float(kr), S=complex(Sv*sc), D=complex(Dv*sc), prop=True, vg=float(vg), norm=int(mpmath.sign(dens)),
                         Om=float(om - vv*kr), z=complex(mpmath.exp(1j*kr*hh)))
            else:
                m.update(prop=False)
            out.append(m)
    return out


def classify(omega, bgL, bgR, gamma, h):
    mL = lattice_modes(omega, *bgL, gamma, h, 'L'); mR = lattice_modes(omega, *bgR, gamma, h, 'R')
    Lin = [m for m in mL if m['prop'] and m['vg'] > 0]
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0) or (not m['prop'] and abs(m['z']) > 1)]
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0) or (not m['prop'] and abs(m['z']) < 1)]
    return Lin, Rin, Laway, Raway, mL, mR


# ---------------- scattering ----------------
def scatter(x, bg, omega, gamma, inject, kind='eqref', consistent=True):
    N = len(x); h = x[1] - x[0]
    A, B, K, a0 = blocks(x, bg, gamma, consistent)
    rs = 3 if gamma else 2; rd = 2; npt = rs + rd
    bgL = (bg['n'][0], bg['v'][0], a0[0]); bgR = (bg['n'][-1], bg['v'][-1], a0[-1])
    Lin, Rin, Laway, Raway, _, _ = classify(omega, bgL, bgR, gamma, h)
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
    I = sp.identity(N, format='csr')
    M = sp.bmat([[-1j*B - omega*I, A], [A + 2*K, -1j*B - omega*I]], format='csr')
    nT = 2*N + 2*npt; offL, offR = 2*N, 2*N + npt
    keep = np.zeros(2*N, bool); keep[2:N-2] = True; keep[N+rs:2*N-rs] = True
    rows_int = np.where(keep)[0]
    Mi = M[rows_int, :].tocoo()
    r_ = list(rows_int[Mi.row]); c_ = list(Mi.col); d_ = list(Mi.data)       # indices de fila ORIGINALES
    rhs = np.zeros(nT, dtype=complex)
    free = sorted(set(range(nT)) - set(rows_int))
    # matching: s en 2*rs nodos y d en 2*rd nodos por lado
    eqs = []
    for p in range(2*rs): eqs.append(('L', 's', p))
    for p in range(2*rd): eqs.append(('L', 'd', p))
    for p in range(2*rs): eqs.append(('R', 's', N-1-p))
    for p in range(2*rd): eqs.append(('R', 'd', N-1-p))
    assert len(free) == len(eqs), (len(free), len(eqs))
    for row, (side, fld, p) in zip(free, eqs):
        away, off, inc = (Laway, offL, list(zip(incL, Lin))) if side == 'L' else (Raway, offR, list(zip(incR, Rin)))
        pref = 0 if side == 'L' else N-1
        col = p if fld == 's' else N + p
        r_.append(row); c_.append(col); d_.append(1.0)
        for j, m in enumerate(away):
            r_.append(row); c_.append(off + j); d_.append(-(m['S'] if fld == 's' else m['D'])*m['z']**(p - pref))
        rhs[row] = sum(a*(m['S'] if fld == 's' else m['D'])*m['z']**(p - pref) for a, m in inc)
    Afull = sp.csc_matrix((d_, (r_, c_)), shape=(nT, nT))
    if kind == 'debug':
        return dict(A=Afull, rhs=rhs, Laway=Laway, Raway=Raway, Lin=Lin, Rin=Rin)
    z = Dg.solve(Afull, rhs, kind)
    aL = z[offL:offL+npt]; aR = z[offR:offR+npt]
    outs = ([(m['norm'], 'L', abs(aL[j])**2) for j, m in enumerate(Laway) if m['prop']] +
            [(m['norm'], 'R', abs(aR[j])**2) for j, m in enumerate(Raway) if m['prop']])
    nL = sum(a for (q, s_, a) in outs if q > 0 and s_ == 'L')
    defect = abs(sum(q*a for (q, s_, a) in outs) - eta_in)
    return dict(ok=True, chan=chan, outs=outs, n_L=nL, defect=defect, N=N)


def grid(w, h, Lfac=1.0):
    L = max(40.0, 8*w + 20.0)*Lfac; N = int(round(2*L/h)) + 1
    return np.linspace(-L, L, N)


def nL_horizon(w, omega, gamma, h, kind='eqref', Lfac=1.0):
    x = grid(w, h, Lfac); bg = background(x, w)
    return scatter(x, bg, omega, gamma, ('R', 'neg'), kind)
