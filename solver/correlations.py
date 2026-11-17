"""
RUMBO #1 — Correlacion densidad-densidad de Hawking (la "lengua" de Steinhauer).

Observable: C(x,x') = <delta n(x) delta n(x')>_connected en el vacio-in a T=0.
Formula (Recati-Pavloff-Carusotto, PRA 80 043603 (2009), arXiv:0907.4305, Ec.32):

   n G0^(2)(omega,x,x') = (1/4pi)[ sum_{I=u,d1} v*_I(x) r_I(x')
                                   + u*_d2(x) r_d2(x') + c.c. ]
   r_I = u_I + v_I ,   {u,d1,d2} = 3 canales ENTRANTES (d2 = norma NEGATIVA).

Los modos entrantes se normalizan a |J|=1 (corriente de norma) -> absorben el
1/4pi (modo_paper = modo_nuestro/sqrt(4pi)). C debe ser REAL y SIMETRICA en x<->x'
(las densidades conmutan a tiempo igual) -> simetrizamos explicitamente.

Reusa el operador BdG FD validado (fd_bdg_operator.build_MBdG) y la logica de
matching modal exacto de mode_matching.py (JUEZ independiente, validado). Aqui se
PARAMETRIZA en omega y se devuelve el CAMPO COMPLETO (u(x),v(x)) de cada canal.

knob trans-planckiano (lam): dispersion (omega-vk)^2 = c^2 k^2 + lam*k^4/4, con
c^2=gn FIJO (metrica acustica y kappa intactas). lam=1 = fisico; lam<1 hidrodinamico;
lam>1 mas dispersivo (UV). Implementado con kinetic eps_k=sqrt(lam) k^2/2 y
acoplo G=gn/sqrt(lam) (ver build_op_lambda) -> conveccion omega-vk intacta.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import os, sys, contextlib
with contextlib.redirect_stdout(open(os.devnull, 'w')):   # silencia auto-test import
    from fd_bdg_operator import g, build_MBdG, plane_wave_mode, quartic_roots

n_u, mach_u = 1.0, 0.5
vLb = mach_u*np.sqrt(g*n_u); jcur = n_u*vLb; mu = g*n_u + 0.5*vLb**2


# --------------------------------------------------------------------------
# Fondo waterfall (tanh), identico a mode_matching.bg
# --------------------------------------------------------------------------
def bg(x, w=2.0, mode='horizon'):
    if mode == 'homog':
        return np.full_like(x, n_u), np.full_like(x, vLb), np.zeros_like(x)
    md = 2.0 if mode == 'horizon' else 0.8
    n_d = (jcur/(md*np.sqrt(g)))**(2/3)
    t = np.tanh(x/w); s = 0.5*(1+t); Dn = n_d-n_u
    n = n_u+Dn*s; n1 = Dn*(1-t**2)/(2*w); n2 = -Dn*t*(1-t**2)/w**2
    v = jcur/n; Qp = 0.5*(n2/(2*n)-n1**2/(4*n**2))
    return n, v, mu-g*n-0.5*v**2+Qp


# --------------------------------------------------------------------------
# Dispersion y modos plano con knob lam
#   (omega - v k)^2 = c^2 k^2 + lam k^4 / 4 ,  c^2 = g n
# --------------------------------------------------------------------------
def quartic_roots_lam(omega, v, c, lam=1.0):
    # lam/4 k^4 + (c^2 - v^2) k^2 + 2 omega v k - omega^2 = 0
    return np.roots([lam/4.0, 0.0, c**2 - v**2, 2*omega*v, -omega**2])


def _symbol_lam(k, n, v, lam, omega):
    """Simbolo 2x2 del operador modificado (homogeneo):
       L+_sym = t k^2 + v k + G ,  L-_sym = t k^2 - v k + G ,
       t=sqrt(lam)/2 , G=g n/sqrt(lam).  M=[[L+,G],[-G,-L-]]."""
    t = np.sqrt(lam)/2.0; G = g*n/np.sqrt(lam)
    Lp = t*k*k + v*k + G
    Lm = t*k*k - v*k + G
    return np.array([[Lp, G], [-G, -Lm]], dtype=complex)


def plane_mode_lam(k, n, v, lam=1.0, omega=0.0):
    """(U,V) = autovector del SIMBOLO del operador modificado para la raiz k,
    con autovalor ~ omega (consistente con build_op_lambda). Normaliza |U^2-V^2|=1
    (signo del modo lo da el cociente); la normalizacion a |J|=1 se hace fuera."""
    M = _symbol_lam(k, n, v, lam, omega)
    vals, vecs = np.linalg.eig(M)
    i = int(np.argmin(np.abs(vals - omega)))
    U, V = vecs[:, i]
    d = U*U - V*V
    s = 1.0/np.sqrt(d) if abs(d) > 1e-14 else 1.0
    return np.array([U*s, V*s], dtype=complex)


# --------------------------------------------------------------------------
# Operador BdG FD con knob lam  (eps_k = sqrt(lam) k^2/2 -> kinetic coef sqrt(lam)/... )
#   row u:  [ -sqrt(lam)/2 (d+iv)^2 + W_lam ] u + G v = omega u
#   row v:  -G u + [ sqrt(lam)/2 (d-iv)^2 - W_lam ] v = omega v
#   con eps = -sqrt(lam)/2 d^2 (kinetic), conveccion -i v d (coef 1, NO escalada),
#   G = g n / sqrt(lam),  W_lam = 2 G + V_ext - mu.
# Chequeo: homogeneo da (omega - v k)^2 = c^2 k^2 + lam k^4/4 con c^2 = g n.
# --------------------------------------------------------------------------
def build_op_lambda(x, n, v, V, lam=1.0):
    """Operador BdG con dispersion modificada (omega-vk)^2 = c^2 k^2 + lam k^4/4,
    c^2 = g n FIJO (metrica acustica, kappa, T_H intactas). Solo cambia la
    estructura UV (k^4); el fondo fisico (incl. presion cuantica Qp) se CONSERVA.
       t = sqrt(lam)/2 (cinetico) ,  G = g n / sqrt(lam) (acoplo) ,  Omega = omega-vk.
    diagonal = G + Qp  con  Qp = V - mu + g n + v^2/2 (presion cuantica del fondo,
    ->0 asintoticamente). En lam=1: t=1/2, G=gn, diag=gn+Qp = build_MBdG EXACTO
    (reproduce beta^2=0.01008). Conveccion simetrizada coef 1 (Omega=omega-vk)."""
    N = len(x); dx = x[1]-x[0]
    t = np.sqrt(lam)/2.0
    G = g*n/np.sqrt(lam)
    Qp = V - mu + g*n + 0.5*v**2                   # presion cuantica del fondo (lam-indep)
    dv = np.gradient(v, x)
    e = np.ones(N)
    D1 = sp.diags([-e[:-1], e[:-1]], [-1, 1], (N, N)).tolil()/(2*dx)
    D2 = sp.diags([e[:-1], -2*e, e[:-1]], [-1, 0, 1], (N, N)).tolil()/dx**2
    conv = sp.diags(v)@D1 + 0.5*sp.diags(dv)       # (v d + v'/2) simetrizado
    diag = sp.diags(G + Qp)
    Gd = sp.diags(G)
    Lp = -t*D2 - 1j*conv + diag
    Lm = -t*D2 + 1j*conv + diag
    M = sp.bmat([[Lp, Gd], [-Gd, -Lm]], format='csr')
    return M


def _classify(n, v, V, c, omega, lam=1.0, op='validated'):
    """4 modos homogeneos clasificados. op='validated' usa los primitivos VALIDADOS
    (solo lam=1); op='lambda' usa el operador modificado (cualquier lam).
    vg con dispersion lam: vg = v + (c^2 k + lam k^3/2)/Omega."""
    out = []
    roots = quartic_roots(omega, v, c) if op == 'validated' \
        else quartic_roots_lam(omega, v, c, lam)
    for k in roots:
        if op == 'validated':
            (U, Vv), _ = plane_wave_mode(k, n, v, V, mu, omega)
        else:
            U, Vv = plane_mode_lam(k, n, v, lam, omega)
        prop = abs(k.imag) < 1e-7
        if prop:
            kr = k.real; Om = omega - v*kr
            vg = v + (c**2*kr + lam*kr**3/2)/Om
            # corriente de norma con dispersion lam: J=(sqrt(lam)k+v)|U|^2+(sqrt(lam)k-v)|V|^2
            sl = np.sqrt(lam)
            J = (sl*kr+v)*abs(U)**2 + (sl*kr-v)*abs(Vv)**2
            s = 1.0/np.sqrt(abs(J))
            out.append(dict(k=kr, U=U*s, V=Vv*s, prop=True, vg=vg,
                            norm=int(np.sign(Om))))
        else:
            out.append(dict(k=k, U=U, V=Vv, prop=False))
    return out


# --------------------------------------------------------------------------
# Scattering: inyecta UN canal entrante y devuelve el CAMPO COMPLETO (u(x),v(x))
#   channel: 'u'  = entrante upstream (vg>0, norma+)
#            'd1' = entrante downstream norma +  (vg<0)
#            'd2' = entrante downstream norma -  (vg<0)
# --------------------------------------------------------------------------
def scatter_field(L, N, omega, channel, w=2.0, lam=1.0, mode='horizon', op=None):
    x = np.linspace(-L, L, N)
    n, v, V = bg(x, w, mode)
    cL = np.sqrt(g*n[0]); cR = np.sqrt(g*n[-1])
    if op is None:
        op = 'validated' if lam == 1.0 else 'lambda'
    if op == 'validated':
        Mop = (build_MBdG(x, n, v, V, mu) - omega*sp.identity(2*N)).tolil()
    else:
        Mop = (build_op_lambda(x, n, v, V, lam) - omega*sp.identity(2*N)).tolil()
    mL = _classify(n[0], v[0], V[0], cL, omega, lam, op)
    mR = _classify(n[-1], v[-1], V[-1], cR, omega, lam, op)

    Lin = next(m for m in mL if m['prop'] and m['vg'] > 0)         # entrante upstream
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0)
             or (not m['prop'] and m['k'].imag < 0)]
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0)
             or (not m['prop'] and m['k'].imag > 0)]
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]             # entrantes downstream
    # clasificar Rin por norma
    Rin_pos = next((m for m in Rin if m['norm'] > 0), None)
    Rin_neg = next((m for m in Rin if m['norm'] < 0), None)

    # si el canal pedido no existe (fuera de ventana de Hawking), devolver vacio
    if (channel == 'd1' and Rin_pos is None) or (channel == 'd2' and Rin_neg is None):
        return dict(x=x, n=n, v=v, u=np.zeros(N, complex), v_mode=np.zeros(N, complex),
                    Lin=Lin, Rin_pos=Rin_pos, Rin_neg=Rin_neg, nchan=len(Rin),
                    omega=omega, out_amp={}, flux_out=0.0, eta_in=0)

    # amplitudes entrantes
    incL = 0.0; incR = [0.0]*len(Rin)
    if channel == 'u':
        incL = 1.0
    elif channel == 'd1':
        incR[Rin.index(Rin_pos)] = 1.0
    elif channel == 'd2':
        incR[Rin.index(Rin_neg)] = 1.0
    else:
        raise ValueError(channel)

    mLn, mRn = len(Laway), len(Raway); nT = 2*N + mLn + mRn
    iu = lambda i: i; iv = lambda i: N+i
    iaL = lambda j: 2*N+j; iaR = lambda j: 2*N+mLn+j
    A = sp.lil_matrix((nT, nT), dtype=complex); rhs = np.zeros(nT, dtype=complex)
    for i in range(1, N-1):
        A.rows[iu(i)] = list(Mop.rows[iu(i)]); A.data[iu(i)] = list(Mop.data[iu(i)])
        A.rows[iv(i)] = list(Mop.rows[iv(i)]); A.data[iv(i)] = list(Mop.data[iv(i)])

    Lterms = [(incL, Lin)]
    Rterms = list(zip(incR, Rin))

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

    for p in (0, 1):
        set_row(iu(p) if p == 0 else 2*N, iu(p), Laway, iaL, Lterms, p)
        set_row(iv(p) if p == 0 else 2*N+1, iv(p), Laway, iaL, Lterms, p)
    for p in (N-1, N-2):
        set_row(iu(p) if p == N-1 else 2*N+2, iu(p), Raway, iaR, Rterms, p)
        set_row(iv(p) if p == N-1 else 2*N+3, iv(p), Raway, iaR, Rterms, p)

    z = spla.spsolve(A.tocsc(), rhs)
    u_f = z[:N]; v_f = z[N:2*N]
    aL = z[2*N:2*N+mLn]; aR = z[2*N+mLn:2*N+mLn+mRn]
    # amplitudes salientes propagantes clasificadas por norma
    out_amp = {}
    for j, m in enumerate(Laway):
        if m['prop']: out_amp[('L', m['norm'])] = aL[j]
    for j, m in enumerate(Raway):
        if m['prop']: out_amp[('R', m['norm'])] = aR[j]
    # pseudo-unitaridad: sum_out eta|a|^2 debe = eta_in del canal inyectado
    flux_out = (sum(abs(aL[j])**2*m['norm'] for j, m in enumerate(Laway) if m['prop'])
                + sum(abs(aR[j])**2*m['norm'] for j, m in enumerate(Raway) if m['prop']))
    eta_in = (Lin['norm'] if channel == 'u'
              else (Rin_pos['norm'] if channel == 'd1' else Rin_neg['norm']))
    return dict(x=x, n=n, v=v, u=u_f, v_mode=v_f,
                Lin=Lin, Rin_pos=Rin_pos, Rin_neg=Rin_neg,
                nchan=len(Rin), omega=omega, out_amp=out_amp,
                flux_out=flux_out, eta_in=eta_in)


# --------------------------------------------------------------------------
# Correlacion densidad-densidad (Ec. 32), integrada en omega, simetrizada.
# --------------------------------------------------------------------------
def density_correlation(L, N, omegas, w=2.0, lam=1.0, mode='horizon',
                        ncoarse=320, terms='all', taper=True, op=None, verbose=False):
    """Devuelve xc (subgrid), C(xc,xc') simetrica real, y n(xc).
    terms: 'all' (Ec.32 completa), 'd2' (solo canal anomalo = Hawking puro),
           'ud1' (canales ordinarios norma+). taper: ventana suave en bordes de
           banda (reduce ringing de Gibbs del corte en omega)."""
    xfull = np.linspace(-L, L, N)
    idx = np.linspace(0, N-1, ncoarse).astype(int)
    xc = xfull[idx]
    n_full, _, _ = bg(xfull, w, mode)
    nc = n_full[idx]
    acc = np.zeros((ncoarse, ncoarse), dtype=complex)
    dω = omegas[1]-omegas[0]
    # ventana Hann suave en los extremos de la banda
    if taper:
        win = 0.5*(1-np.cos(np.pi*np.minimum(np.arange(len(omegas)),
                   len(omegas)-1-np.arange(len(omegas)))/(0.15*len(omegas))))
        win = np.clip(win, 0, 1)
    else:
        win = np.ones(len(omegas))
    for iw, om in enumerate(omegas):
        modes = {}
        ok = True
        for ch in ('u', 'd1', 'd2'):
            r = scatter_field(L, N, om, ch, w, lam, mode, op)
            if r['Rin_neg'] is None or r['nchan'] != 2:
                ok = False; break
            u = r['u'][idx]; vv = r['v_mode'][idx]
            modes[ch] = dict(u=u, v=vv, r=u+vv)
        if not ok:
            if verbose: print(f"  omega={om:.3f}: fuera de ventana - skip")
            continue
        br = np.zeros((ncoarse, ncoarse), dtype=complex)
        if terms in ('all', 'ud1'):
            for I in ('u', 'd1'):
                br += np.outer(np.conj(modes[I]['v']), modes[I]['r'])
        if terms in ('all', 'd2'):
            br += np.outer(np.conj(modes['d2']['u']), modes['d2']['r'])
        br = br + np.conj(br.T)               # "+ c.c." simetrizado (x<->x')
        acc += br*dω*win[iw]
        if verbose:
            print(f"  omega={om:.3f}: ok  |bracket|_max={np.abs(br).max():.3e}")
    pref = np.sqrt(np.outer(nc, nc))/(4*np.pi)
    C = np.real(pref*acc)
    return xc, C, nc
