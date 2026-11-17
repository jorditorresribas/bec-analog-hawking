"""
DIAGNOSTICO (no gate) — suelo de redondeo del solver de matching modal y efecto del equilibrado.
Rama diag-condicionamiento (2026-11-17). NO modifica el solver ni los gates: re-ensambla el MISMO sistema
lineal que gate_E1_lrcp.solve_d2 (replica literal, verificada bit a bit contra solve_d2: control R0) y lo
resuelve de varias formas.

Solvers comparados (mismo sistema A z = rhs):
  base  : scipy.sparse.linalg.spsolve (lo que usan v1.0 y E1)
  eq    : equilibrado de Ruiz (filas y columnas, 5 iteraciones, norma infinito) + spsolve
  ref   : base + refinamiento iterativo en precision mixta (residuo r = rhs - A z en long double de 80 bits,
          eps ~1e-19; correccion con la misma LU en double; 4 iteraciones)
  eqref : eq + refinamiento
'ref' y 'eqref' aproximan la solucion EXACTA del sistema almacenado en double -> referencia para medir el
error de redondeo del solve (err_solve = |n_L(solver) - n_L(ref)| / n_L(ref)).

Sensibilidad a los DATOS (independiente del solver): n_L con L -> L*(1+2e-16) (1 ulp), resuelto con 'ref'.
Si esa sensibilidad es grande incluso con 'ref', el sistema esta mal condicionado como problema (no se
arregla con mejor algebra lineal); si solo es grande con 'base', es redondeo del solve y el equilibrado o
el refinamiento lo quitan.

Configuraciones: (E1) waterfall de Larre et al. 2012, m_u=0.5, L=30, dx en (0.02, 0.01, 0.005), ventana W1
(12 omegas), ajuste cubico de omega*n_L -> a, Richardson en dx como en el gate; (T) perfil tanh, el caso
mas exigente de diag_scaling_equivalence (w=1, L=20, N=14400, 9 omegas).
Salida: out/diag_conditioning.csv (por punto) y out/diag_conditioning_a.csv (coeficiente a de E1).
"""
import os, sys, io, csv, time, contextlib
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver')); sys.path.insert(0, os.path.join(HERE, 'gates'))
with contextlib.redirect_stdout(io.StringIO()):
    import gate_E1_lrcp as E                      # main() no se ejecuta al importar
    import fd_bdg_operator as F
    import correlations as C
g = F.g


# ---------------- ensamblado (replica literal de gate_E1_lrcp.solve_d2) ----------------
def assemble(x, n, v, V, mu, omega):
    N = len(x)
    Mop = (F.build_MBdG(x, n, v, V, mu) - omega*sp.identity(2*N)).tolil()
    mL = E.modes(n[0], v[0], V[0], mu, omega); mR = E.modes(n[-1], v[-1], V[-1], mu, omega)
    Lin = [m for m in mL if m['prop'] and m['vg'] > 0]
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0) or (not m['prop'] and m['k'].imag < 0)]
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0) or (not m['prop'] and m['k'].imag > 0)]
    assert len(Lin) == 1 and len(Rin) == 2 and len(Laway) == 2 and len(Raway) == 2
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
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
    Lterms = [(0.0, Lin[0])]; Rterms = list(zip(incR, Rin))
    for p in (0, 1):
        set_row(iu(p) if p == 0 else 2*N, iu(p), Laway, iaL, Lterms, p)
        set_row(iv(p) if p == 0 else 2*N+1, iv(p), Laway, iaL, Lterms, p)
    for p in (N-1, N-2):
        set_row(iu(p) if p == N-1 else 2*N+2, iu(p), Raway, iaR, Rterms, p)
        set_row(iv(p) if p == N-1 else 2*N+3, iv(p), Raway, iaR, Rterms, p)
    def nL_of(z):
        aL = z[2*N:2*N+mLn]
        return sum(abs(aL[j])**2 for j, m in enumerate(Laway) if m['prop'] and m['norm'] > 0)
    return A.tocsc(), rhs, nL_of


# ---------------- solvers ----------------
def ruiz(A, iters=5):
    A = A.tocsr().astype(complex); m, n = A.shape
    dr = np.ones(m); dc = np.ones(n)
    for _ in range(iters):
        ab = abs(A)
        r = np.sqrt(ab.max(axis=1).toarray().ravel()); c = np.sqrt(ab.max(axis=0).toarray().ravel())
        r[r == 0] = 1; c[c == 0] = 1
        A = sp.diags(1/r) @ A @ sp.diags(1/c); dr /= r; dc /= c
    return A.tocsc(), dr, dc


def resid_ld(A, z, rhs):
    """r = rhs - A z con productos y sumas en long double (80 bits)."""
    Ac = A.tocsr(); zl = z.astype(np.clongdouble)
    prod = Ac.data.astype(np.clongdouble)*zl[Ac.indices]
    rows = np.repeat(np.arange(Ac.shape[0]), np.diff(Ac.indptr))
    Az = np.zeros(Ac.shape[0], dtype=np.clongdouble); np.add.at(Az, rows, prod)
    return (rhs.astype(np.clongdouble) - Az)


def solve(A, rhs, kind):
    if kind == 'base':
        return spla.spsolve(A, rhs)
    if kind in ('eq', 'eqref'):
        As, dr, dc = ruiz(A); bs = dr*rhs
        if kind == 'eq':
            return dc*spla.spsolve(As, bs)
        lu = spla.splu(As); y = lu.solve(bs)
        for _ in range(4):
            r = resid_ld(As, y, bs); y = (y.astype(np.clongdouble) + lu.solve(r.astype(complex))).astype(complex)
        return dc*y
    if kind == 'ref':
        lu = spla.splu(A); z = lu.solve(rhs)
        for _ in range(4):
            r = resid_ld(A, z, rhs); z = (z.astype(np.clongdouble) + lu.solve(r.astype(complex))).astype(complex)
        return z
    raise ValueError(kind)


KINDS = ('base', 'eq', 'ref', 'eqref')


def point(builder, Lval, om):
    """n_L por solver y n_L(ref) con L perturbado 1 ulp."""
    x, n, v, V, mu = builder(Lval)
    A, rhs, nL_of = assemble(x, n, v, V, mu, om)
    res = {k: nL_of(solve(A, rhs, k)) for k in KINDS}
    x2, n2, v2, V2, mu2 = builder(Lval*(1+2e-16))
    A2, rhs2, nL2 = assemble(x2, n2, v2, V2, mu2, om)
    res['base_ulp'] = nL2(solve(A2, rhs2, 'base'))
    res['ref_ulp'] = nL2(solve(A2, rhs2, 'ref'))
    return res


def main():
    t0 = time.time(); rows = []
    # ---- R0: la replica de ensamblado + 'base' reproduce gate_E1_lrcp.solve_d2 bit a bit ----
    m = 0.5; a_th, b_th, T57 = E.lrcp_theory(m)
    x, n, v, V, mu = E.bg_lrcp(E.grid(30.0, 0.02), 30.0, m); om = T57*0.1
    nl_gate, _, _ = E.solve_d2(x, n, v, V, mu, om)
    A, rhs, nL_of = assemble(x, n, v, V, mu, om)
    R0 = abs(nL_of(solve(A, rhs, 'base')) - nl_gate)
    print(f"R0 replica vs gate solve_d2: |dn_L| = {R0:.1e}  -> {'OK' if R0 == 0 else 'STOP'}")
    if R0 != 0:
        return 2
    # residuo relativo de cada solver (backward error) en ese punto
    for k in KINDS:
        z = solve(A, rhs, k); r = resid_ld(A, z, rhs)
        print(f"   backward error {k:5s}: ||rhs-Az||/||rhs|| = {float(np.abs(r).max()/np.abs(rhs).max()):.1e}")

    # ---- (E1) waterfall LRCP ----
    oms = T57*np.geomspace(*E.W1, E.NW)
    a_rows = []
    for dx in E.DXS:
        N = E.grid(30.0, dx)
        builder = lambda Lv, N=N: E.bg_lrcp(N, Lv, m)
        per = {k: [] for k in KINDS + ('base_ulp', 'ref_ulp')}
        for om in oms:
            r = point(builder, 30.0, om)
            for k in per: per[k].append(r[k])
            rows.append(dict(config='E1_m0.5', dx=dx, N=N, omega=om, **{f'nL_{k}': r[k] for k in per}))
        per = {k: np.array(vv) for k, vv in per.items()}
        ref = per['ref']
        msg = '  '.join(f"{k}:{np.max(np.abs(per[k]-ref)/ref):.1e}" for k in ('base', 'eq', 'eqref'))
        print(f"E1 dx={dx}: err_solve max vs ref -> {msg} | ulp(L): base {np.max(np.abs(per['base_ulp']-per['base'])/per['base']):.1e}"
              f"  ref {np.max(np.abs(per['ref_ulp']-ref)/ref):.1e}   [{time.time()-t0:.0f} s]", flush=True)
        for k in KINDS + ('base_ulp', 'ref_ulp'):
            a, b = E.fit_ab(oms, per[k], 3)
            a_rows.append(dict(dx=dx, solver=k, a=a, b=b, a_rel_err=(a-a_th)/a_th))
    # Richardson en dx de a (como el gate) por solver
    print("\nE1 coeficiente a (W1, cubico) — a_th = %.9f" % a_th)
    summ = []
    for k in KINDS + ('base_ulp', 'ref_ulp'):
        seq = [next(q['a'] for q in a_rows if q['dx'] == dx and q['solver'] == k) for dx in E.DXS]
        d1, d2 = abs(seq[0]-seq[1]), abs(seq[1]-seq[2]); p = np.log2(d1/d2)
        rich = seq[2] + (seq[2]-seq[1])/(2**p-1)
        rich2 = seq[2] + (seq[2]-seq[1])/3.0
        print(f"  {k:8s}: a(dx)=" + ', '.join(f'{s:.9f}' for s in seq) +
              f"  p={p:.3f}  Richardson(p obs)={rich:.9f} err {abs(rich-a_th)/a_th:.2e}  (p=2: err {abs(rich2-a_th)/a_th:.2e})")
        summ.append(dict(dx='richardson', solver=k, a=rich, b=np.nan, a_rel_err=(rich-a_th)/a_th, p_obs=p,
                         a_rich_p2=rich2, a_rel_err_p2=(rich2-a_th)/a_th))
    # ---- (T) tanh, caso mas exigente ----
    lam = 4.0; s = np.sqrt(lam); Lt, Nt = 40.0/s, 14400
    ref_T2 = [q for q in csv.DictReader(open(os.path.join(OUT, 'diag_T2_nL.csv')))
              if q['variante'] == 'QP_LAM' and float(q['lam']) == lam and int(q['N']) == Nt]
    def builder_t(Lv):
        xx = np.linspace(-Lv, Lv, Nt)
        with contextlib.redirect_stdout(io.StringIO()):
            nn, vv, VV = C.bg(xx, 2.0/s, 'horizon')
        return xx, nn, vv, VV, C.mu
    per = {k: [] for k in KINDS + ('base_ulp', 'ref_ulp', 'T2')}
    for q in ref_T2:
        om = s*float(q['omega']); r = point(builder_t, Lt, om)
        for k in KINDS + ('base_ulp', 'ref_ulp'): per[k].append(r[k])
        per['T2'].append(float(q['n_L']))
        rows.append(dict(config='tanh_w1_N14400', dx=2*Lt/(Nt-1), N=Nt, omega=om, **{f'nL_{k}': r[k] for k in KINDS + ('base_ulp', 'ref_ulp')}))
    per = {k: np.array(vv) for k, vv in per.items()}; ref = per['ref']
    print(f"\nTANH w=1 N=14400: err_solve max vs ref -> " +
          '  '.join(f"{k}:{np.max(np.abs(per[k]-ref)/ref):.1e}" for k in ('base', 'eq', 'eqref')) +
          f" | ulp(L): base {np.max(np.abs(per['base_ulp']-per['base'])/per['base']):.1e}  ref {np.max(np.abs(per['ref_ulp']-ref)/ref):.1e}"
          f" | equivalencia lam=4: T2 vs base {np.max(np.abs(per['T2']-per['base'])/per['base']):.1e}, T2 vs ref {np.max(np.abs(per['T2']-ref)/ref):.1e}")
    with open(os.path.join(OUT, 'diag_conditioning.csv'), 'w', newline='') as f:
        wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
    with open(os.path.join(OUT, 'diag_conditioning_a.csv'), 'w', newline='') as f:
        keys = ['dx', 'solver', 'a', 'b', 'a_rel_err', 'p_obs', 'a_rich_p2', 'a_rel_err_p2']
        wr = csv.DictWriter(f, fieldnames=keys); wr.writeheader()
        for q in a_rows + summ: wr.writerow({k: q.get(k, '') for k in keys})
    print(f"\nfin ({time.time()-t0:.0f} s). CSV: out/diag_conditioning.csv, out/diag_conditioning_a.csv")
    return 0


if __name__ == '__main__':
    sys.exit(main())
