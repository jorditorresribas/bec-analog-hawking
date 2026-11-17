"""
GATE EXTERNO E1 (v2: correcciones de auditoria — I3 sin extrapolacion tautologica, I4 exige salientes propagantes) — Larre, Recati, Carusotto, Pavloff, PRA 85, 013621 (2012): configuracion WATERFALL.
Diseno: residual_anatomy/GATE_EXT_LRCP_DESIGN.md (incl. I3 del referee y §7 fijado antes de implementar).

Compara a = lim omega*n_L y b (coef. lineal de omega*n_L) con las formulas analiticas A8/A14 de LRCP.
Instrumento: operador v1.0 SIN MODIFICAR (fd_bdg_operator.build_MBdG, plane_wave_mode, quartic_roots);
driver de matching modal = replica literal del algoritmo de mode_matching.solve_inject con fondo generico.
Exit: 0 terminado (con veredicto) | 2 STOP integridad (sin veredicto) | 1 excepcion | 3 rehusa sobrescribir.
"""
import os, sys, io, csv, json, time, hashlib, platform, contextlib, traceback
import numpy as np, scipy, scipy.sparse as sp, scipy.sparse.linalg as spla

HERE = os.path.dirname(os.path.abspath(__file__)); RA = os.path.dirname(HERE); ROOT = os.path.dirname(RA)
OUT = os.path.join(ROOT, 'out'); sys.path.insert(0, os.path.join(ROOT, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import fd_bdg_operator as F
    import mode_matching as MM
g = F.g

# ---------------- parametros fijados en el diseno §7 ----------------
MUS = (0.5, 0.7)
DXS = (0.02, 0.01, 0.005); L_REF = 30.0
LS = (20.0, 30.0, 40.0); DX_L = 0.01
W1 = (0.1, 0.6); W2 = (0.05, 0.3); NW = 12
TOL = dict(a_floor=1e-3, a_ceil=5e-3, b_floor=5e-3, b_ceil=2e-2, reg=1e-8, p_lo=0.8, p_hi=2.4)

def lrcp_theory(m):
    a = 2*m*(1-m)**1.5*(1+m*m)**1.5/((1+m)**0.5*(1+m+m*m)**2)                 # |f_{u,d2}|^2  (A8)
    b = 2*(-m*(2+m*m)*(1+2*m*m)/((1+m)**2*(1+m+m*m)**2))                    # 2Re(f*h)      (A14)
    T57 = 0.5*(1-m**4)**1.5/((2+m*m)*(1+2*m*m))                              # Ec. 57
    return a, b, T57

def bg_lrcp(N, L, m):
    """Waterfall LRCP en unidades hbar=m=g=n_u=1 (c_u=xi_u=1). x=0 en el nodo central (N impar)."""
    assert N % 2 == 1
    x = np.linspace(-L, L, N); i0 = (N-1)//2; x[i0] = 0.0
    ct = np.sqrt(1-m*m)
    n = np.where(np.arange(N) < i0, ct**2*np.tanh(x*ct)**2 + m*m, m*m)
    v = m/n                                                                  # j = n_u V_u = m
    U0 = m*m/2 + 1/(2*m*m) - 1                                               # Ec. 17
    V = np.where(np.arange(N) < i0, 0.0, -U0); V[i0] = -U0/2                # interfaz: promedio (§7.2)
    mu = g*1.0 + 0.5*m*m
    return x, n, v, V, mu

# ---------------- matching modal (replica de mode_matching.modes_full / solve_inject) ----------------
def modes(n0, v0, V0, mu, omega):
    c = np.sqrt(g*n0); out = []
    for k in F.quartic_roots(omega, v0, c):
        (U, Vv), _ = F.plane_wave_mode(k, n0, v0, V0, mu, omega)
        if abs(k.imag) < 1e-7:
            kr = k.real; Om = omega - v0*kr; vg = v0 + (c**2*kr + kr**3/2)/Om
            J = (kr+v0)*abs(U)**2 + (kr-v0)*abs(Vv)**2
            U, Vv = U/np.sqrt(abs(J)), Vv/np.sqrt(abs(J))
            out.append(dict(k=k, U=U, V=Vv, prop=True, vg=vg, norm=int(np.sign(Om))))
        else:
            out.append(dict(k=k, U=U, V=Vv, prop=False))
    return out

def solve_d2(x, n, v, V, mu, omega):
    """Inyecta el entrante downstream de norma negativa (d2|in) con amplitud 1. Devuelve n_L, defecto, canales."""
    N = len(x)
    Mop = (F.build_MBdG(x, n, v, V, mu) - omega*sp.identity(2*N)).tolil()
    mL = modes(n[0], v[0], V[0], mu, omega); mR = modes(n[-1], v[-1], V[-1], mu, omega)
    Lin = [m for m in mL if m['prop'] and m['vg'] > 0]
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0) or (not m['prop'] and m['k'].imag < 0)]
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0) or (not m['prop'] and m['k'].imag > 0)]
    chan = dict(Lin=len(Lin), Rin=len(Rin), Lout=sum(m['prop'] for m in Laway), Rout=sum(m['prop'] for m in Raway),
                Rin_neg=sum(m['norm'] < 0 for m in Rin))
    if not (chan['Lin'] == 1 and chan['Rin'] == 2 and chan['Rin_neg'] == 1 and len(Laway) == 2 and len(Raway) == 2
            and chan['Lout'] == 1 and chan['Rout'] == 2):            # v2: salientes PROPAGANTES (auditoria)
        return None, None, chan
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
    z = spla.spsolve(A.tocsc(), rhs)
    if not np.all(np.isfinite(z)):
        return None, None, chan
    aL = z[2*N:2*N+mLn]; aR = z[2*N+mLn:2*N+mLn+mRn]
    outs = ([(m['norm'], 'L', abs(aL[j])**2) for j, m in enumerate(Laway) if m['prop']] +
            [(m['norm'], 'R', abs(aR[j])**2) for j, m in enumerate(Raway) if m['prop']])
    nL = sum(a for (q, s_, a) in outs if q > 0 and s_ == 'L')
    defect = abs(sum(q*a for (q, s_, a) in outs) + 1.0)          # |sum eta|a|^2 - eta_in|, eta_in = -1
    return nL, defect, chan

def fit_ab(oms, nL, order):
    c = np.polyfit(oms, oms*nL, order)
    return float(c[-1]), float(c[-2])                            # a (constante), b (lineal)

def grid(L, dx):
    N = int(round(2*L/dx)) + 1
    return N if N % 2 == 1 else N+1

# ======================================================================
def main():
    t0 = time.time()
    names = ['e1_lrcp_nL.csv', 'e1_lrcp_fits.csv', 'e1_lrcp_integrity.csv', 'e1_lrcp_manifest.json']
    if any(os.path.exists(os.path.join(OUT, f)) for f in names):
        print("Rehuso sobrescribir resultados E1 existentes."); return 3
    rows_nL, rows_fit, rows_int = [], [], []
    stop = None; verdict = {}

    # ---- I1: regresion del driver contra v1.0 (solve_inject), fondo tanh, omega=0.1 ----
    MM.omega = 0.1; MM.w = 2.0
    xr = np.linspace(-40.0, 40.0, 9600)
    with contextlib.redirect_stdout(io.StringIO()):
        nr, vr, Vr = MM.bg(xr, 'horizon')
        Rin_r = [m for m in MM.modes_full(nr[-1], vr[-1], Vr[-1]) if m['prop'] and m['vg'] < 0]
        out_r, _ = MM.solve_inject(L=40.0, N=9600, mode='horizon', incL=0.0,
                                   incR=[1.0 if m['norm'] < 0 else 0.0 for m in Rin_r])
    nL_ref = sum(a for (q, s_, a) in out_r if q > 0 and s_ == 'L')
    nL_new, _, _ = solve_d2(xr, nr, vr, Vr, MM.mu, 0.1)
    I1 = abs(nL_new - nL_ref)/abs(nL_ref)
    rows_int.append(dict(prueba='I1_regresion_v1.0', m_u='', valor=I1, criterio=f"<= {TOL['reg']}", ok=I1 <= TOL['reg']))
    print(f"I1 regresion driver vs v1.0 solve_inject (omega=0.1): n_L {nL_new:.12e} vs {nL_ref:.12e}  rel {I1:.2e}")
    if not I1 <= TOL['reg']:
        stop = 'I1'

    results = {}
    for m in (MUS if stop is None else ()):
        a_th, b_th, T57 = lrcp_theory(m)
        oms1 = T57*np.geomspace(*W1, NW); oms2 = T57*np.geomspace(*W2, NW)
        print(f"\n=== m_u={m}  a_th={a_th:.6f}  b_th={b_th:.6f}  T_H={T57:.6f} ===")
        # ---- I2: residuo GP del fondo (interior), 3 mallas ----
        gp = []
        for dx in DXS:
            x, n, v, V, mu = bg_lrcp(grid(L_REF, dx), L_REF, m); s = np.sqrt(n); h = x[1]-x[0]
            lap = (s[2:] - 2*s[1:-1] + s[:-2])/h**2
            r = -0.5*lap + (V[1:-1] + g*n[1:-1] + 0.5*v[1:-1]**2 - mu)*s[1:-1]
            gp.append(float(np.max(np.abs(r))))
        p_gp = [np.log2(gp[i]/gp[i+1]) for i in range(2)]
        ok_gp = gp[0] > gp[1] > gp[2]
        rows_int.append(dict(prueba='I2_residuo_GP_max', m_u=m, valor=';'.join(f'{q:.3e}' for q in gp),
                             criterio=f"decrece; ordenes {p_gp[0]:.2f},{p_gp[1]:.2f}", ok=ok_gp))
        print(f"I2 residuo GP (max interior) dx={DXS}: {gp}  ordenes {p_gp}")
        if not ok_gp:
            stop = f'I2 m_u={m}'; break
        # ---- barridos ----
        store = {}
        def run(L, dx, tag):
            N = grid(L, dx); x, n, v, V, mu = bg_lrcp(N, L, m)
            res = {}
            for wname, oms in (('W1', oms1), ('W2', oms2)):
                nLs, dfs = [], []
                for om in oms:
                    nl, df, ch = solve_d2(x, n, v, V, mu, om)
                    if nl is None:
                        raise RuntimeError(f"I4 canales/solve fallo: m_u={m} L={L} dx={dx} omega={om} {ch}")
                    nLs.append(nl); dfs.append(df)
                    rows_nL.append(dict(m_u=m, L=L, dx=dx, N=N, ventana=wname, omega=om, omega_over_TH=om/T57,
                                        n_L=nl, defecto_pu=df))
                    rows_int.append(dict(prueba='I4_canales_por_omega', m_u=m, valor=(f"L={L} dx={dx} {wname} omega={om:.6e} "
                                         f"Lin={ch['Lin']} Rin={ch['Rin']} Rin_neg={ch['Rin_neg']} Lout={ch['Lout']} Rout={ch['Rout']}"),
                                         criterio='Lin=1 Rin=2 Rin_neg=1 Lout=1 Rout=2', ok=True))
                res[wname] = (np.array(oms), np.array(nLs), np.array(dfs))
            store[(L, dx)] = res
            return res
        try:
            for dx in DXS: run(L_REF, dx, 'dx')
            for L in LS:
                if L != L_REF: run(L, DX_L, 'L')
        except RuntimeError as e:
            rows_int.append(dict(prueba='I4_canales', m_u=m, valor=str(e), criterio='3 in / 3 out', ok=False))
            stop = f'I4 m_u={m}'; break
        rows_int.append(dict(prueba='I4_canales', m_u=m, valor='todas las omega', criterio='3 in / 3 out', ok=True))
        # ---- I3: convergencia del defecto de pseudo-unitariedad (max sobre W1) ----
        # v2 (auditoria): sin extrapolacion tautologica. Ordenes grueso (0.02->0.01) y fino (0.01->0.005) por
        # ventana; extrapolacion SOLO con p=2 declarado: D_fino - (D_medio - D_fino)/3.
        ok3 = True
        for wname in ('W1', 'W2'):
            dfm = [float(np.max(store[(L_REF, dx)][wname][2])) for dx in DXS]
            pc = np.log2(dfm[0]/dfm[1]) if dfm[1] > 0 and dfm[0] > 0 else np.nan
            pf = np.log2(dfm[1]/dfm[2]) if dfm[2] > 0 and dfm[1] > 0 else np.nan
            ext2 = dfm[2] - (dfm[1]-dfm[2])/3.0
            okw = (dfm[0] > dfm[1] > dfm[2]) and all(np.isfinite(q) and TOL['p_lo'] <= q <= TOL['p_hi'] for q in (pc, pf))
            ok3 &= bool(okw)
            rows_int.append(dict(prueba=f'I3_defecto_pu_max_{wname}', m_u=m,
                                 valor=';'.join(f'{q:.3e}' for q in dfm) + f'; p_grueso={pc:.3f}; p_fino={pf:.3f}; extrap_p2_declarado={ext2:.2e}',
                                 criterio=f"monotono y p_grueso, p_fino in [{TOL['p_lo']},{TOL['p_hi']}]", ok=bool(okw)))
            print(f"I3 {wname} defecto pu (max) dx={DXS}: {[f'{q:.3e}' for q in dfm]}  p_grueso={pc:.3f} p_fino={pf:.3f}  "
                  f"extrap(p=2 declarado)={ext2:.2e}  {'OK' if okw else 'FALLO'}")
        if not ok3:
            stop = f'I3 m_u={m}'; break
        # ---- ajustes ----
        def ab(L, dx, wname, order):
            o, nl, _ = store[(L, dx)][wname]; return fit_ab(o, nl, order)
        for (L, dx) in store:
            for wname in ('W1', 'W2'):
                for order in (2, 3):
                    a, b = ab(L, dx, wname, order)
                    rows_fit.append(dict(m_u=m, L=L, dx=dx, ventana=wname, orden=order, a=a, b=b,
                                         a_rel=(a-a_th)/a_th, b_rel=(b-b_th)/abs(b_th)))
        out = {}
        for q, idx, th in (('a', 0, a_th), ('b', 1, b_th)):
            seq = [ab(L_REF, dx, 'W1', 3)[idx] for dx in DXS]
            d1, d2 = abs(seq[0]-seq[1]), abs(seq[1]-seq[2])
            p = np.log2(d1/d2) if d2 > 0 and d1 > 0 else np.nan
            rich = seq[2] + (seq[2]-seq[1])/(2**p-1) if (np.isfinite(p) and 0.5 < p < 4) else seq[2]
            dN = abs(seq[2]-rich)
            dL = max(abs(ab(L, DX_L, 'W1', 3)[idx] - ab(L_REF, DX_L, 'W1', 3)[idx]) for L in LS)
            ref_f = ab(L_REF, DXS[-1], 'W1', 3)[idx]
            dE = max(abs(ab(L_REF, DXS[-1], w_, o_)[idx] - ref_f) for (w_, o_) in (('W1', 2), ('W2', 2), ('W2', 3)))
            delta = float(np.sqrt(dN**2 + dL**2 + dE**2))
            rel_delta = delta/abs(th); diff = rich - th; rel = abs(diff)/abs(th)
            floor, ceil = TOL[f'{q}_floor'], TOL[f'{q}_ceil']
            if rel_delta > ceil:
                st = 'INCONCLUSO'
            elif rel <= max(3*rel_delta, floor):
                st = 'PASA'
            else:
                st = 'FALLA'
            out[q] = dict(teoria=th, numerico=float(rich), seq_dx=seq, p_obs=float(p), delta_N=dN, delta_L=dL,
                          delta_est=dE, delta_rel=rel_delta, error_rel=rel, estado=st)
            print(f"  {q}: num={rich:.6f}  th={th:.6f}  err_rel={rel:.2e}  delta_rel={rel_delta:.2e} "
                  f"(N {dN:.1e}, L {dL:.1e}, est {dE:.1e}; p={p:.2f})  -> {st}")
        a_num, b_num = out['a']['numerico'], out['b']['numerico']
        out['Gamma_num'] = -2*b_num; out['TH_num'] = a_num/(-2*b_num); out['TH_ec57'] = T57
        print(f"  derivados: Gamma={-2*b_num:.5f} (th {-2*b_th:.5f})  T_H={a_num/(-2*b_num):.5f} (Ec.57 {T57:.5f})")
        results[m] = out

    # ---- clasificacion ----
    if stop is None:
        sts = [(results[m]['a']['estado'], results[m]['b']['estado']) for m in MUS]
        if all(a == 'PASA' and b == 'PASA' for a, b in sts):
            cls = 'PASS'
        elif any('FALLA' in (a, b) for a, b in sts):
            cls = 'FALLA'
        elif all(a == 'PASA' for a, b in sts):
            cls = 'PARCIAL'
        else:
            cls = 'INCONCLUSO'
    else:
        cls = None
    secs = time.time()-t0
    os.makedirs(OUT, exist_ok=True)
    def wcsv(name, rr):
        if rr:
            with open(os.path.join(OUT, name), 'w', newline='') as f:
                w_ = csv.DictWriter(f, fieldnames=list(rr[0].keys())); w_.writeheader(); w_.writerows(rr)
    wcsv(names[0], rows_nL); wcsv(names[1], rows_fit); wcsv(names[2], rows_int)
    sha = lambda p: hashlib.sha256(open(p, 'rb').read()).hexdigest()
    man = dict(gate='E1 externo LRCP waterfall', timestamp=time.strftime('%Y-%m-%d %H:%M:%S %Z'), segundos=secs,
               integridad='PASS' if stop is None else f'STOP ({stop})', clasificacion=cls,
               resultados={str(k): v for k, v in results.items()}, tolerancias=TOL,
               parametros=dict(MUS=MUS, DXS=DXS, L_REF=L_REF, LS=LS, DX_L=DX_L, W1=W1, W2=W2, NW=NW,
                               interfaz='U(0)=-U0/2', unidades='hbar=m=g=n_u=1'),
               referencia='Larre, Recati, Carusotto, Pavloff, PRA 85, 013621 (2012), arXiv:1110.4464v2, Ecs. 17, 44, 56, 57, A8, A14',
               sha256={**{p: sha(os.path.join(ROOT, p)) for p in
                           ('solver/fd_bdg_operator.py', 'solver/mode_matching.py', 'residual_anatomy/gates/gate_E1_lrcp.py',
                            'residual_anatomy/GATE_EXT_LRCP_DESIGN.md')},
                       **{f: sha(os.path.join(OUT, f)) for f in names[:3] if os.path.exists(os.path.join(OUT, f))}},
               entorno=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                            platform=platform.platform()))
    json.dump(man, open(os.path.join(OUT, names[3]), 'w'), indent=1, default=float)
    print(f"\nINTEGRIDAD: {man['integridad']}   CLASIFICACION: {cls}   ({secs:.0f} s)")
    return 0 if stop is None else 2

if __name__ == '__main__':
    try:
        sys.exit(main())
    except Exception:
        traceback.print_exc(); sys.exit(1)
