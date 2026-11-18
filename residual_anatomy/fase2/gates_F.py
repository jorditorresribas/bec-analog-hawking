"""
FASE 2 — gates tecnicos F0-F3 (criterios en docs/FASE2_DISENO.md §5, en commit antes de implementarlos).
Uso: python3 residual_anatomy/fase2/gates_F.py [F0 F1 F2 F3]   (sin argumentos: todos, en orden; si uno falla, se
informan los siguientes que no dependen de el). Salidas: out/fase2_F*.csv. No ejecuta la medicion fisica.
"""
import os, sys, io, csv, time, types, contextlib
import numpy as np, scipy.sparse as sp
from scipy.optimize import minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__)); RA = os.path.dirname(HERE); ROOT = os.path.dirname(RA)
OUT = os.path.join(ROOT, 'out'); sys.path.insert(0, HERE)
import sextic as S
with contextlib.redirect_stdout(io.StringIO()):
    import mode_matching as MM
import diag_conditioning as D
t0 = time.time()
VERD = {}


def wcsv(name, rows):
    if not rows: return
    keys = list(dict.fromkeys(k for r in rows for k in r))
    with open(os.path.join(OUT, name), 'w', newline='') as f:
        w = csv.DictWriter(f, fieldnames=keys); w.writeheader(); w.writerows(rows)


def order(a, b, f=2.0):
    return float(np.log(a/b)/np.log(f)) if a > 0 and b > 0 else float('nan')


def T0(w):
    x = np.linspace(-max(40., 8*w+20.), max(40., 8*w+20.), 80000); n, v, _, _ = S.background(x, w)
    f = np.sqrt(S.g*n) - S.JCUR/n; i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    xH = x[i] + f[i]/(f[i]-f[i+1])*(x[i+1]-x[i]); return abs(np.interp(xH, x, np.gradient(f, x)))/(2*np.pi)


OMF = np.array([0.05, 0.065, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.27])     # banda del benchmark, en unidades de T0


# ============================== F0 ==============================
def F0():
    print("F0 — fondo consistente y Goldstone (w = 2, L = 40, interior |x| <= 30)", flush=True)
    rows = []; res = {}
    for N in (2000, 4000, 8000):
        x = np.linspace(-40, 40, N); n, v, V0, n4 = S.background(x, 2.0); s = np.sqrt(n)
        Phi = np.concatenate([s, -s]).astype(complex); mk = np.abs(x) <= 30.0; mk2 = np.concatenate([mk, mk])
        for gam in (0.0, 1.0, 10.0):
            for cons in ((True, False) if gam > 0 else (True,)):
                r = S.operator(x, n, v, V0, n4, gam, cons) @ Phi
                val = float(np.linalg.norm(r[mk2])/np.linalg.norm(Phi[mk2]))
                res[(N, gam, cons)] = val
                rows.append(dict(N=N, gamma=gam, fondo='consistente' if cons else 'NO consistente (V0)', residuo_rel=val))
    ok = True
    for gam in (1.0, 10.0):
        da = max(abs(res[(N, gam, True)] - res[(N, 0.0, True)])/res[(N, 0.0, True)] for N in (2000, 4000, 8000))
        pc = [order(res[(2000, gam, True)], res[(4000, gam, True)]), order(res[(4000, gam, True)], res[(8000, gam, True)])]
        cr = min(res[(N, gam, False)]/res[(N, gam, True)] for N in (2000, 4000, 8000))
        pn = order(res[(4000, gam, False)], res[(8000, gam, False)])
        a = da <= 1e-12; b = all(1.8 <= p <= 2.2 for p in pc); c = cr >= 100 and pn < 0.5
        ok &= a and b and c
        print(f"  gamma={gam}: (a) dif. con v1.0 {da:.1e} [{'OK' if a else 'FALLA'}]  (b) ordenes {pc[0]:.2f}, {pc[1]:.2f} "
              f"[{'OK' if b else 'FALLA'}]  (c) NO consistente / consistente >= {cr:.0f}, orden {pn:.2f} [{'OK' if c else 'FALLA'}]")
    for N in (2000, 4000, 8000):
        print(f"  N={N}: residuo v1.0 {res[(N, 0.0, True)]:.3e} | gamma=10 consistente {res[(N, 10.0, True)]:.3e}"
              f" | NO consistente {res[(N, 10.0, False)]:.3e}")
    wcsv('fase2_F0.csv', rows); VERD['F0'] = ok
    print(f"  F0: {'PASA' if ok else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


# ============================== F1 ==============================
def omega_max_hom(n, v, gam):
    c = np.sqrt(S.g*n); lo, hi = 1e-4, 3.0
    for _ in range(60):
        mid = 0.5*(lo+hi); nreal = sum(abs(k.imag) < 1e-7 for k in S.roots(mid, v, c, gam))
        lo, hi = (mid, hi) if nreal >= 4 else (lo, mid)
    return lo


def F1():
    print("F1 — homogeneo (subsonico Mach 0,5 y supersonico Mach 2), gamma = 1 y 10", flush=True)
    nd = (S.JCUR/2.0)**(2/3); vd = S.JCUR/nd; Vd = S.MU - S.g*nd - 0.5*vd**2
    regs = {'sub': (1.0, 0.5, 0.0), 'sup': (nd, vd, Vd)}
    rows = []; ok = True
    for gam in (1.0, 10.0):
        om_max = omega_max_hom(nd, vd, gam)
        for reg, (n0, v0, V0) in regs.items():
            for om in (0.05, 0.2*om_max):
                ms = S.modes(om, n0, v0, V0, gam)
                nreal = sum(m['prop'] for m in ms); ncpx = len(ms) - nreal
                want = (2, 4) if reg == 'sub' else (4, 2)
                cnt = (nreal, ncpx) == want
                # (b) normas
                nrm = all(int(np.sign(abs(m['U'])**2 - abs(m['V'])**2)) == m['norm'] and
                          abs(abs(m['vg']*(abs(m['U'])**2 - abs(m['V'])**2)) - 1) <= 1e-12 for m in ms if m['prop'])
                symr = max(m['sym_res'] for m in ms)/max(abs(om), 1e-12)
                # (a) residuo FD de cada onda plana, orden en N
                fdres = {}
                for N in (2000, 4000, 8000):
                    x = np.linspace(-20, 20, N); z = np.zeros(N)
                    M = S.operator(x, np.full(N, n0), np.full(N, v0), np.full(N, V0), z, gam)
                    sl = slice(N//2-50, N//2+50); worst = 0.0
                    for m in ms:
                        e = np.exp(1j*m['k']*(x - x[N//2])); Psi = np.concatenate([m['U']*e, m['V']*e])
                        r = M @ Psi - om*Psi
                        worst = max(worst, np.max(np.abs(np.concatenate([r[:N][sl], r[N:][sl]]))) /
                                    np.max(np.abs(np.concatenate([Psi[:N][sl], Psi[N:][sl]]))))
                    fdres[N] = worst
                pfd = order(fdres[4000], fdres[8000])
                a = 1.8 <= order(fdres[2000], fdres[4000]) <= 2.2 and 1.8 <= pfd <= 2.2 and symr <= 1e-10
                # (d) transmision
                N = 8000; x = np.linspace(-20, 20, N); z = np.zeros(N)
                injs = [('L', 0)] + ([('R', 'pos'), ('R', 'neg')] if reg == 'sup' else [])
                dworst = 0.0; Tworst = 0.0
                for inj in injs:
                    r = S.scatter(x, np.full(N, n0), np.full(N, v0), np.full(N, V0), z, om, gam, inj)
                    if not r['ok']: dworst = np.inf; break
                    amps = sorted((a_ for (_, _, a_) in r['outs']), reverse=True)
                    Tworst = max(Tworst, abs(amps[0] - 1)); dworst = max(dworst, sum(amps[1:]))
                d_ = Tworst <= 1e-6 and dworst <= 1e-6
                okc = a and nrm and cnt and d_; ok &= okc
                rows.append(dict(gamma=gam, region=reg, omega=om, n_reales=nreal, n_complejas=ncpx, recuento_ok=cnt,
                                 normas_ok=nrm, sym_res=symr, fd_res_2000=fdres[2000], fd_res_4000=fdres[4000],
                                 fd_res_8000=fdres[8000], orden_fd=pfd, T_menos_1=Tworst, R_y_otros=dworst, ok=okc))
                print(f"  gamma={gam} {reg} omega={om:.4f}: raices {nreal}+{ncpx} [{'OK' if cnt else 'FALLA'}]  normas "
                      f"[{'OK' if nrm else 'FALLA'}]  FD {fdres[2000]:.1e}->{fdres[8000]:.1e} orden {pfd:.2f} "
                      f"[{'OK' if a else 'FALLA'}]  |T-1|={Tworst:.1e} resto={dworst:.1e} [{'OK' if d_ else 'FALLA'}]", flush=True)
    # (e) pseudo-unitariedad en el horizonte
    for gam in (1.0, 10.0):
        T = T0(2.0); defs = []
        for f in (1.0, 1.5, 2.25):
            dmax = max(S.nL_horizon(2.0, om, gam, f)['defect'] for om in (0.05*T, 0.15*T, 0.27*T))
            defs.append(dmax)
        pe = order(defs[0], defs[1], 1.5); pe2 = order(defs[1], defs[2], 1.5)
        e = 1.5 <= pe <= 2.5 and 1.5 <= pe2 <= 2.5 and defs[-1] <= 1e-3; ok &= e
        rows.append(dict(gamma=gam, region='horizonte w=2', defecto_N0=defs[0], defecto_1p5N0=defs[1], defecto_2p25N0=defs[2],
                         orden=pe2, ok=e))
        print(f"  (e) horizonte gamma={gam}: defecto pseudo-u {defs[0]:.1e} -> {defs[1]:.1e} -> {defs[2]:.1e}, ordenes "
              f"{pe:.2f}, {pe2:.2f} [{'OK' if e else 'FALLA'}]", flush=True)
    wcsv('fase2_F1.csv', rows); VERD['F1'] = ok
    print(f"  F1: {'PASA' if ok else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


# ============================== F2 ==============================
def nL_v10(w, om, f=1.0):
    """v1.0 (mode_matching.solve_inject), con su spsolve sustituido por el solve refinado (como diag_convergence)."""
    L = max(40., 8*w+20.); N = int(round(L*240*f))
    base = MM.spla; MM.spla = types.SimpleNamespace(spsolve=lambda A, b: D.solve(A.tocsc(), b, 'eqref'))
    try:
        MM.w = w; MM.omega = om; x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
        Rin = [q for q in MM.modes_full(n[-1], v[-1], V[-1]) if q['prop'] and q['vg'] < 0]
        with contextlib.redirect_stdout(io.StringIO()):
            out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=[1.0 if q['norm'] < 0 else 0.0 for q in Rin])
    finally:
        MM.spla = base
    return sum(a for (q, s, a) in out if q > 0 and s == 'L')


def F2():
    print("F2 — limite gamma -> 0 (w = 2, banda del benchmark, solve refinado en los dos lados)", flush=True)
    T = T0(2.0); rows = []; worst = 0.0
    for om in OMF*T:
        a = S.nL_horizon(2.0, om, 0.0)['n_L']; b = nL_v10(2.0, om)
        d = abs(a-b)/abs(b); worst = max(worst, d); rows.append(dict(parte='a', omega=om, nL_nuevo=a, nL_v10=b, dif_rel=d))
    okA = worst <= 1e-8
    print(f"  (a) gamma = 0: driver nuevo frente a v1.0, max dif rel = {worst:.2e} (tol 1e-8) [{'OK' if okA else 'FALLA'}]")
    if VERD.get('F0') is False:
        print("  (b) NO se ejecuta: usa el horizonte con gamma > 0, que depende de F0 (fallido)")
        wcsv('fase2_F2.csv', rows); VERD['F2'] = okA
        print(f"  F2: {'PASA' if okA else 'FALLA'} (solo la parte a)   [{time.time()-t0:.0f} s]\n", flush=True); return
    om = 0.10*T; n0 = S.nL_horizon(2.0, om, 0.0, 2.25)['n_L']; prev = None
    print("  (b) informativo (N = 2,25·N0 para resolver el par evanescente): n_L(gamma) - n_L(0)")
    for gam in (1e-2, 3e-3, 1e-3, 3e-4):
        d = S.nL_horizon(2.0, om, gam, 2.25)['n_L'] - n0
        p = order(abs(prev[1]), abs(d), prev[0]/gam) if prev else float('nan')
        print(f"      gamma={gam:.0e}: {d:+.3e}" + (f"   orden en gamma {p:.2f}" if prev else ''))
        rows.append(dict(parte='b', gamma=gam, omega=om, dif=d, orden=p)); prev = (gam, d)
    wcsv('fase2_F2.csv', rows); VERD['F2'] = okA
    print(f"  F2: {'PASA' if okA else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


# ============================== F3 ==============================
def F3():
    print("F3 — no degeneracion: (gamma=10, w=2) frente a la familia de reescalados n_L^(0, w/s)(s·omega)", flush=True)
    w = 2.0; T = T0(w); oms = OMF*T; rows = []
    ng = np.array([S.nL_horizon(w, om, 10.0)['n_L'] for om in oms])
    ng15 = np.array([S.nL_horizon(w, om, 10.0, 1.5)['n_L'] for om in oms])
    noise = float(np.max(np.abs(ng - ng15)/ng))
    fam = lambda s, om_arr: np.array([S.nL_horizon(w/s, s*om, 0.0)['n_L'] for om in om_arr])
    Dfun = lambda s, target, om_arr: float(np.max(np.abs(target - fam(s, om_arr))/target))
    ss = np.linspace(0.5, 2.5, 21); Ds = [Dfun(s, ng, oms) for s in ss]
    j = int(np.argmin(Ds)); r = minimize_scalar(lambda s: Dfun(s, ng, oms), bounds=(ss[max(j-1, 0)], ss[min(j+1, 20)]),
                                                method='bounded', options=dict(xatol=1e-4))
    Dmin, smin = float(r.fun), float(r.x)
    okA = Dmin >= 100*noise
    print(f"  gamma=10: ruido (N frente a 1,5N) = {noise:.2e}; D_min = {Dmin:.3e} en s = {smin:.4f}; "
          f"D_min/ruido = {Dmin/noise:.0f} (criterio >= 100) [{'OK' if okA else 'FALLA'}]", flush=True)
    for s, d in zip(ss, Ds): rows.append(dict(caso='gamma10', s=s, D=d))
    rows.append(dict(caso='gamma10_min', s=smin, D=Dmin, ruido=noise))
    # control positivo: lambda = 4 (fase 1), desde out/diag_T2_nL.csv (QP_LAM)
    t2 = list(csv.DictReader(open(os.path.join(OUT, 'diag_T2_nL.csv'))))
    get = lambda N: {float(q['omega']): float(q['n_L']) for q in t2 if q['variante'] == 'QP_LAM' and float(q['lam']) == 4.0 and q['N'] == N}
    a9, a14 = get('9600'), get('14400'); om4 = np.array(sorted(a9)); nl4 = np.array([a9[o] for o in om4])
    noise4 = float(np.max(np.abs(nl4 - np.array([a14[o] for o in om4]))/nl4))
    r4 = minimize_scalar(lambda s: Dfun(s, nl4, om4), bounds=(1.8, 2.2), method='bounded', options=dict(xatol=1e-5))
    okC = r4.fun <= 10*noise4
    print(f"  control lambda=4 (degenerado por construccion): ruido = {noise4:.2e}; D_min = {r4.fun:.3e} en s = {r4.x:.4f} "
          f"(esperado s = 2); D_min/ruido = {r4.fun/noise4:.1f} (criterio <= 10) [{'OK' if okC else 'FALLA'}]", flush=True)
    rows.append(dict(caso='control_lambda4', s=float(r4.x), D=float(r4.fun), ruido=noise4))
    wcsv('fase2_F3.csv', rows); VERD['F3'] = okA and okC
    print(f"  F3: {'PASA' if VERD['F3'] else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


if __name__ == '__main__':
    todo = sys.argv[1:] or ['F0', 'F1', 'F2', 'F3']
    for gname in todo:
        if gname == 'F3' and VERD.get('F0') is False:
            print("F3: NO se ejecuta (usa el horizonte con gamma > 0, que depende de F0, fallido)\n"); VERD['F3'] = None; continue
        globals()[gname]()
    print("RESUMEN: " + '  '.join(f"{k} {'NO EJECUTADO' if v is None else ('PASA' if v else 'FALLA')}" for k, v in VERD.items())
          + f"   [{time.time()-t0:.0f} s]")
    sys.exit(0 if all(v is True for v in VERD.values()) else 1)
