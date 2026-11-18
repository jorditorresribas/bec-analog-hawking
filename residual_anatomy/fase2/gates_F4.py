"""
FASE 2 — gates revisados F0'-F3' (criterios: docs/FASE2_DISENO.md §6.4, en commit 39f7607 antes de implementar).
Formulacion: residual_anatomy/fase2/sextic4.py ((s, d), diferencias finitas de 4o orden, modos exactos de la red).
Orden: F1' (homogeneo) -> F0' -> F1'(e) (horizonte) -> F2' -> F3'. Si F1' falla, se para.
Regla del suelo de redondeo: medido como sensibilidad a 1 ulp de L, con el mismo solve; debe quedar >= 100x por debajo de
la tolerancia del gate. Salidas: out/fase2p_*.csv. No ejecuta la medicion fisica.
"""
import os, sys, io, csv, time, contextlib
import numpy as np
from scipy.optimize import brentq, minimize_scalar

HERE = os.path.dirname(os.path.abspath(__file__)); RA = os.path.dirname(HERE); ROOT = os.path.dirname(RA)
OUT = os.path.join(ROOT, 'out'); sys.path.insert(0, HERE)
import sextic4 as S4
import sextic as S2                          # solo para las raices del CONTINUO (S2.roots) y omega_max
from gates_F import nL_v10, T0, OMF, wcsv
t0 = time.time(); VERD = {}; HS = (0.04, 0.03, 0.02); ULP = 1 + 2.2e-16
nd = (S4.JCUR/2.0)**(2/3); vd = S4.JCUR/nd


def p_err(e, hs=HS):
    """orden observado a partir de errores (pares consecutivos)."""
    return [float(np.log(e[i]/e[i+1])/np.log(hs[i]/hs[i+1])) if e[i] > 0 and e[i+1] > 0 else float('nan') for i in range(len(e)-1)]


def p_three(q, hs=HS):
    """orden observado a partir de tres valores q(h) (diferencias), resolviendo (h1^p-h2^p)/(h2^p-h3^p) = (q1-q2)/(q2-q3)."""
    r = (q[0]-q[1])/(q[1]-q[2]) if q[1] != q[2] else np.inf
    f = lambda p: (hs[0]**p - hs[1]**p)/(hs[1]**p - hs[2]**p) - r
    try: return float(brentq(f, 0.2, 12))
    except Exception: return float('nan')


def rich(q, p, hs=HS):
    return q[-1] + (q[-1]-q[-2])/((hs[-2]/hs[-1])**p - 1)


def omega_max_hom(n, v, gam):
    c = np.sqrt(n); lo, hi = 1e-4, 3.0
    for _ in range(60):
        mid = 0.5*(lo+hi); nreal = sum(abs(k.imag) < 1e-7 for k in S2.roots(mid, v, c, gam))
        lo, hi = (mid, hi) if nreal >= 4 else (lo, mid)
    return lo


# ============================== F1' (homogeneo) ==============================
def F1p():
    print("F1' — homogeneo, gamma = 1 y 10, h = 0.04/0.03/0.02", flush=True)
    rows = []; ok = True; floor_worst = 0.0
    for gam in (1.0, 10.0):
        om_max = omega_max_hom(nd, vd, gam)
        for reg, (n0, v0) in (('sub', (1.0, 0.5)), ('sup', (nd, vd))):
            for om in (0.05, 0.2*om_max):
                kc = [k for k in S2.roots(om, v0, np.sqrt(n0), gam)]
                errs = []; nrm_ok = True; cnt_ok = True; Tw = 0.0; Rw = 0.0
                for h in HS:
                    ms = S4.lattice_modes(om, n0, v0, 0.0, gam, h, 'L')
                    # (a) raices fisicas de la red -> continuo
                    errs.append(max(min(abs(m['k'] - k)/abs(k) for m in ms) for k in kc))
                    # (b) normas
                    for m in ms:
                        if m['prop']:
                            dens = np.real(np.conj(m['S'])*m['D'])
                            nrm_ok &= int(np.sign(dens)) == int(np.sign(m['Om'])) and abs(abs(m['vg']*dens) - 1) <= 1e-12
                    # (c) recuento
                    x = np.linspace(-20, 20, int(round(40/h)) + 1); bg = S4.homog(x, n0, v0)
                    injs = [('L', 0)] + ([('R', 'pos'), ('R', 'neg')] if reg == 'sup' else [])
                    for inj in injs:
                        r = S4.scatter(x, bg, om, gam, inj)
                        if not r['ok'] or len(ms) != 10:
                            cnt_ok = False; continue
                        amps = sorted((a for _, _, a in r['outs']), reverse=True)
                        Tw = max(Tw, abs(amps[0]-1)); Rw = max(Rw, sum(amps[1:]))
                        if h == HS[-1]:      # (e) suelo de redondeo: 1 ulp en L
                            x2 = np.linspace(-20*ULP, 20*ULP, len(x)); r2 = S4.scatter(x2, S4.homog(x2, n0, v0), om, gam, inj)
                            a2 = sorted((a for _, _, a in r2['outs']), reverse=True)
                            floor_worst = max(floor_worst, abs(a2[0]-amps[0]))
                pa = p_err(errs)
                a = all(3.5 <= p <= 4.5 for p in pa); d = Tw <= 1e-8 and Rw <= 1e-8
                okc = a and nrm_ok and cnt_ok and d; ok &= okc
                rows.append(dict(gamma=gam, region=reg, omega=om, err_k=';'.join(f'{e:.2e}' for e in errs), orden_k=';'.join(f'{p:.2f}' for p in pa),
                                 normas=nrm_ok, recuento=cnt_ok, T_menos_1=Tw, resto=Rw, ok=okc))
                print(f"  gamma={gam} {reg} omega={om:.4f}: k red->cont {errs[0]:.1e}->{errs[-1]:.1e} ordenes {pa[0]:.2f},{pa[1]:.2f} "
                      f"[{'OK' if a else 'FALLA'}]  normas [{'OK' if nrm_ok else 'FALLA'}]  recuento [{'OK' if cnt_ok else 'FALLA'}]  "
                      f"|T-1|={Tw:.1e} resto={Rw:.1e} [{'OK' if d else 'FALLA'}]", flush=True)
    e = floor_worst <= 1e-10; ok &= e
    print(f"  (e) suelo de redondeo (1 ulp de L) en T: {floor_worst:.1e} (<= 1e-10) [{'OK' if e else 'FALLA'}]")
    rows.append(dict(gamma='todas', region='suelo', omega='', T_menos_1=floor_worst, ok=e))
    wcsv('fase2p_F1.csv', rows); VERD["F1'"] = ok
    print(f"  F1' (a-e homogeneo): {'PASA' if ok else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


# ============================== F0' ==============================
def F0p():
    print("F0' — Goldstone en (s, d): s = 0, d = 2 sqrt(n); w = 2, L = 40, interior |x| <= 30", flush=True)
    res = {}; rows = []
    for h in HS:
        x = S4.grid(2.0, h); bg = S4.background(x, 2.0); N = len(x)
        Phi = np.concatenate([np.zeros(N), 2*np.sqrt(bg['n'])]).astype(complex)
        mk = np.abs(x) <= 30.0; mk2 = np.concatenate([mk, mk])
        for gam in (0.0, 1.0, 10.0):
            for cons in ((True, False) if gam > 0 else (True,)):
                r = S4.operator(x, bg, gam, cons) @ Phi
                res[(h, gam, cons)] = float(np.linalg.norm(r[mk2])/np.linalg.norm(Phi[mk2]))
                rows.append(dict(h=h, gamma=gam, fondo='consistente' if cons else 'NO consistente', residuo=res[(h, gam, cons)]))
    da = max(abs(res[(h, g_, True)] - res[(h, 0.0, True)])/res[(h, 0.0, True)] for h in HS for g_ in (1.0, 10.0))
    pb = p_err([res[(h, 0.0, True)] for h in HS])
    okc = all(res[(HS[0], g_, False)] >= 100*res[(HS[0], g_, True)] and p_err([res[(h, g_, False)] for h in HS])[-1] < 0.5 for g_ in (1.0, 10.0))
    a = da <= 1e-14; b = all(3.5 <= p <= 4.5 for p in pb); ok = a and b and okc
    print(f"  (a) dif. entre gamma = 0, 1, 10: {da:.1e} [{'OK' if a else 'FALLA'}]   (b) residuo " +
          ' -> '.join(f"{res[(h, 0.0, True)]:.2e}" for h in HS) + f", ordenes {pb[0]:.2f}, {pb[1]:.2f} [{'OK' if b else 'FALLA'}]")
    print(f"  (c) NO consistente: " + ', '.join(f"gamma={g_}: {res[(HS[0], g_, False)]:.2e}" for g_ in (1.0, 10.0)) +
          f" [{'OK' if okc else 'FALLA'}]")
    wcsv('fase2p_F0.csv', rows); VERD["F0'"] = ok
    print(f"  F0': {'PASA' if ok else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


# ============================== F1'(e) horizonte ==============================
def F1e():
    print("F1'(e) — pseudo-unitariedad en el horizonte (w = 2), gamma = 1 y 10", flush=True)
    T = T0(2.0); rows = []; ok = True
    for gam in (1.0, 10.0):
        defs = [max(S4.nL_horizon(2.0, om, gam, h)['defect'] for om in (0.05*T, 0.15*T, 0.27*T)) for h in HS]
        pe = p_err(defs)
        fl = max(abs(S4.nL_horizon(2.0, om, gam, HS[-1], Lfac=ULP)['n_L'] - S4.nL_horizon(2.0, om, gam, HS[-1])['n_L']) /
                 S4.nL_horizon(2.0, om, gam, HS[-1])['n_L'] for om in (0.05*T, 0.15*T, 0.27*T))
        e = all(3 <= p <= 5 for p in pe) and defs[-1] <= 1e-5 and fl <= 1e-7; ok &= e
        rows.append(dict(gamma=gam, defectos=';'.join(f'{q:.2e}' for q in defs), ordenes=';'.join(f'{p:.2f}' for p in pe), suelo_nL=fl, ok=e))
        print(f"  gamma={gam}: defecto " + ' -> '.join(f'{q:.1e}' for q in defs) + f", ordenes {pe[0]:.2f}, {pe[1]:.2f}; "
              f"suelo de n_L {fl:.1e} (<= 1e-7) [{'OK' if e else 'FALLA'}]", flush=True)
    wcsv('fase2p_F1e.csv', rows); VERD["F1'(e)"] = ok
    print(f"  F1'(e): {'PASA' if ok else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


# ============================== F2' ==============================
def F2p():
    print("F2' — gamma = 0 en el limite del continuo frente a v1.0 (w = 2, 9 omegas, solve refinado)", flush=True)
    T = T0(2.0); rows = []; worst = 0.0; pw = []
    for om in OMF*T:
        q = [S4.nL_horizon(2.0, om, 0.0, h)['n_L'] for h in HS]; p = p_three(q); pw.append(p)
        l4 = rich(q, 4.0)
        v = [nL_v10(2.0, om, f) for f in (1.0, 1.5, 2.25)]; l2 = v[2] + (v[2]-v[1])/(1.5**2 - 1)
        d = abs(l4 - l2)/abs(l2); worst = max(worst, d)
        rows.append(dict(parte='a', omega=om, nuevo=';'.join(f'{z:.12f}' for z in q), orden=p, lim_nuevo=l4,
                         v10=';'.join(f'{z:.12f}' for z in v), lim_v10=l2, dif_rel=d))
    a = worst <= 1e-6 and all(3.5 <= p <= 4.5 for p in pw)
    print(f"  (a) orden observado {min(pw):.2f}-{max(pw):.2f}; |lim(nuevo) - lim(v1.0)|/lim max = {worst:.2e} (<= 1e-6) "
          f"[{'OK' if a else 'FALLA'}]", flush=True)
    om = 0.10*T; n0 = S4.nL_horizon(2.0, om, 0.0, HS[-1])['n_L']; prev = None
    print("  (b) informativo: n_L(gamma) - n_L(0) en omega = 0.1 T0, h = 0.02")
    for gam in (0.1, 0.03, 0.01):
        d = S4.nL_horizon(2.0, om, gam, HS[-1])['n_L'] - n0
        p = float(np.log(abs(prev[1])/abs(d))/np.log(prev[0]/gam)) if prev else float('nan')
        print(f"      gamma={gam}: {d:+.4e}" + (f"   orden en gamma {p:.2f}" if prev else ''))
        rows.append(dict(parte='b', gamma=gam, dif=d, orden=p)); prev = (gam, d)
    wcsv('fase2p_F2.csv', rows); VERD["F2'"] = a
    print(f"  F2': {'PASA' if a else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


# ============================== F3' ==============================
def F3p():
    print("F3' — no degeneracion (gamma = 10, w = 2) frente a la familia n_L^(0, w/s)(s omega), h = 0.02", flush=True)
    import csv as _csv
    w = 2.0; T = T0(w); oms = OMF*T; rows = []; h = HS[-1]
    ng = np.array([S4.nL_horizon(w, om, 10.0, h)['n_L'] for om in oms])
    ng3 = np.array([S4.nL_horizon(w, om, 10.0, HS[-2])['n_L'] for om in oms])
    noise = float(np.max(np.abs(ng - ng3)/ng))
    floor = float(max(abs(S4.nL_horizon(w, om, 10.0, h, Lfac=ULP)['n_L'] - n_)/n_ for om, n_ in zip(oms[:3], ng[:3])))
    fam = lambda s, om_arr: np.array([S4.nL_horizon(w/s, s*om, 0.0, h)['n_L'] for om in om_arr])
    Dfun = lambda s, target, om_arr: float(np.max(np.abs(target - fam(s, om_arr))/target))
    ss = np.linspace(0.5, 2.5, 21); Ds = [Dfun(s, ng, oms) for s in ss]; j = int(np.argmin(Ds))
    r = minimize_scalar(lambda s: Dfun(s, ng, oms), bounds=(ss[max(j-1, 0)], ss[min(j+1, 20)]), method='bounded', options=dict(xatol=1e-4))
    Dmin = float(r.fun)
    okA = Dmin >= 100*noise and floor <= noise and floor <= Dmin/100
    print(f"  gamma=10: ruido (h 0.03 frente a 0.02) {noise:.2e}; suelo {floor:.1e}; D_min = {Dmin:.3e} en s = {r.x:.4f}; "
          f"D_min/ruido = {Dmin/noise:.0f} (>= 100) [{'OK' if okA else 'FALLA'}]", flush=True)
    for s_, d_ in zip(ss, Ds): rows.append(dict(caso='gamma10', s=s_, D=d_))
    rows.append(dict(caso='gamma10_min', s=float(r.x), D=Dmin, ruido=noise, suelo=floor))
    t2 = list(_csv.DictReader(open(os.path.join(OUT, 'diag_T2_nL.csv'))))
    get = lambda N: {float(q['omega']): float(q['n_L']) for q in t2 if q['variante'] == 'QP_LAM' and float(q['lam']) == 4.0 and q['N'] == N}
    a9, a14 = get('9600'), get('14400'); om4 = np.array(sorted(a9)); nl4 = np.array([a9[o] for o in om4])
    noise4 = float(np.max(np.abs(nl4 - np.array([a14[o] for o in om4]))/nl4))
    r4 = minimize_scalar(lambda s: Dfun(s, nl4, om4), bounds=(1.8, 2.2), method='bounded', options=dict(xatol=1e-5))
    okC = r4.fun <= 10*noise4
    print(f"  control lambda=4: ruido {noise4:.2e}; D_min = {r4.fun:.3e} en s = {r4.x:.4f} (esperado 2); D_min/ruido = "
          f"{r4.fun/noise4:.1f} (<= 10) [{'OK' if okC else 'FALLA'}]", flush=True)
    rows.append(dict(caso='control_lambda4', s=float(r4.x), D=float(r4.fun), ruido=noise4))
    wcsv('fase2p_F3.csv', rows); VERD["F3'"] = okA and okC
    print(f"  F3': {'PASA' if VERD[chr(70)+chr(51)+chr(39)] else 'FALLA'}   [{time.time()-t0:.0f} s]\n", flush=True)


if __name__ == '__main__':
    F1p()
    if not VERD["F1'"]:
        print("F1' fallo: se para aqui (regla del diseño §6.4).");
    else:
        F0p(); F1e(); F2p(); F3p()
    print("RESUMEN: " + '  '.join(f"{k} {'PASA' if v else 'FALLA'}" for k, v in VERD.items()) + f"   [{time.time()-t0:.0f} s]")
    sys.exit(0 if all(VERD.values()) and len(VERD) == 5 else 1)
