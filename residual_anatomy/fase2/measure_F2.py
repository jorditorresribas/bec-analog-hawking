"""
FASE 2 — MEDICION FISICA pre-registrada (residual_anatomy/FASE2_PREREG.md + FASE2_PREREG_ADDENDUM.md).
NO EJECUTAR antes de: (1) revision del revisor externo; (2) gates F0'-F3' en verde (hecho: b2a9f1e); (3) este script en
commit. El script comprueba (3) por si mismo (I0) y se niega a ejecutar si no se cumple.

Medicion (fijada en el pre-registro):
  w in {2, 3, 4, 6, 8}, gamma in {0, 1, 10}; T_up = T_0(n_L) con el estimador del benchmark (banda [0.05, 0.27]*kappa/2pi,
  9 omegas, T_w = w/ln(1 + 1/n_L), ajuste cuadratico, termino constante) en h in {0.04, 0.03, 0.02} y Richardson p = 4;
  convergencia exigida |T(0.02) - T_R|/T_R <= 1e-4 (si no, STOP sin veredicto).
  ell*(w, gamma): anchura de la ventana simetrica centrada en el horizonte (la de FP_d_sensitivity.kbar, anchura TOTAL
  d_xi; la semianchura es ell*/2 y el exponente no cambia) con kbar(ell*)/2pi = T_up; se toma la PRIMERA solucion al
  crecer ell desde 0. Si no hay solucion: "sin d*".
  p(gamma): -pendiente de log ell* frente a log kappa_local (minimos cuadrados, 5 anchuras); bootstrap 1000, semilla 0.
Criterios: K (control), P1 (confirma), P0 (refuta), resto INCONCLUSO; ver FASE2_PREREG.md §4.
Uso: python3 residual_anatomy/fase2/measure_F2.py [--selftest]
Exit: 0 terminado con veredicto | 2 STOP (I0, convergencia, K o d*) | 3 rehusa sobrescribir | 1 excepcion.
"""
import os, sys, io, csv, json, time, hashlib, platform, subprocess, contextlib, traceback
import numpy as np, scipy

HERE = os.path.dirname(os.path.abspath(__file__)); RA = os.path.dirname(HERE); ROOT = os.path.dirname(RA)
OUT = os.path.join(ROOT, 'out'); sys.path.insert(0, HERE)
import sextic4 as S4

WS = (2.0, 3.0, 4.0, 6.0, 8.0); GAMMAS = (0.0, 1.0, 10.0); HS = (0.04, 0.03, 0.02)
OMF = np.array([0.05, 0.065, 0.08, 0.10, 0.12, 0.15, 0.18, 0.22, 0.27])
TOL = dict(conv=1e-4, K_T=1e-4, K_p=(0.25, 0.42), P1_dev=0.05, P1_sep=0.07, P0=0.03, min_dstar=4)
TRACKED = [os.path.join(RA, 'FASE2_PREREG.md'), os.path.join(RA, 'FASE2_PREREG_ADDENDUM.md'),
           os.path.abspath(__file__), os.path.join(HERE, 'sextic4.py')]
NAMES = ['fase2_medida_T.csv', 'fase2_medida_dstar.csv', 'fase2_medida_resultado.csv', 'fase2_medida_manifest.json']
g, n_u, jc = 1.0, 1.0, 0.5; n_d = (jc/(2.0*np.sqrt(g)))**(2/3)


def sha(p): return hashlib.sha256(open(p, 'rb').read()).hexdigest()
def git(*a):
    r = subprocess.run(['git', '-C', ROOT] + list(a), capture_output=True, text=True); return r.returncode, r.stdout.strip()


def check_I0():
    msgs = []
    for p in TRACKED:
        rel = os.path.relpath(p, ROOT)
        if git('ls-files', '--error-unmatch', rel)[0] != 0: msgs.append(f'{rel}: no esta en git'); continue
        if git('diff', '--quiet', 'HEAD', '--', rel)[0] != 0: msgs.append(f'{rel}: cambios sin commit')
    return not msgs, msgs, git('rev-parse', 'HEAD')[1]


# ---------------- geometria (identica a benchmarks/fp_v1/FP_d_sensitivity.kbar) ----------------
def geom(w):
    L = max(40., 8*w+20.); x = np.linspace(-L, L, 80000); n = n_u+(n_d-n_u)*0.5*(1+np.tanh(x/w))
    c = np.sqrt(n); v = jc/n; f = c-v; i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    xH = x[i]+f[i]/(f[i]-f[i+1])*(x[i+1]-x[i]); cH = np.interp(xH, x, c)
    kK = abs(np.interp(xH, x, np.gradient(f, x)))
    return x, f, xH, cH, kK


def kbar_width(G, ell):
    x, f, xH, _, _ = G
    return abs(np.interp(xH+0.5*ell, x, f)-np.interp(xH-0.5*ell, x, f))/ell


def ell_star(G, T):
    """primera anchura ell con kbar(ell)/2pi = T, al crecer ell desde ~0; None si no hay."""
    target = 2*np.pi*T; ells = np.geomspace(1e-3, 200.0, 4000)
    vals = np.array([kbar_width(G, e) for e in ells]) - target
    idx = np.where(np.sign(vals[:-1]) != np.sign(vals[1:]))[0]
    if len(idx) == 0: return None
    a, b = ells[idx[0]], ells[idx[0]+1]
    for _ in range(80):
        m = 0.5*(a+b)
        if np.sign(kbar_width(G, m) - target) == np.sign(kbar_width(G, a) - target): a = m
        else: b = m
    return 0.5*(a+b)


def T_up(w, gamma, h, kK):
    oms = OMF*kK/(2*np.pi); nL = []
    for om in oms:
        r = S4.nL_horizon(w, om, gamma, h)
        if not r['ok']: raise RuntimeError(f"canales: w={w} gamma={gamma} h={h} omega={om:.4e} {r['chan']}")
        nL.append(r['n_L'])
    Tw = oms/np.log(1+1/np.array(nL))
    return float(np.polyfit(oms, Tw, 2)[-1])


def fit_p(kappas, ells, seed=0, nboot=1000):
    lk, le = np.log(kappas), np.log(ells)
    p = -float(np.polyfit(lk, le, 1)[0]); rng = np.random.default_rng(seed); bs = []
    for _ in range(nboot):
        j = rng.integers(0, len(lk), len(lk))
        if len(set(j)) < 2: continue
        bs.append(-float(np.polyfit(lk[j], le[j], 1)[0]))
    return p, float(np.percentile(bs, 16)), float(np.percentile(bs, 84))


def selftest():
    """Solo geometria: kbar con la anchura de FP a d = 1.58 reproduce FP_d_sensitivity / la referencia (no calcula T_up)."""
    ref = {float(r['w']): r for r in csv.DictReader(open(os.path.join(ROOT, 'reference', 'fp_benchmark_reference.csv')))}
    worst = 0.0
    for w in WS:
        G = geom(w); _, _, _, cH, kK = G; dxi = 1.58*cH/(kK*(2*cH**2)**2)**(1/3)
        val = float(ref[w]['T_upstream'])/(kbar_width(G, dxi)/(2*np.pi))
        worst = max(worst, abs(val - float(ref[w]['Tup_over_kbar2pi'])))
    ok = worst < 1e-6
    i0, msgs, head = check_I0()
    print(f"selftest geometria: max dif con la referencia {worst:.1e} (tol 1e-6) -> {'OK' if ok else 'FALLA'}")
    print(f"I0 (pre-registro, addendum, script y sextic4 en commit): {'OK' if i0 else 'NO: ' + '; '.join(msgs)}  HEAD={head[:12]}")
    print("selftest: no se ha calculado ninguna temperatura ni ningun exponente.")
    return 0 if ok else 1


def main():
    t0 = time.time()
    if any(os.path.exists(os.path.join(OUT, f)) for f in NAMES):
        print("Rehuso sobrescribir resultados de la medicion."); return 3
    i0, msgs, head = check_I0()
    print(f"I0: {'OK' if i0 else 'NO'}  HEAD={head}")
    if not i0:
        for m in msgs: print('   ' + m)
        print("STOP: nada calculado."); return 2
    ref = {float(r['w']): float(r['Tup_over_kloc2pi']) for r in csv.DictReader(open(os.path.join(OUT, 'fp_benchmark.csv')))}
    rowsT, rowsD = [], []; T = {}; E = {}; stop = None
    for gam in GAMMAS:
        for w in WS:
            G = geom(w); kK = G[4]
            q = [T_up(w, gam, h, kK) for h in HS]
            TR = q[-1] + (q[-1]-q[-2])/((HS[-2]/HS[-1])**4 - 1); conv = abs(q[-1]-TR)/abs(TR)
            T[(gam, w)] = TR
            rowsT.append(dict(gamma=gam, w=w, kappa_local=kK, T_h004=q[0], T_h003=q[1], T_h002=q[2], T_R=TR,
                              T_R_over_kloc2pi=TR/(kK/(2*np.pi)), conv=conv))
            print(f"gamma={gam:<4} w={w}: T_up/(kappa/2pi) = {TR/(kK/(2*np.pi)):.6f}  conv {conv:.1e}   [{time.time()-t0:.0f} s]", flush=True)
            if conv > TOL['conv']: stop = stop or f'convergencia gamma={gam} w={w}'
            e = ell_star(G, TR); E[(gam, w)] = e
            rowsD.append(dict(gamma=gam, w=w, kappa_local=kK, ell_star=e if e else '', d_FP=(e*(kK*(2*G[3]**2)**2)**(1/3)/G[3]) if e else ''))
    res = []; verdict = None
    Kt = max(abs(T[(0.0, w)]/(geom(w)[4]/(2*np.pi)) - ref[w]) for w in WS)
    P = {}
    for gam in GAMMAS:
        ok_w = [w for w in WS if E[(gam, w)]]
        P[gam] = fit_p(np.array([geom(w)[4] for w in ok_w]), np.array([E[(gam, w)] for w in ok_w])) if len(ok_w) >= 2 else (np.nan, np.nan, np.nan)
        res.append(dict(gamma=gam, n_dstar=len(ok_w), p=P[gam][0], p16=P[gam][1], p84=P[gam][2]))
        print(f"p(gamma={gam}) = {P[gam][0]:.4f}  [16-84 %: {P[gam][1]:.4f}, {P[gam][2]:.4f}]  (d* en {len(ok_w)}/5 anchuras)")
    nd0 = sum(1 for w in WS if E[(0.0, w)]); nd10 = sum(1 for w in WS if E[(10.0, w)])
    K = Kt <= TOL['K_T'] and TOL['K_p'][0] <= P[0.0][0] <= TOL['K_p'][1]
    print(f"K: max |T_up/(k/2pi) - benchmark| = {Kt:.1e} (<= 1e-4); p(0) = {P[0.0][0]:.4f} en [0.25, 0.42] -> {'OK' if K else 'FALLA'}")
    if stop: verdict = f'SIN VEREDICTO ({stop})'
    elif nd0 < TOL['min_dstar'] or nd10 < TOL['min_dstar']: verdict = 'SIN VEREDICTO (d* en menos de 4 anchuras)'
    elif not K: verdict = 'SIN VEREDICTO (falla el control K)'
    else:
        p0, p10 = P[0.0][0], P[10.0][0]
        if abs(p10 - 0.2) <= TOL['P1_dev'] and p0 - p10 >= TOL['P1_sep']: verdict = 'CONFIRMA (P1)'
        elif abs(p10 - p0) <= TOL['P0']: verdict = 'REFUTA (P0)'
        else: verdict = 'INCONCLUSO'
    print(f"VEREDICTO: {verdict}")
    def wcsv(name, rows):
        with open(os.path.join(OUT, name), 'w', newline='') as f:
            wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
    wcsv(NAMES[0], rowsT); wcsv(NAMES[1], rowsD); wcsv(NAMES[2], res + [dict(gamma='K', n_dstar='', p=Kt, p16='', p84=''),
                                                                        dict(gamma='veredicto', n_dstar='', p=verdict, p16='', p84='')])
    man = dict(medicion='fase2', veredicto=verdict, git_HEAD=head, tolerancias={k: (list(v) if isinstance(v, tuple) else v) for k, v in TOL.items()},
               sha256={os.path.relpath(p, ROOT): sha(p) for p in TRACKED + [os.path.join(OUT, n) for n in NAMES[:3]]},
               python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__, segundos=round(time.time()-t0, 1),
               fecha=time.strftime('%Y-%m-%dT%H:%M:%S'))
    json.dump(man, open(os.path.join(OUT, NAMES[3]), 'w'), indent=2, default=str)
    return 0 if not verdict.startswith('SIN') else 2


if __name__ == '__main__':
    try: sys.exit(selftest() if '--selftest' in sys.argv else main())
    except Exception: traceback.print_exc(); sys.exit(1)
