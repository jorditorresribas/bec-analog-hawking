"""
GATE E1-bis (PRE-REGISTRADO) — Larre, Recati, Carusotto, Pavloff, PRA 85, 013621 (2012): coeficiente a (A8) con
redondeo, discretizacion y extrapolacion en omega controlados.
Criterios: residual_anatomy/E1BIS_PREREG.md (deben estar en commit ANTES de ejecutar; se comprueba en I0).

Instrumento: operador v1.0 SIN MODIFICAR (fd_bdg_operator via gate_E1_lrcp); ensamblado = diag_conditioning.assemble
(replica bit a bit de gate_E1_lrcp.solve_d2); solve = diag_conditioning.solve(..., 'eqref').
Exit: 0 terminado con veredicto | 2 STOP integridad (sin veredicto) | 1 excepcion | 3 rehusa sobrescribir.
Modo --selftest: reproduce el ajuste ya publicado de E1 (m_u=0.5, W1, orden 3, dx=0.02); no toca ninguna
configuracion del gate y no escribe nada.
"""
import os, sys, io, csv, json, time, hashlib, platform, subprocess, contextlib, traceback
import numpy as np, scipy

HERE = os.path.dirname(os.path.abspath(__file__)); RA = os.path.dirname(HERE); ROOT = os.path.dirname(RA)
OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver')); sys.path.insert(0, RA); sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import gate_E1_lrcp as E
    import diag_conditioning as D
    import fd_bdg_operator as F

# ---------------- parametros fijados en E1BIS_PREREG.md §3-§5 ----------------
MU_CONF = (0.6, 0.7)              # confirmatorias (deciden el veredicto)
MU_REPL = (0.5,)                  # replica del diagnostico a posteriori (NO cuenta)
L_GATE = 30.0
DXS = (0.01, 0.005, 0.0025)
WINS = {'W2': (0.05, 0.3), 'W3': (0.025, 0.15)}; W_MAIN = 'W2'; NW = 12
ORDER = 4
SOLVER = 'eqref'
TOL = dict(I1=1e-7, p_lo=1.5, p_hi=2.5, C1=1e-6, C2=1e-6, C3=1e-7, C4=1e-3)
PREREG = os.path.join(RA, 'E1BIS_PREREG.md')
TRACKED = [PREREG, os.path.abspath(__file__), os.path.join(RA, 'diag_conditioning.py')]
NAMES = ['e1bis_nL.csv', 'e1bis_fits.csv', 'e1bis_result.csv', 'e1bis_manifest.json']


def sha(p):
    with open(p, 'rb') as f:
        return hashlib.sha256(f.read()).hexdigest()


def git(*args):
    r = subprocess.run(['git', '-C', ROOT] + list(args), capture_output=True, text=True)
    return r.returncode, r.stdout.strip()


def check_I0():
    """Los archivos de TRACKED estan en git y sin cambios respecto a HEAD."""
    msgs = []
    for p in TRACKED:
        rel = os.path.relpath(p, ROOT)
        rc, _ = git('ls-files', '--error-unmatch', rel)
        if rc != 0:
            msgs.append(f'{rel}: no esta en git'); continue
        rc, _ = git('diff', '--quiet', 'HEAD', '--', rel)
        if rc != 0:
            msgs.append(f'{rel}: cambios sin commit')
    rc, head = git('rev-parse', 'HEAD')
    _, pre_commit = git('log', '-1', '--format=%H %cI', '--', os.path.relpath(PREREG, ROOT))
    return (len(msgs) == 0), msgs, head, pre_commit


def channels(x, n, v, V, mu, om):
    """Mismo recuento de canales que gate_E1_lrcp.solve_d2 (v2: salientes propagantes)."""
    mL = E.modes(n[0], v[0], V[0], mu, om); mR = E.modes(n[-1], v[-1], V[-1], mu, om)
    Lin = [m for m in mL if m['prop'] and m['vg'] > 0]
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    Laway = [m for m in mL if (m['prop'] and m['vg'] < 0) or (not m['prop'] and m['k'].imag < 0)]
    Raway = [m for m in mR if (m['prop'] and m['vg'] > 0) or (not m['prop'] and m['k'].imag > 0)]
    ch = dict(Lin=len(Lin), Rin=len(Rin), Rin_neg=sum(m['norm'] < 0 for m in Rin),
              Lout=sum(m['prop'] for m in Laway), Rout=sum(m['prop'] for m in Raway), nLa=len(Laway), nRa=len(Raway))
    ok = (ch['Lin'] == 1 and ch['Rin'] == 2 and ch['Rin_neg'] == 1 and ch['nLa'] == 2 and ch['nRa'] == 2
          and ch['Lout'] == 1 and ch['Rout'] == 2)
    return ok, ch


def nL_point(m, dx, om, kind):
    x, n, v, V, mu = E.bg_lrcp(E.grid(L_GATE, dx), L_GATE, m)
    ok, ch = channels(x, n, v, V, mu, om)
    if not ok:
        return None, ch
    A, rhs, nL_of = D.assemble(x, n, v, V, mu, om)
    return float(nL_of(D.solve(A, rhs, kind))), ch


def selftest():
    """Reproduce out/e1_lrcp_fits.csv (m_u=0.5, L=30, dx=0.02, W1, orden 3) con 'base' y comprueba 'eqref'."""
    ref = None
    with open(os.path.join(OUT, 'e1_lrcp_fits.csv')) as f:
        for r in csv.DictReader(f):
            if (float(r['m_u']) == 0.5 and float(r['L']) == 30.0 and float(r['dx']) == 0.02
                    and r['ventana'] == 'W1' and int(r['orden']) == 3):
                ref = float(r['a'])
    _, _, T57 = E.lrcp_theory(0.5)
    oms = T57*np.geomspace(*E.W1, E.NW)
    res = {}
    for kind in ('base', 'eqref'):
        nl = []
        for om in oms:
            v_, ch = nL_point(0.5, 0.02, om, kind)
            if v_ is None:
                print(f"selftest: canales inesperados {ch}"); return 1
            nl.append(v_)
        res[kind] = E.fit_ab(oms, np.array(nl), 3)[0]
    d_base = abs(res['base'] - ref)/abs(ref); d_eq = abs(res['eqref'] - res['base'])/abs(ref)
    ok = d_base <= 1e-12 and d_eq <= 1e-7
    print(f"selftest: a(base) = {res['base']:.16f}   E1 publicado = {ref:.16f}   dif rel {d_base:.1e} (<= 1e-12)")
    print(f"          a(eqref) = {res['eqref']:.16f}   dif rel frente a base {d_eq:.1e} (<= 1e-7)")
    i0, msgs, head, pc = check_I0()
    print(f"          I0 (criterios en commit): {'OK' if i0 else 'NO: ' + '; '.join(msgs)}   HEAD={head[:12]}")
    print(f"selftest {'OK' if ok else 'FALLA'} (no se ha calculado ninguna configuracion del gate ni escrito nada)")
    return 0 if ok else 1


def main():
    t0 = time.time()
    if any(os.path.exists(os.path.join(OUT, f)) for f in NAMES):
        print("Rehuso sobrescribir resultados E1-bis existentes."); return 3
    rows_int = []
    # ---- I0 ----
    i0, msgs, head, pre_commit = check_I0()
    rows_int.append(dict(prueba='I0_criterios_en_commit', m_u='', valor='; '.join(msgs) or f'HEAD={head}',
                         criterio='archivos en git y sin cambios', ok=i0))
    print(f"I0 criterios en commit: {'OK' if i0 else 'NO'}  HEAD={head}  ultimo commit del pre-registro: {pre_commit}")
    if not i0:
        for s in msgs: print("   " + s)
        print("STOP: haz commit de los criterios antes de ejecutar. No se ha calculado nada."); return 2

    # ---- I1 ----
    _, _, T57 = E.lrcp_theory(0.6); om1 = T57*WINS[W_MAIN][0]
    nb, _ = nL_point(0.6, DXS[0], om1, 'base'); ne, _ = nL_point(0.6, DXS[0], om1, SOLVER)
    I1 = abs(nb - ne)/abs(nb) if (nb is not None and ne is not None) else np.inf
    rows_int.append(dict(prueba='I1_eqref_vs_spsolve', m_u=0.6, valor=I1, criterio=f"<= {TOL['I1']}", ok=I1 <= TOL['I1']))
    print(f"I1 eqref frente a spsolve (m_u=0.6, dx={DXS[0]}, omega={om1:.4e}): {I1:.2e}")
    if not I1 <= TOL['I1']:
        print("STOP I1"); return finish(t0, rows_int, [], [], {}, head, pre_commit, 'STOP I1')

    rows_nL, rows_fit = [], []; A = {}
    for m in MU_CONF + MU_REPL:
        a_th, b_th, T57 = E.lrcp_theory(m)
        print(f"\n=== m_u={m} ({'confirmatoria' if m in MU_CONF else 'replica, no cuenta'})  a_A8={a_th:.10f}  b_A14={b_th:.10f} ===")
        for wname, (lo, hi) in WINS.items():
            oms = T57*np.geomspace(lo, hi, NW)
            for dx in DXS:
                nl = []
                for om in oms:
                    v_, ch = nL_point(m, dx, om, SOLVER)
                    if v_ is None:
                        rows_int.append(dict(prueba='I2_canales', m_u=m, valor=f"dx={dx} {wname} omega={om:.6e} {ch}",
                                             criterio='Lin=1 Rin=2 Rin_neg=1 Lout=1 Rout=2', ok=False))
                        print(f"STOP I2: {ch}"); return finish(t0, rows_int, rows_nL, rows_fit, {}, head, pre_commit, 'STOP I2')
                    nl.append(v_)
                    rows_nL.append(dict(m_u=m, ventana=wname, dx=dx, N=E.grid(L_GATE, dx), omega=om, omega_over_TH=om/T57, n_L=v_))
                a, b = E.fit_ab(oms, np.array(nl), ORDER); A[(m, wname, dx)] = (a, b)
                rows_fit.append(dict(m_u=m, ventana=wname, dx=dx, orden=ORDER, a=a, b=b,
                                     a_rel=(a-a_th)/a_th, b_rel=(b-b_th)/abs(b_th)))
                print(f"   {wname} dx={dx}: a err {(a-a_th)/a_th:+.3e}  b err {(b-b_th)/abs(b_th):+.3e}   [{time.time()-t0:.0f} s]", flush=True)
        rows_int.append(dict(prueba='I2_canales', m_u=m, valor='todas las omega', criterio='Lin=1 Rin=2 Rin_neg=1 Lout=1 Rout=2', ok=True))

    # ---- I3 y criterios ----
    def rich(m, w, idx):
        q = [A[(m, w, dx)][idx] for dx in DXS]
        return q[2] + (q[2] - q[1])/3.0, q
    summary = {}; stop = None
    for m in MU_CONF + MU_REPL:
        a_th, b_th, _ = E.lrcp_theory(m)
        aR, q = rich(m, W_MAIN, 0)
        den = abs(q[1] - q[2])
        p_obs = float(np.log2(abs(q[0] - q[1])/den)) if den > 0 else np.nan
        okp = bool(np.isfinite(p_obs) and TOL['p_lo'] <= p_obs <= TOL['p_hi'])
        rows_int.append(dict(prueba='I3_orden_observado_W2', m_u=m, valor=p_obs, criterio=f"in [{TOL['p_lo']},{TOL['p_hi']}]", ok=okp))
        if m in MU_CONF and not okp:
            stop = f'I3 m_u={m}'
        aR3, _ = rich(m, 'W3', 0); bR, _ = rich(m, W_MAIN, 1)
        summary[m] = dict(eps_a=abs(aR - a_th)/a_th, win=abs(aR - aR3)/a_th, eps_b=abs(bR - b_th)/abs(b_th),
                          aR=aR, aR3=aR3, bR=bR, p_obs=p_obs)
    if stop:
        print(f"STOP {stop}"); return finish(t0, rows_int, rows_nL, rows_fit, summary, head, pre_commit, f'STOP {stop}')

    crit = []
    for m in MU_CONF:
        s = summary[m]
        crit += [dict(criterio='C1', m_u=m, valor=s['eps_a'], tol=TOL['C1'], ok=s['eps_a'] <= TOL['C1']),
                 dict(criterio='C2', m_u=m, valor=s['win'], tol=TOL['C2'], ok=s['win'] <= TOL['C2']),
                 dict(criterio='C3', m_u=m, valor=s['eps_a'], tol=TOL['C3'], ok=s['eps_a'] <= TOL['C3']),
                 dict(criterio='C4', m_u=m, valor=s['eps_b'], tol=TOL['C4'], ok=s['eps_b'] <= TOL['C4'])]
    allok = lambda c: all(r['ok'] for r in crit if r['criterio'] == c)
    verdict = 'PASS' if (allok('C1') and allok('C2')) else 'FAIL'
    print("\n================ RESULTADO E1-bis ================")
    for m in MU_CONF + MU_REPL:
        s = summary[m]
        print(f"m_u={m} {'(confirmatoria)' if m in MU_CONF else '(replica, NO cuenta)':22s} eps_a={s['eps_a']:.2e}  "
              f"|W2-W3|={s['win']:.2e}  eps_b={s['eps_b']:.2e}  p_obs={s['p_obs']:.2f}")
    print(f"C1 (eps_a <= {TOL['C1']:.0e}): {allok('C1')}   C2 (ventana <= {TOL['C2']:.0e}): {allok('C2')}")
    print(f"C3 (nivel fuerte <= {TOL['C3']:.0e}, solo se informa): {allok('C3')}   C4 (b <= {TOL['C4']:.0e}, secundario): {allok('C4')}")
    print(f"VEREDICTO E1-bis: {verdict}")
    return finish(t0, rows_int, rows_nL, rows_fit, summary, head, pre_commit, verdict, crit)


def finish(t0, rows_int, rows_nL, rows_fit, summary, head, pre_commit, verdict, crit=()):
    os.makedirs(OUT, exist_ok=True)
    def wcsv(name, rows):
        if not rows: return
        with open(os.path.join(OUT, name), 'w', newline='') as f:
            w = csv.DictWriter(f, fieldnames=list(rows[0].keys())); w.writeheader(); w.writerows(rows)
    wcsv('e1bis_nL.csv', rows_nL); wcsv('e1bis_fits.csv', rows_fit)
    wcsv('e1bis_result.csv', [dict(tipo='integridad', **r) for r in rows_int] +
         [dict(tipo='criterio', prueba=r['criterio'], m_u=r['m_u'], valor=r['valor'], criterio=f"<= {r['tol']}", ok=r['ok']) for r in crit])
    files = {os.path.relpath(p, ROOT): sha(p) for p in TRACKED + [os.path.join(HERE, 'gate_E1_lrcp.py'),
                                                                  os.path.join(ROOT, 'solver', 'fd_bdg_operator.py')]}
    for n_ in NAMES[:3]:
        p = os.path.join(OUT, n_)
        if os.path.exists(p): files['out/' + n_] = sha(p)
    man = dict(gate='E1-bis', veredicto=verdict, git_HEAD=head, commit_pre_registro=pre_commit, sha256=files,
               resumen={str(k): v for k, v in summary.items()}, tolerancias=TOL, python=platform.python_version(),
               numpy=np.__version__, scipy=scipy.__version__, plataforma=platform.platform(),
               segundos=round(time.time() - t0, 1), fecha=time.strftime('%Y-%m-%dT%H:%M:%S'))
    with open(os.path.join(OUT, 'e1bis_manifest.json'), 'w') as f:
        json.dump(man, f, indent=2, default=float)
    print(f"fin ({time.time()-t0:.0f} s)")
    return 0 if verdict in ('PASS', 'FAIL') else 2


if __name__ == '__main__':
    try:
        sys.exit(selftest() if '--selftest' in sys.argv else main())
    except Exception:
        traceback.print_exc(); sys.exit(1)
