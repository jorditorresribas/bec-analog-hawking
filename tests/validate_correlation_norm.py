"""
TEST de normalizacion ABSOLUTA de la correlacion densidad-densidad (correlations.py), anadido en la auditoria del
2026-11-18 (origen: residual_anatomy/diag_recati_norm.py y diag_recati_direct.py).

En un condensado homogeneo subsonico (Mach 0,5) la correlacion de vacio es exacta:
    <dn(x) dn(x')>_conexa sin termino propio = n * int dk/2pi [S(k) - 1] e^{ik(x-x')},   S = eps_k / E_k.
Se construye C con el MISMO codigo que density_correlation (scatter_field, mode_bracket, CORR_PREF; canales u y d1,
los unicos que existen) y se compara con la formula exacta sobre el mismo rango de k que cubren las omegas.
Criterios (fijados antes de ejecutar esta version del test):
  - factor de ajuste codigo/exacto en [1 - 1e-3, 1 + 1e-3]   (con el prefactor antiguo 1/(4 pi) daba 0,500001)
  - residuo de forma tras el factor <= 1e-3 (relativo al maximo)
  - control negativo: con 1/(4 pi) el factor tiene que salir fuera de tolerancia
Sale con codigo 1 si falla algo. ~1 min.
"""
import os, sys, io, contextlib
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
sys.path.insert(0, os.path.join(ROOT, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C

TOL_FACTOR, TOL_SHAPE = 1e-3, 1e-3
L, N, NC = 40.0, 4000, 241
idx = np.linspace(0, N-1, NC).astype(int); x = np.linspace(-L, L, N)[idx]
OM = np.linspace(0.004, 1.2, 150); wts = np.full(len(OM), OM[1]-OM[0]); wts[[0, -1]] *= 0.5
acc = np.zeros((NC, NC), complex); kp, km = [], []
for om, wt in zip(OM, wts):
    modes = {}
    for ch in ('u', 'd1'):
        r = C.scatter_field(L, N, om, ch, 2.0, 1.0, 'homog')
        if r['Rin_neg'] is not None:
            print("FALLA: aparece un canal de norma negativa en el homogeneo subsonico"); sys.exit(1)
        (kp if ch == 'u' else km).append((r['Lin'] if ch == 'u' else r['Rin_pos'])['k'])
        u = r['u'][idx]; vv = r['v_mode'][idx]; modes[ch] = dict(u=u, v=vv, r=u+vv)
    acc += C.mode_bracket(modes, 'all')*wt
n0, g = C.n_u, C.g
Ccode = np.real(n0*C.CORR_PREF*acc)


def exact(d):
    out = 0.0
    for a, b in ((min(kp), max(kp)), (min(km), max(km))):
        kk = np.linspace(a, b, 20001); e = kk**2/2; E = np.sqrt(e*(e + 2*g*n0))
        f = (e/E - 1)*np.cos(kk*d)
        out += np.sum(0.5*(f[1:] + f[:-1])*np.diff(kk))/(2*np.pi)
    return n0*out


ic = NC//2; seps = np.arange(0, 40)
cc = np.array([Ccode[ic, ic+s] for s in seps]); ee = np.array([exact(x[ic+s] - x[ic]) for s in seps])
fac = float(np.sum(cc*ee)/np.sum(ee*ee)); shape = float(np.max(np.abs(cc - fac*ee))/np.max(np.abs(ee)))
fac_old = fac*(1/(4*np.pi))/C.CORR_PREF
ok1 = abs(fac - 1) <= TOL_FACTOR; ok2 = shape <= TOL_SHAPE; ok3 = abs(fac_old - 1) > TOL_FACTOR
print(f"CORR_PREF = {C.CORR_PREF:.6f} (1/(2 pi) = {1/(2*np.pi):.6f})")
print(f"[{'OK' if ok1 else 'FALLA'}] factor codigo/exacto = {fac:.6f}  (tol |f-1| <= {TOL_FACTOR})")
print(f"[{'OK' if ok2 else 'FALLA'}] residuo de forma = {shape:.1e}  (tol <= {TOL_SHAPE})")
print(f"[{'OK' if ok3 else 'FALLA'}] control negativo: con 1/(4 pi) el factor seria {fac_old:.6f} (debe quedar fuera)")
sys.exit(0 if (ok1 and ok2 and ok3) else 1)
