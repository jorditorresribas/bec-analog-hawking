"""
DIAGNOSTICO (auditoria fase 1, 2026-11-17) — equivalencia de escala del knob lam con fondo consistente.

Afirmacion: el operador BdG con dispersion (omega - v k)^2 = c^2 k^2 + lam k^4/4 y fondo CONSISTENTE
(Qp_lam = sqrt(lam)*Qp, variante QP_LAM de diag_T2_qp_consistent.py) sobre el perfil tanh de anchura w
es EXACTAMENTE el operador fisico v1.0 (lam = 1) sobre el perfil de anchura w/sqrt(lam), reescalado:
    x' = x/sqrt(lam),  k' = sqrt(lam) k,  omega' = sqrt(lam) omega,  M_1(x') = sqrt(lam) M_lam(x).
Comprobacion por terminos: cinetico sqrt(lam)/2 k^2 -> k'^2/2; acoplo g n/sqrt(lam) -> g n; conveccion
v k -> v k'; Qp_lam = sqrt(lam) Qp -> lam Qp = Qp' (Qp escala como 1/w^2). v, c y n no cambian, y
T/T_H es adimensional. En la malla, con L' = L/sqrt(lam) y el MISMO N: D2' = lam D2, D1' = sqrt(lam) D1,
asi que la identidad tambien es exacta en el operador DISCRETO (salvo redondeo).
Suelo de redondeo: el sistema lineal esta mal condicionado; perturbar L en 1 ulp (2e-16 relativo) cambia
n_L hasta ~7e-6 relativo con N=14400 (medido en la auditoria y re-medido aqui, ver NOISE). Por eso la
tolerancia es 1e-5 y no precision de maquina: diferencias por debajo del suelo de redondeo son ruido.
DECLARACION: la primera ejecucion usaba TOL = 1e-8 (fijada antes de medir el suelo de redondeo) y FALLO
(max 3.2e-6). La tolerancia se cambio DESPUES, tras medir el suelo (~7e-6). Es un cambio post-hoc
justificado por una medida independiente, pero post-hoc: el revisor debe valorarlo como tal.

Prueba: para cada (lam, N, omega) de out/diag_T2_nL.csv (variante QP_LAM, w=2, L=40) se recalcula n_L con
el solver v1.0 (correlations.scatter_field, op='validated', lam=1 -> fd_bdg_operator.build_MBdG, SIN usar
build_op_lambda) en w'=2/sqrt(lam), L'=40/sqrt(lam), mismo N, omega'=sqrt(lam)*omega, y se compara.
Consecuencia: con fondo consistente, el knob lam NO aisla fisica ultravioleta con la metrica fija: es el
mismo parametro que la anchura del horizonte en unidades de la longitud de curacion.

Salida: out/diag_scaling_equivalence.csv. Codigo de salida 1 si max diferencia relativa > TOL.
"""
import os, sys, io, csv, time, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C

TOL = 1e-5          # identidad exacta; tolerancia = suelo de redondeo del solver (ver NOISE abajo)
W, L = 2.0, 40.0    # condiciones de diag_T2_qp_consistent.py (box(2))

t0 = time.time()
ref = [r for r in csv.DictReader(open(os.path.join(OUT, 'diag_T2_nL.csv'))) if r['variante'] == 'QP_LAM']
rows, worst = [], 0.0
for r in ref:
    lam, N, om, nL_T2 = float(r['lam']), int(r['N']), float(r['omega']), float(r['n_L'])
    s = np.sqrt(lam)
    with contextlib.redirect_stdout(io.StringIO()):
        out = C.scatter_field(L/s, N, s*om, 'd2', w=W/s, lam=1.0, op='validated')
    nL_v1 = abs(out['out_amp'].get(('L', 1), 0.0))**2
    rel = abs(nL_v1 - nL_T2)/abs(nL_T2)
    worst = max(worst, rel)
    rows.append(dict(lam=lam, N=N, omega_T2=om, w_equiv=W/s, L_equiv=L/s, omega_equiv=s*om,
                     n_L_T2_QP_LAM=nL_T2, n_L_v10_rescaled=nL_v1, rel_diff=rel))
by = {}
for q in rows:
    by.setdefault((q['lam'], q['N']), []).append(q['rel_diff'])
for (lam, N), v in sorted(by.items()):
    print(f"  lam={lam:4.2f}  N={N:5d}  w'={W/np.sqrt(lam):.4f}  max rel diff n_L (9 omegas) = {max(v):.2e}")
with open(os.path.join(OUT, 'diag_scaling_equivalence.csv'), 'w', newline='') as f:
    wr = csv.DictWriter(f, fieldnames=list(rows[0].keys())); wr.writeheader(); wr.writerows(rows)
# NOISE: sensibilidad del propio solver v1.0 a 1 ulp en L, en el caso mas exigente (lam=4, N=14400, 1a omega)
q = next(q for q in rows if q['lam'] == 4.0 and q['N'] == 14400)
with contextlib.redirect_stdout(io.StringIO()):
    o2 = C.scatter_field(q['L_equiv']*(1+2e-16), 14400, q['omega_equiv'], 'd2', w=q['w_equiv'], lam=1.0, op='validated')
noise = abs(abs(o2['out_amp'].get(('L', 1), 0.0))**2 - q['n_L_v10_rescaled'])/q['n_L_v10_rescaled']
print(f"NOISE suelo de redondeo (v1.0 consigo mismo, L*(1+2e-16), lam=4 N=14400) = {noise:.2e}")
ok = worst < TOL
print(f"max diferencia relativa global = {worst:.2e}  (tol {TOL:.0e})  -> {'PASA' if ok else 'FALLA'}"
      f"   [{len(rows)} puntos, {time.time()-t0:.0f} s]  CSV: out/diag_scaling_equivalence.csv")
sys.exit(0 if ok else 1)
