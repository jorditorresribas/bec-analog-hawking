"""
BLOQUE 2 — B2.4c: VERIFICACION SVD de la contradiccion estimador/adjunto a X0 grande.

Contradiccion (B2.4b): ||y||=4.79e8 en M^dag y = c implica sigma_min <= 2e-9, PERO el estimador dio
sigma_min=6.9e-5. ChatGPT: ||y|| solo acota sigma_min si el adjunto esta BIEN RESUELTO (r_adj pequeno).

Para X0=20,24,28 (dx=0.02), sobre la MISMA Mrep (fisica) usada en el adjunto:
  - sigma_min, sigma_max, cond por SVD densa (verdad de terreno).
  - residuo PRIMAL  r_pri = ||Mrep z - brep|| / ||brep||.
  - residuo ADJUNTO r_adj = ||Mrep^dag y - crep|| / ||crep||,  y  ||y||, ||crep||.
  - comparacion con el estimador (sin modificarlo).
  - pico RSS del proceso.

Tabla de decision (ChatGPT):
  SVD ~2e-9                      -> el estimador fallo.
  SVD >> 2e-9  y r_adj grande    -> el adjunto no resolvio el sistema.
  SVD >> 2e-9  y r_adj pequeno   -> inconsistencia de representacion/escalado/c.

NO modifica B2.1/2/3/v1.0. NO toca el log.
"""
import resource
import numpy as np
import scipy.linalg as la
import gate_B2_1_algebraic as g1
import gate_B2_2_conditioning as g2
import gate_B2_4_X0_convergence as b24


def verify(w, X0, dx=0.02, out_sel=0, in_prob=0):
    N = int(round(2*X0/dx)); ndof = 4*(N+1)
    print("="*96)
    print(f"[B2.4c] w={w}  X0={X0}  N={N}  ndof={ndof}")
    print("="*96)
    M, info = g1.assemble(w, N, X0)
    dr, dc = g2.physical_scaling(info)
    Mrep = dr[:, None] * M * dc[None, :]
    del M
    # --- estimador (igual que B2.4, sin strict para que reporte aunque "falle") ---
    est = b24.estim_metrics(Mrep, info["ndof"], strict=False)
    # --- primal + adjunto con una sola LU ---
    lu = la.lu_factor(Mrep)
    sin = np.zeros(3); sin[in_prob] = 1.0
    brep = dr * g1.rhs_for(info, sin)
    z = la.lu_solve(lu, brep)
    r_pri = la.norm(Mrep @ z - brep) / la.norm(brep)
    outs = g2.out_indices(info)
    c = np.zeros(info["ndof"], dtype=complex); c[outs[out_sel]] = 1.0
    crep = np.conj(dc) * c
    y = la.lu_solve(lu, crep, trans=2)
    r_adj = la.norm(Mrep.conj().T @ y - crep) / la.norm(crep)
    ny = la.norm(y)
    del lu
    # --- SVD densa (verdad de terreno); overwrite para minimizar memoria ---
    sv = la.svdvals(Mrep, overwrite_a=True)
    smin_svd, smax_svd = sv[-1], sv[0]; cond_svd = smax_svd/smin_svd
    bound = 1.0/ny                                   # sigma_min <= ||c||/||y|| (||c||=1)
    peak = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss / 1e9
    print(f"  SVD (verdad):   sigma_min={smin_svd:.4e}  sigma_max={smax_svd:.4e}  cond={cond_svd:.4e}")
    print(f"  Estimador:      sigma_min={est['sigma_min']:.4e}  cond={est['cond']:.4e}  "
          f"zona={est['zone']}  convergido={est['converged']}")
    print(f"  Cota por ||y||: sigma_min <= ||c||/||y|| = 1/{ny:.3e} = {bound:.4e}")
    print(f"  Residuos:       r_primal={r_pri:.2e}   r_adjunto={r_adj:.2e}   ||y||={ny:.3e}  ||c||={la.norm(crep):.3e}")
    # veredicto automatico segun la tabla
    if smin_svd <= 3e-9:
        verdict = "SVD confirma sigma_min minusculo -> EL ESTIMADOR FALLO (cond real explota)."
    elif r_adj > 1e-6:
        verdict = "SVD grande y r_adj GRANDE -> el adjunto NO resolvio; ||y|| NO acota sigma_min."
    else:
        verdict = "SVD grande y r_adj pequeno -> inconsistencia de representacion/escalado/c (investigar)."
    print(f"  >>> {verdict}")
    print(f"  pico RSS ~ {peak:.2f} GB")
    del Mrep, sv
    return dict(X0=X0, smin_svd=smin_svd, cond_svd=cond_svd, smin_est=est['sigma_min'],
                r_pri=r_pri, r_adj=r_adj, ny=ny, verdict=verdict)


if __name__ == "__main__":
    res = []
    for X0 in (20, 24, 28):                 # ascendente: si 28 fallara por memoria, ya tenemos 20,24
        res.append(verify(0.10, X0, dx=0.02))
    print("\n" + "="*96)
    print("RESUMEN  w=0.10")
    print(f"  {'X0':>5} {'sigma_min_SVD':>14} {'sigma_min_est':>14} {'cond_SVD':>12} {'r_adj':>10} {'||y||':>11}")
    for r in res:
        print(f"  {r['X0']:5.1f} {r['smin_svd']:14.4e} {r['smin_est']:14.4e} {r['cond_svd']:12.4e} "
              f"{r['r_adj']:10.2e} {r['ny']:11.3e}")
    print(">>> B2.4c COMPLETO")
