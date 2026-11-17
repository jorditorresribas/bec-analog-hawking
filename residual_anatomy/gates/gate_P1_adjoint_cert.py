"""
P1 — CERTIFICACION de la ruta adjunta (plan CORRECTION_PLAN.md).

Certifica M_rep^dag y = c_rep con >=2 rutas DIRECTAS independientes (R1 LU, R2 QR) que coincidan;
R3 (GMRES con precond. ADJUNTO) como tercera comprobacion; diagnostico del fallo de la.solve(M^H)
SIN atribuir causa no demostrada. NO recalcula kappa_out. Aditivo: importa B2.1/B2.2, no los modifica.

Controles (ChatGPT): integridad de entradas (hash LOGICO de valores vs FISICO de memoria, disposicion,
mutacion); residuos abs/rel del sistema ORIGINAL; identidad s=c^dag x vs y^dag b; cota SVD para solucion
APROXIMADA ||y||<=(||c||+||r||)/sigma_min y error ||y-y*||<=||r||/sigma_min; certificacion solo donde
la SVD es verdad-terreno, evidencia de consistencia (no certificacion) en tamanos mayores; manifiesto.
"""
import hashlib, platform, time, sys
import numpy as np
import scipy.linalg as la
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import gate_B2_1_algebraic as g1
import gate_B2_2_conditioning as g2
_OUTDIR = str(__import__("pathlib").Path(__file__).resolve().parents[2] / "out")  # out/ en la raíz del repo

EPS = np.finfo(np.complex128).eps


# ----------------------------- integridad -------------------------------------
def integrity(arr):
    """hash LOGICO (valores, orden C canonico) vs FISICO (buffer tal como se almacena), disposicion.
    SHA-256 COMPLETO (trazabilidad documental)."""
    logical = hashlib.sha256(np.ascontiguousarray(arr).tobytes(order='C')).hexdigest()
    physical = hashlib.sha256(arr.tobytes(order='A')).hexdigest()
    return dict(logical=logical, physical=physical, shape=arr.shape, dtype=str(arr.dtype),
                C=arr.flags['C_CONTIGUOUS'], F=arr.flags['F_CONTIGUOUS'])


def same_values(a, b):
    """igualdad MATEMATICA independiente de disposicion C/F."""
    return a.shape == b.shape and a.dtype == b.dtype and np.array_equal(a, b)


# ----------------------------- sistema ----------------------------------------
def build_system(w, X0, dx, in_prob=0, out_sel=0):
    N = int(round(2*X0/dx))
    M, info = g1.assemble(w, N, X0)
    dr, dc = g2.physical_scaling(info)
    Mrep = dr[:, None] * M * dc[None, :]
    sin = np.zeros(3); sin[in_prob] = 1.0
    brep = (dr * g1.rhs_for(info, sin)).astype(np.complex128)
    outs = g2.out_indices(info)
    c = np.zeros(info["ndof"], dtype=np.complex128); c[outs[out_sel]] = 1.0
    crep = np.conj(dc) * c                                   # c_rep = D_c^dag e_out
    return np.ascontiguousarray(Mrep), brep, crep, info, N


def resid_adj(Mrep, y, crep):
    r = Mrep.conj().T @ y - crep
    ra = la.norm(r)
    return ra, ra/la.norm(crep)


# ----------------------------- rutas ------------------------------------------
def R1_lu(Mrep, crep):
    lu = la.lu_factor(Mrep.copy())                          # copia independiente; sin overwrite
    y = la.lu_solve(lu, crep.copy(), trans=2)               # resuelve M^H y = c
    return y, lu


def R2_qr(Mrep, crep):
    # Q,R = qr(Mrep); M^H = R^H Q^H ; R^H w = c ; y = Q w. (Householder via scipy; Q economico.)
    Q, R = la.qr(Mrep.copy(), mode='economic')
    w = la.solve_triangular(R.conj().T, crep.copy(), lower=True)
    y = Q @ w
    return y


def R3_gmres(Mrep, crep, fill_factor=10, drop_tol=1e-4, restart=200, maxiter=2000, rtol=1e-12):
    """GMRES sobre M^H y = c con precond. ADJUNTO de ILU(M_rep) (trans='H' respeta permutaciones)."""
    n = Mrep.shape[0]
    Msp = sp.csc_matrix(Mrep)
    ilu = spla.spilu(Msp, fill_factor=fill_factor, drop_tol=drop_tol)
    # A = M^H sin formar M^H:  M^H x = conj(M^T conj(x))
    def mvA(x):  return np.conj(Mrep.T @ np.conj(x))
    A = spla.LinearOperator((n, n), matvec=mvA, dtype=np.complex128)
    # Precond ~ (M^H)^{-1}: resolucion ADJUNTA de los factores ILU (respeta permutaciones)
    def mvP(x):  return ilu.solve(x, trans='H')
    P = spla.LinearOperator((n, n), matvec=mvP, dtype=np.complex128)
    res_pre = []
    y, infoc = spla.gmres(A, crep.copy(), M=P, rtol=rtol, atol=0.0, restart=restart,
                          maxiter=maxiter, callback=lambda rk: res_pre.append(float(rk)),
                          callback_type='pr_norm')
    return y, dict(info=int(infoc), iters=len(res_pre), nnz_ilu=ilu.L.nnz + ilu.U.nnz)


# -------------------- diagnostico del fallo (sin atribuir causa) --------------
def diagnose_lasolve(Mrep, crep, lu_R1):
    """D-old vs D-MH + factorizacion/resolucion/memoria. Resultado: identificado / acotado / NO identif."""
    out = {}
    # D-old: la.solve sobre M^H (vista conj.T) — ruta defectuosa historica
    y_old = la.solve(Mrep.conj().T, crep.copy())
    out['D_old'] = dict(ny=float(la.norm(y_old)), r_rel=float(resid_adj(Mrep, y_old, crep)[1]))
    # D-MH: formar M^H contiguo y factorizar/solucionar explicitamente
    MH = np.ascontiguousarray(Mrep.conj().T)
    lu_mh, piv_mh = la.lu_factor(MH.copy())
    y_mh = la.lu_solve((lu_mh, piv_mh), crep.copy())
    # SEPARAR factorizacion vs resolucion del M^H FORMADO (convencion de pivoteo SciPy: swaps secuenciales)
    Lh = np.tril(lu_mh, -1) + np.eye(lu_mh.shape[0]); Uh = np.triu(lu_mh)
    Ph = np.eye(lu_mh.shape[0])
    for i, p in enumerate(piv_mh):
        Ph[[i, p]] = Ph[[p, i]]
    r_fact_MH = float(la.norm(Ph @ MH - Lh @ Uh)/la.norm(MH))     # ||P M^H - L U||_F / ||M^H||_F
    r_solve_MH = float(resid_adj(Mrep, y_mh, crep)[1])            # ||M^H y - c|| / ||c||
    out['D_MH'] = dict(ny=float(la.norm(y_mh)), r_rel=r_solve_MH,
                       r_fact=r_fact_MH, r_solve=r_solve_MH)
    # control de factorizacion: residuo de reconstruccion P A = L U de la LU de M (R1)
    lu, piv = lu_R1
    Lm = np.tril(lu, -1) + np.eye(lu.shape[0]); Um = np.triu(lu)
    P = np.eye(lu.shape[0])
    for i, p in enumerate(piv):
        P[[i, p]] = P[[p, i]]
    out['fact_PA_minus_LU'] = float(la.norm(P @ Mrep - Lm @ Um)/la.norm(Mrep))
    # memoria: repetir D-MH con disposicion F
    MH_F = np.asfortranarray(Mrep.conj().T)
    y_mhF = la.lu_solve(la.lu_factor(MH_F.copy()), crep.copy())
    out['D_MH_Forder'] = dict(ny=float(la.norm(y_mhF)), r_rel=float(resid_adj(Mrep, y_mhF, crep)[1]))
    # veredicto (prudente): solo clasifica lo que los controles distinguen
    old_bad = out['D_old']['r_rel'] > 1e-6
    mh_bad  = out['D_MH']['r_rel'] > 1e-6
    if old_bad and mh_bad:
        out['verdict'] = "fallo en AMBAS (la.solve y lu de M^H formado) -> NO atribuible solo a la.solve; etapa no distinguida -> ACOTADO a resolver/factorizar M^H formado"
    elif old_bad and not mh_bad:
        out['verdict'] = "la.solve(M^H) falla pero lu_factor+lu_solve de M^H acierta -> diferencia de ruta; causa exacta NO identificada (copia/params)"
    elif not old_bad:
        out['verdict'] = "D-old NO reproduce el fallo en este tamano/proceso -> dependiente de escala/entorno; NO identificado aqui"
    return out


# ----------------------------- verdad-terreno ---------------------------------
def svd_bounds(Mrep, crep, y, feasible=True):
    if not feasible:
        return None
    sv = la.svdvals(Mrep)
    smin = sv[-1]
    r = resid_adj(Mrep, y, crep)[0]
    bound_y = (la.norm(crep) + r)/smin                      # ||y|| <= (||c||+||r||)/smin  (sol aprox)
    err_bound = r/smin                                      # ||y - y*|| <= ||r||/smin
    return dict(sigma_min=float(smin), ny=float(la.norm(y)), bound_y=float(bound_y),
                within=bool(la.norm(y) <= bound_y*(1+1e-9)), err_bound=float(err_bound))


# ----------------------------- orquestacion -----------------------------------
def run(w=0.10, X0=12.0, dx=0.02, in_prob=0, out_sel=0, do_qr=True, do_svd=True):
    print("="*98)
    print(f"[P1] w={w} X0={X0} dx={dx} in_prob={in_prob} out_sel={out_sel}")
    print("="*98)
    Mrep, brep, crep, info, N = build_system(w, X0, dx, in_prob, out_sel)
    ndof = info["ndof"]
    print(f"  N={N} ndof={ndof}")
    # --- integridad de entradas (copias independientes por ruta) ---
    iM, ib, ic = integrity(Mrep), integrity(brep), integrity(crep)
    print(f"  integridad Mrep: log={iM['logical']} phys={iM['physical']} {iM['shape']} C={iM['C']} F={iM['F']}")
    before = (iM['logical'], ib['logical'], ic['logical'])
    # --- R1 (LU) ---
    y1, lu1 = R1_lu(Mrep, crep)
    ra1_abs, ra1_rel = resid_adj(Mrep, y1, crep)
    x1 = la.lu_solve(lu1, brep.copy())                      # primal para identidad
    s_dir = np.vdot(crep, x1); s_adj = np.vdot(y1, brep)
    print(f"  [R1 LU ] ||y||={la.norm(y1):.4e} r_adj_rel={ra1_rel:.2e} r_adj_abs={ra1_abs:.2e}  "
          f"s=c^dag x={s_dir:+.4e} y^dag b={s_adj:+.4e} |dif|={abs(s_dir-s_adj):.1e}")
    # --- R2 (QR) solo donde viable ---
    y2 = None
    if do_qr:
        y2 = R2_qr(Mrep, crep)
        ra2_abs, ra2_rel = resid_adj(Mrep, y2, crep)
        print(f"  [R2 QR ] ||y||={la.norm(y2):.4e} r_adj_rel={ra2_rel:.2e} r_adj_abs={ra2_abs:.2e}  "
              f"||y2-y1||={la.norm(y2-y1):.2e}")
    else:
        print(f"  [R2 QR ] OMITIDA (N grande, QR densa no viable) -> solo evidencia por R1/R3")
    # --- R3 (GMRES precond adjunto) ---
    y3, g3 = R3_gmres(Mrep, crep)
    ra3_abs, ra3_rel = resid_adj(Mrep, y3, crep)
    print(f"  [R3 GMR] ||y||={la.norm(y3):.4e} r_adj_rel={ra3_rel:.2e} (ORIGINAL) info={g3['info']} "
          f"iters={g3['iters']} ||y3-y1||={la.norm(y3-y1):.2e}")
    # --- cota SVD (verdad-terreno) donde viable ---
    sb = svd_bounds(Mrep, crep, y1, feasible=do_svd)
    if sb:
        print(f"  [SVD  ] sigma_min={sb['sigma_min']:.3e}  ||y1||={sb['ny']:.3e} <= bound={sb['bound_y']:.3e}? "
              f"{sb['within']}  err_bound(||y-y*||)={sb['err_bound']:.2e}")
    else:
        print(f"  [SVD  ] OMITIDA -> fuera de verdad-terreno: EVIDENCIA de consistencia, NO certificacion")
    # --- diagnostico del fallo ---
    diag = diagnose_lasolve(Mrep, crep, lu1)
    print(f"  [DIAG ] D-old r_rel={diag['D_old']['r_rel']:.2e} | D-MH r_rel={diag['D_MH']['r_rel']:.2e} "
          f"(r_fact={diag['D_MH']['r_fact']:.1e} r_solve={diag['D_MH']['r_solve']:.2e}) | "
          f"D-MH(F) r_rel={diag['D_MH_Forder']['r_rel']:.2e} | PA-LU(R1)={diag['fact_PA_minus_LU']:.1e}")
    print(f"  [DIAG ] veredicto: {diag['verdict']}")
    # --- mutacion de entradas tras resolver ---
    after = (integrity(Mrep)['logical'], integrity(brep)['logical'], integrity(crep)['logical'])
    print(f"  [MUT  ] entradas sin mutar tras resolver: {before == after}")
    ra2a, ra2r = (np.nan, np.nan) if y2 is None else resid_adj(Mrep, y2, crep)
    rec = dict(w=w, X0=X0, N=N, ndof=ndof, in_prob=in_prob, out_sel=out_sel,
               R1_radj_abs=ra1_abs, R1_radj_rel=ra1_rel, R1_ny=float(la.norm(y1)),
               R2_radj_rel=(np.nan if y2 is None else ra2r), R2_diff_y1=(np.nan if y2 is None else float(la.norm(y2-y1))),
               R3_radj_rel=ra3_rel, R3_info=g3['info'], R3_iters=g3['iters'], R3_diff_y1=float(la.norm(y3-y1)),
               r_primal=float(la.norm(Mrep@x1-brep)/la.norm(brep)),
               s_dir=complex(s_dir), s_adj=complex(s_adj), s_dif=float(abs(s_dir-s_adj)),
               sigma_min=(np.nan if sb is None else sb['sigma_min']),
               bound_y=(np.nan if sb is None else sb['bound_y']), within_bound=(None if sb is None else sb['within']),
               err_bound=(np.nan if sb is None else sb['err_bound']),
               Dold_rrel=diag['D_old']['r_rel'], DMH_rrel=diag['D_MH']['r_rel'],
               DMH_rfact=diag['D_MH']['r_fact'], DMH_rsolve=diag['D_MH']['r_solve'],
               DMH_F_rrel=diag['D_MH_Forder']['r_rel'], R1_PA_LU=diag['fact_PA_minus_LU'],
               diag_verdict=diag['verdict'], inputs_mutated=(before != after),
               hash_M=iM['logical'], hash_b=ib['logical'], hash_c=ic['logical'],
               Mrep_C=iM['C'], Mrep_F=iM['F'])
    return rec


def manifest():
    print("\n--- MANIFIESTO ---")
    print(f"  entorno: Python {platform.python_version()} numpy {np.__version__} scipy "
          f"{__import__('scipy').__version__} {platform.platform()}")
    for f in ('gate_P1_adjoint_cert.py', 'gate_B2_1_algebraic.py', 'gate_B2_2_conditioning.py'):
        try:
            h = hashlib.sha256(open(f, 'rb').read()).hexdigest()[:16]
            print(f"  hash {f}: {h}")
        except Exception as e:
            print(f"  hash {f}: (no leido: {e})")


def _sha(path):
    try: return hashlib.sha256(open(path, 'rb').read()).hexdigest()
    except Exception as e: return f"(no leido {e})"


def write_csv_and_manifest(recs, path_csv, path_man, cmd, exit_code=0, expected_rows=None):
    import csv, os
    os.makedirs(os.path.dirname(path_csv), exist_ok=True)
    keys = list(recs[0].keys())
    with open(path_csv, 'w', newline='') as f:
        wtr = csv.DictWriter(f, fieldnames=keys); wtr.writeheader()
        for r in recs: wtr.writerow(r)
    csv_sha = _sha(path_csv)
    try:
        blas = str(np.__config__.CONFIG.get('Build Dependencies', {}).get('blas', {}))
    except Exception:
        blas = "(no disponible)"
    npz = _OUTDIR + "/waterfall_state.npz"
    with open(path_man, 'w') as f:
        f.write(f"P1 MANIFIESTO — {time.strftime('%Y-%m-%d %H:%M:%S %Z')}\n")
        f.write(f"comando: {cmd}\n")
        f.write(f"exit_code: {exit_code}\n")
        f.write(f"python {platform.python_version()}  numpy {np.__version__}  scipy "
                f"{__import__('scipy').__version__}\n")
        f.write(f"plataforma: {platform.platform()}\n")
        f.write(f"BLAS/LAPACK: {blas}\n")
        f.write(f"filas CSV: {len(recs)}" + (f" (esperado {expected_rows})" if expected_rows else "") + "\n")
        f.write(f"sha256 {os.path.basename(path_csv)} (CSV): {csv_sha}\n")
        f.write(f"sha256 waterfall_state.npz (datos de entrada): {_sha(npz)}\n")
        for fn in ('gate_P1_adjoint_cert.py', 'gate_B2_1_algebraic.py', 'gate_B2_2_conditioning.py'):
            f.write(f"sha256 {fn}: {_sha(fn)}\n")
        f.write("nota: hashes de M,b,c en el CSV son SHA-256 COMPLETOS; salida stdout conservada aparte.\n")
    print(f"\n  CSV -> {path_csv}  ({len(recs)} filas)  sha256={csv_sha[:16]}...\n  MANIFIESTO -> {path_man}")


def run_selectors(w=0.10, X0=12.0, dx=0.02):
    """Pase POR SELECTOR a un tamano representativo (A de ChatGPT): out_sel 0,1,2 con R1/R2/R3/SVD."""
    print("="*98); print(f"### POR SELECTOR (representativo) w={w} X0={X0}  out_sel=0,1,2 ###"); print("="*98)
    recs = [run(w, X0, dx, in_prob=0, out_sel=s, do_qr=True, do_svd=True) for s in (0, 1, 2)]
    write_csv_and_manifest(recs, _OUTDIR + "/p1_selectors_X012.csv",
                           _OUTDIR + "/p1_selectors_manifest.txt",
                           "python3 -c 'import gate_P1_adjoint_cert as p; p.run_selectors()'",
                           exit_code=0, expected_rows=3)
    return recs


if __name__ == "__main__":
    manifest()
    recs = []
    print("\n### CERTIFICACION (X0 con SVD + QR viables): X0=12, 20 ###")
    recs.append(run(0.10, X0=12.0, do_qr=True,  do_svd=True))
    recs.append(run(0.10, X0=20.0, do_qr=True,  do_svd=True))
    print("\n### EVIDENCIA DE CONSISTENCIA (N mayor, QR densa costosa): X0=28 ###")
    recs.append(run(0.10, X0=28.0, do_qr=False, do_svd=True))
    write_csv_and_manifest(
        recs, _OUTDIR + "/p1_adjoint_cert.csv",
        _OUTDIR + "/p1_adjoint_manifest.txt",
        "python3 residual_anatomy/gates/gate_P1_adjoint_cert.py")
    print("\n>>> P1 COMPLETO")
