"""
P2 — RE-REFERENCIA DE FASE DE S (implementacion del diseno CONGELADO residual_anatomy/P2_DESIGN.md).

Pone a prueba H0 (HIPOTESIS): a dx fijo, la variacion de S_ext con X0 es SOLO un re-faseo por canal con
k_box=(2/dx)arctan(k dx/2). El gate esta disenado para que H0 pueda FALLAR.

Alcance: solver LOCAL lambda=1 (A_at) de B2.1/B2.2, importados en SOLO LECTURA. NO valida la S fisica
(C1/P5), ni el Pade, ni n_L, ni kappa_out. NO calcula kappa_out ni adjuntos de salida: lu_solve(trans=2)
solo aparece DENTRO del estimador de sigma_min de B2.4 (ruta certificada en P1).

Estructura (P2_DESIGN.md):
  §1.1 base modal congelada (bit a bit) + deteccion de permutaciones (NUNCA modifica Xi)
  §3.1 K1 (formula) / K2 (autovalores de la caja) / K3-syn (control homogeneo sintetico, decisivo) /
       K3 real (solo diagnostico) / reconstruccion
  §3.2 barrido X0 a dx=0.02: invariantes (modulos + fases plaqueta), S_box, S_cont, proxy eps_S
  §3.3 barrido dx (X0=12,16; N=500/1000/2000): Richardson por elemento/Frobenius, S_ref, E_ref, conjunto M
  §4   integridad I-C0/I-C1/I-C4 (STOP sin clasificacion) -> ciencia C2/C3/C5 -> clasificacion A/B/C/D/D*
       + ADDENDUM_1: estados de fase NINGUNA/FALLO/OK/PARCIAL; C2 INCONCLUSO -> clase I (prevalece sobre C3/C5)
  §6   CSV en out/ (raíz del repo) + manifiesto JSON

Uso:  cd residual_anatomy/gates && python3 gate_P2_phase_reference.py [--overwrite]
Codigos de salida: 0 = completado (cualquier clasificacion); 2 = STOP por integridad; 1 = excepcion;
                   3 = rehusa sobrescribir artefactos existentes sin --overwrite.
"""
import os, sys, json, csv, time, hashlib, platform, resource, traceback
import numpy as np
import scipy
import scipy.linalg as la

HERE = os.path.dirname(os.path.abspath(__file__))
if HERE not in sys.path:
    sys.path.insert(0, HERE)
import gate_B2_1_algebraic as g1        # ensamblado auditado (solo lectura)
import gate_B2_2_conditioning as g2     # representacion fisica + extraccion de S (solo lectura)
import gate_B2_4_X0_convergence as g4   # estimador cond (sigma_max potencia / sigma_min inversa) (solo lectura)

# ----------------------------------------------------------------------------- parametros (§3, §7)
OMEGAS = (0.10, 0.05)
DX_X0SWEEP = 0.02
X0_LIST = (8.0, 10.0, 12.0, 14.0, 16.0, 20.0, 24.0, 28.0)
X_REF = 16.0
X0_PRE = (8.0, 10.0, 12.0, 14.0)          # tramo donde C2 exige decrecimiento monotono
X0_POST = (16.0, 20.0, 24.0, 28.0)        # tramo X0 >= 16
DX_SWEEP_X0 = (12.0, 16.0)
DX_SWEEP_N = (500, 1000, 2000)
K3_NPAIRS = 10
K3SYN_NCELLS = 20
K3_FILTER = 1e-3
PHASE_DEF = 1e-3                          # |S_ij| >= 1e-3 ||S||_F para interpretar fase / definir plaqueta / material
SEED_K3SYN = 20261116

THR = dict(
    eta_back=1e-12,          # I-C0
    mu_K12=1e-13,            # I-C1 |mu_K2-mu_K1|
    dxk_K12=1e-13,           # I-C1 dx|kbox_K2-kbox_K1|
    absmu=1e-13,             # I-C1 ||mu|-1| propagantes
    k3syn_factor=10.0,       # I-C1 K3-syn <= 10 eps_K3 (margen empirico)
    nyquist=np.pi/4,         # I-C0 |k| dx
    recon=1e-14,             # §3.1 reconstruccion (registrada; ver nota en informe)
    ident=1e-12,             # I-C4
    C2=1e-6, C3=1e-6,
    C3_ratio=100.0, C3_floor_factor=10.0,
    p_lo=1.8, p_hi=2.2, p_floor=1e-10,
    E_ref=5e-5,
    degen=1e-3,
    match_gap_warn=0.1,      # R1: aviso (NO decide) si best-second overlap < 0.1
    k3syn_min_samples=1,     # C1(audit): muestras validas minimas por (extremo, modo) en K3-syn
)
U = np.finfo(float).eps
OUT_DIR = os.path.join(os.path.dirname(os.path.dirname(HERE)), "out")
NPZ = g1.STATE
F_MODAL = os.path.join(OUT_DIR, "p2_modal_ref.csv")
F_KBOX = os.path.join(OUT_DIR, "p2_kbox_check.csv")
F_X0 = os.path.join(OUT_DIR, "p2_X0_sweep.csv")
F_INV = os.path.join(OUT_DIR, "p2_invariants.csv")
F_DX = os.path.join(OUT_DIR, "p2_dx_sweep.csv")
F_MAN = os.path.join(OUT_DIR, "p2_manifest.json")
F_DESIGN = os.path.join(os.path.dirname(HERE), "P2_DESIGN.md")     # diseno congelado (solo se lee su hash)
F_ADDENDUM = os.path.join(os.path.dirname(HERE), "P2_DESIGN_ADDENDUM_1.md")   # addendum aprobado (solo hash)
ALL_OUT = (F_MODAL, F_KBOX, F_X0, F_INV, F_DX, F_MAN)

SIG_IN = np.array([-1.0, +1.0, +1.0])     # sigma de canales IN  (L, R1, R2)
SIG_OUT = np.array([-1.0, +1.0, +1.0])    # sigma de canales OUT (L, R1, R2)
CH = ("L", "R1", "R2")


class IntegrityStop(Exception):
    """Fallo de integridad (I-C0/I-C1/I-C4): el gate se detiene SIN clasificacion."""


def _sha_file(p):
    with open(p, "rb") as f:
        return hashlib.sha256(f.read()).hexdigest()


def _sha_arr(a):
    return hashlib.sha256(np.ascontiguousarray(a).tobytes(order="C")).hexdigest()


def _wrap(a):
    return (np.asarray(a) + np.pi) % (2*np.pi) - np.pi


def safe_solve(A, b):
    """Resolucion protegida: (x, None) si entrada finita, solve sin excepcion y resultado finito;
       si no, (None, motivo). Nunca lanza: LinAlgError (singular) y ValueError (no finitos) se capturan."""
    if not (np.all(np.isfinite(A)) and np.all(np.isfinite(b))):
        return None, "entrada_no_finita"
    try:
        x = la.solve(A, b)
    except (la.LinAlgError, ValueError) as e:
        return None, f"solve_fallo_{type(e).__name__}"
    if not np.all(np.isfinite(x)):
        return None, "resultado_no_finito"
    return x, None


def _peak_rss_bytes():
    r = resource.getrusage(resource.RUSAGE_SELF).ru_maxrss
    return int(r) if sys.platform == "darwin" else int(r)*1024     # macOS: bytes; Linux: KiB


# ============================================================================ §1.1 base modal congelada
def modal_reference(w):
    """Base de referencia por omega (se calcula UNA vez). Guarda Xi, etiquetas, orden, k, hashes."""
    ref = {}
    for side in ("L", "R"):
        Xi, labels, nfix, cols = g1.build_Xi(side, w)
        ks = np.array([complex(m["k"]) for m in cols])
        sep = min(abs(ks[a]-ks[b]) for a in range(4) for b in range(4) if a != b)
        ref[side] = dict(Xi=Xi.copy(), labels=list(labels), ks=ks, prop=[bool(m["prop"]) for m in cols],
                         sha=_sha_arr(Xi), min_sep=float(sep), near_degenerate=bool(sep < THR["degen"]))
    labL, labR = ref["L"]["labels"], ref["R"]["labels"]
    iL_in, iL_out = labL.index("in"), labL.index("out")
    iR_in = [i for i, L in enumerate(labR) if L == "in"]       # orden = rhs_for (indice ascendente)
    iR_out = [i for i, L in enumerate(labR) if L == "out"]     # orden = g2.out_indices
    ref["k_in"] = np.array([ref["L"]["ks"][iL_in].real, ref["R"]["ks"][iR_in[0]].real, ref["R"]["ks"][iR_in[1]].real])
    ref["k_out"] = np.array([ref["L"]["ks"][iL_out].real, ref["R"]["ks"][iR_out[0]].real, ref["R"]["ks"][iR_out[1]].real])
    ref["col_in"] = [("L", iL_in), ("R", iR_in[0]), ("R", iR_in[1])]
    ref["col_out"] = [("L", iL_out), ("R", iR_out[0]), ("R", iR_out[1])]
    return ref


def check_basis(info, ref):
    """I-C0: Xi bit a bit + etiquetas; solapamiento SOLO para detectar permutaciones (no modifica Xi)."""
    out = {}
    for side, key, lab in (("L", "XiL", "labL"), ("R", "XiR", "labR")):
        Xi = info[key]
        bit = bool(np.array_equal(Xi, ref[side]["Xi"])) and (Xi.dtype == ref[side]["Xi"].dtype)
        labs = (list(info[lab]) == ref[side]["labels"])
        Xr = ref[side]["Xi"]
        O = np.abs(Xr.conj().T @ Xi) / np.outer(la.norm(Xr, axis=0), la.norm(Xi, axis=0))
        perm_ok = bool(np.array_equal(np.argmax(O, axis=0), np.arange(4)))
        out[side] = dict(bitwise=bit, labels=labs, perm_ok=perm_ok, overlap_diag_min=float(np.min(np.diag(O))))
    return out


# ============================================================================ §3.1 k_box: K1, K2, K3-syn, K3 real
def mu_K1(k, dx):
    k = complex(k)
    return (1 + 1j*k*dx/2) / (1 - 1j*k*dx/2), (2.0/dx)*np.arctan(k*dx/2)


def Ahom(side, w):
    return g1.A_homog(g1.nL, g1.vL, g1.VL, w) if side == "L" else g1.A_homog(g1.nR, g1.vR, g1.VR, w)


def k1k2_check(w, dx, ref, tag):
    """K1 vs K2 para los 4 modos de cada extremo. Comparar mu PRIMERO; luego rama de k_box hacia k."""
    rows = []; ok = True; qual = {}
    for side in ("L", "R"):
        A = Ahom(side, w); I4 = np.eye(4)
        T, why = safe_solve(I4 - dx/2*A, I4 + dx/2*A)
        if T is None:                                       # I-C1 falla de forma controlada
            ok = False; qual[side] = dict(error=f"T_caja_{why}")
            rows.append(dict(kind="K1K2", tag=tag, omega=w, dx=dx, side=side, pass_=False, reason=f"T_caja_{why}"))
            continue
        try:
            mus, V = la.eig(T)
        except (la.LinAlgError, ValueError) as e:
            ok = False; qual[side] = dict(error=f"eig_fallo_{type(e).__name__}")
            rows.append(dict(kind="K1K2", tag=tag, omega=w, dx=dx, side=side, pass_=False,
                             reason=f"eig_fallo_{type(e).__name__}"))
            continue
        Xi = ref[side]["Xi"]
        O = np.abs(V.conj().T @ Xi) / np.outer(la.norm(V, axis=0), la.norm(Xi, axis=0))   # [q, m]
        match = np.argmax(O, axis=0)                                                       # modo m -> autovector q
        bij = len(set(match.tolist())) == 4
        ok &= bij
        srt = np.sort(O, axis=0)
        gaps = srt[-1, :] - srt[-2, :]                                                     # best - second por modo
        qual[side] = dict(min_overlap=float(np.min(O[match, np.arange(4)])), min_gap=float(np.min(gaps)),
                          ambiguous=bool(np.min(gaps) < THR["match_gap_warn"]), bijective=bool(bij))
        for m in range(4):
            k = ref[side]["ks"][m]; prop = ref[side]["prop"][m]
            m1, kb1 = mu_K1(k, dx)
            m2 = mus[match[m]]
            dmu = abs(m2 - m1)
            base = (np.angle(m2) - 1j*np.log(abs(m2))) / dx          # -i log(mu)/dx, rama principal
            n = int(np.round((k.real - base.real) * dx / (2*np.pi)))  # rama mas proxima a k
            kb2 = base + 2*np.pi*n/dx
            dxk = dx*abs(kb2 - kb1)
            absmu = abs(abs(m2) - 1.0) if prop else np.nan
            nyq = abs(k)*dx
            r_ok = (dmu <= THR["mu_K12"]) and (dxk <= THR["dxk_K12"]) and (not prop or absmu <= THR["absmu"])
            nyq_ok = (not prop) or (nyq <= THR["nyquist"])
            ok &= bool(r_ok)          # Nyquist se evalua en I-C0, no aqui
            rows.append(dict(kind="K1K2", tag=tag, omega=w, dx=dx, side=side, mode=m, label=ref[side]["labels"][m],
                             prop=prop, k_re=k.real, k_im=k.imag, mu1_re=m1.real, mu1_im=m1.imag,
                             mu2_re=m2.real, mu2_im=m2.imag, dmu=dmu, kbox1_re=kb1.real, kbox1_im=kb1.imag,
                             kbox2_re=kb2.real, kbox2_im=kb2.imag, branch_n=n, dx_dkbox=dxk, absmu_K2=absmu,
                             kdx_over_pi=nyq/np.pi, nyquist_ok=nyq_ok, overlap=float(O[match[m], m]),
                             overlap_gap=float(gaps[m]), match_ambiguous=bool(gaps[m] < THR["match_gap_warn"]),
                             bijective=bij, pass_=bool(r_ok and nyq_ok)))
    return ok, rows, qual


def k3_ratios(Ys, Xi, mus1, labels, tag_common, kind, deltas=None):
    """Procedimiento K3 COMUN (real y sintetico): proyeccion c=Xi^{-1}Y, cociente c_m(j+1)/c_m(j),
       filtro |c_m(j)| >= 1e-3 max|c(j)|, eps_K3 = u cond(Xi) max|c|/|c_m|."""
    rows = []; worst = 0.0; counts = [0, 0, 0, 0]          # muestras VALIDAS por modo (cobertura)
    invalid = []                                           # FAIL-SAFE: cualquier no finito se registra aqui
    if not np.all(np.isfinite(Xi)):                        # sin excepciones: retorno inmediato con invalidez
        invalid.append(dict(pair=None, mode=None, reason="Xi_no_finita"))
        return rows, worst, counts, invalid
    try:
        condXi = np.linalg.cond(Xi)                        # SVD: puede no converger -> LinAlgError
    except (np.linalg.LinAlgError, ValueError):
        invalid.append(dict(pair=None, mode=None, reason="condXi_svd_fallo"))
        return rows, worst, counts, invalid                # filas vacias, cobertura cero, invalidez registrada
    C = []
    for j, Y in enumerate(Ys):                             # proyeccion protegida (nunca lanza)
        cv, why = safe_solve(Xi, Y)
        if cv is None:
            invalid.append(dict(pair=j, mode=None, reason=f"proyeccion_{why}"))
            C.append(np.full(4, np.nan, dtype=complex))    # marcado; cmax=NaN -> se salta abajo
        else:
            C.append(cv)
    if not (np.isfinite(condXi) and condXi > 0):
        invalid.append(dict(pair=None, mode=None, reason="condXi_invalido"))
    for j, cv in enumerate(C):                             # TODAS las proyecciones, incluida la ultima
        if not np.all(np.isfinite(cv)):
            invalid.append(dict(pair=j, mode=None, reason="proyeccion_no_finita"))
    for j in range(len(Ys)-1):
        cj, cj1 = C[j], C[j+1]
        cmax = np.max(np.abs(cj))
        if not np.isfinite(cmax):
            invalid.append(dict(pair=j, mode=None, reason="cmax_no_finito"))
            continue
        for m in range(4):
            if cmax == 0 or abs(cj[m]) < K3_FILTER*cmax:
                continue
            mu3 = cj1[m]/cj[m]
            err = abs(mu3 - mus1[m])
            eps = U*condXi*cmax/abs(cj[m])
            ratio = err/eps if (np.isfinite(eps) and eps > 0) else np.nan
            bad = [nm for nm, v in (("mu3", mu3), ("err", err), ("eps", eps), ("err/eps", ratio))
                   if not np.isfinite(v)]
            if not (np.isfinite(eps) and eps > 0):
                bad.append("eps<=0")
            if bad:                                        # NO se cuenta como muestra valida NI entra en worst
                invalid.append(dict(pair=j, mode=m, reason="no_finito:" + ",".join(sorted(set(bad)))))
                continue
            r = dict(tag_common); r.update(kind=kind, pair=j, mode=m, label=labels[m],
                                           mu3_re=mu3.real, mu3_im=mu3.imag, mu1_re=mus1[m].real,
                                           mu1_im=mus1[m].imag, err=err, eps_K3=eps, condXi=condXi,
                                           ratio_err_eps=err/eps)
            if deltas is not None:
                r.update(delta_10=deltas["d10"], delta_b=deltas["db"], ratio_err_dx_d10=(err/(r["dx"]*deltas["d10"])
                                                                                       if deltas["d10"] > 0 else np.nan))
            rows.append(r); worst = max(worst, ratio); counts[m] += 1      # ratio finito garantizado aqui
    return rows, worst, counts, invalid


def k3_syn(w, dx, ref, side, a, tag_common):
    """K3-syn: Y_0 = Xi a, MISMO esquema de caja con A_j == A_homog, 20 celdas, MISMO procedimiento K3."""
    A = Ahom(side, w); I4 = np.eye(4)
    Bm, Cp = I4 - dx/2*A, I4 + dx/2*A
    Xi = ref[side]["Xi"]
    pre = []                                               # invalideces de la propagacion sintetica
    if not (np.all(np.isfinite(a)) and np.all(np.isfinite(Xi))):
        return [], 0.0, [0, 0, 0, 0], [dict(pair=None, mode=None, reason="entrada_no_finita")]
    Ys = [Xi @ a]
    for j in range(K3SYN_NCELLS):
        y, why = safe_solve(Bm, Cp @ Ys[-1])               # Bm/rhs finitos, solve sin excepcion, y finito
        if y is None:
            pre.append(dict(pair=j, mode=None, reason=f"propagacion_{why}"))
            break                                          # -> secuencia incompleta -> retorno sin analisis
        Ys.append(y)
    if pre:                                                # secuencia INCOMPLETA: no se analiza como completa
        return [], 0.0, [0, 0, 0, 0], pre                  # filas vacias, cobertura cero, invalidez registrada
    mus1 = [mu_K1(k, dx)[0] for k in ref[side]["ks"]]
    rows, worst, counts, invalid = k3_ratios(Ys, Xi, mus1, ref[side]["labels"], tag_common, "K3syn")
    return rows, worst, counts, invalid


def k3_real(w, info, sols, ref):
    """K3 sobre la solucion REAL: SIEMPRE diagnostico. Ultimas 10 parejas junto a cada borde; por modo,
       la columna (problema) donde su amplitud en el borde es mayor."""
    N, dx, X0 = info["N"], info["dx"], info["X0"]
    xs = np.linspace(-X0, X0, N+1)
    rows = []; invalid = []                                # invalidos K3 real: SOLO registro (diagnostico)
    for side in ("L", "R"):
        js = list(range(0, K3_NPAIRS+1)) if side == "L" else list(range(N-K3_NPAIRS, N+1))
        jb = 0 if side == "L" else N
        Xi = ref[side]["Xi"]; Ah = Ahom(side, w)
        try:                                               # norma 2 matricial = SVD (puede no converger)
            dlist = [la.norm(g1.A_at(xs[j], w) - Ah, 2) for j in js]
            deltas = dict(d10=float(max(dlist)), db=float(la.norm(g1.A_at(xs[jb], w) - Ah, 2)))
        except (np.linalg.LinAlgError, ValueError) as e:  # K3 real = diagnostico: registrar y seguir
            deltas = dict(d10=np.nan, db=np.nan)
            invalid.append(dict(side=side, problem=None, target_mode=None, pair=None, mode=None,
                                reason=f"delta_svd_fallo_{type(e).__name__}"))
        mus1 = [mu_K1(k, dx)[0] for k in ref[side]["ks"]]
        cb = []                                            # proyeccion de borde por problema (None si no finita)
        for p_, (_z, x, _b, _M) in enumerate(sols):
            c_, why = safe_solve(Xi, x[4*jb:4*jb+4])      # aceptada SOLO si el resultado es finito
            cb.append(c_)
            if c_ is None:
                invalid.append(dict(side=side, problem=p_, target_mode=None, pair=None, mode=None,
                                    reason=f"proyeccion_borde_invalida_{why}"))
        valid_p = [p_ for p_, c in enumerate(cb) if c is not None]
        if not valid_p:
            continue                                       # ningun problema proyectable en este extremo
        for m in range(4):
            p = max(valid_p, key=lambda q: abs(cb[q][m]))
            x = sols[p][1]
            Ys = [x[4*j:4*j+4] for j in js]
            tc = dict(omega=w, X0=X0, N=N, dx=dx, side=side, problem=p)
            rr, _, _, inv = k3_ratios(Ys, Xi, mus1, ref[side]["labels"], tc, "K3real", deltas)   # diagnostico
            rows += [r for r in rr if r["mode"] == m]
            invalid += [dict(side=side, problem=p, target_mode=m, **d) for d in inv]
    return rows, invalid


# ============================================================================ ejecucion de un caso (w, X0, N)
def run_case(w, X0, N, ref, tag):
    t0 = time.time()
    M, info = g1.assemble(w, N, X0)
    dx = info["dx"]
    labL, labR = info["labL"], info["labR"]
    window = (labL.count("in") == 1 and labL.count("out") == 1 and labL.count("dec") == 1 and
              labL.count("grow") == 1 and labR.count("in") == 2 and labR.count("out") == 2 and
              labR.count("dec") == 0 and labR.count("grow") == 0)
    basis = check_basis(info, ref)
    dr, dc = g2.physical_scaling(info)
    S_ext, sols, outs, Mrep = g2.smatrix_and_solves(info, M, dr, dc)   # ruta primal auditada (la.solve)
    del M
    est = g4.estim_metrics(Mrep, info["ndof"], lu=None, strict=False)  # sigma_min via lu_solve(trans=2) (P1)
    Mnorm = est["sigma_max"]
    etas, rabs = [], []
    for (z, x, brep, _M) in sols:
        r = la.norm(Mrep @ z - brep)
        rabs.append(r); etas.append(r/(Mnorm*la.norm(z) + la.norm(brep)))
    eta = float(max(etas))
    xnorm = float(max(la.norm(s[1]) for s in sols))
    i_dec = info["iaL"] + labL.index("dec")
    a_dec = float(max(abs(s[1][i_dec]) for s in sols))
    # re-referencias (§2)
    k_in, k_out = ref["k_in"], ref["k_out"]
    kb_in = np.array([mu_K1(k, dx)[1].real for k in k_in]); kb_out = np.array([mu_K1(k, dx)[1].real for k in k_out])
    Dout_b, Din_b = np.exp(1j*kb_out*SIG_OUT*X0), np.exp(1j*kb_in*SIG_IN*X0)
    Dout_c, Din_c = np.exp(1j*k_out*SIG_OUT*X0), np.exp(1j*k_in*SIG_IN*X0)
    S_box = (1.0/Dout_b)[:, None]*S_ext*Din_b[None, :]
    S_cont = (1.0/Dout_c)[:, None]*S_ext*Din_c[None, :]
    recon = float(np.max(np.abs(Dout_b[:, None]*S_box*(1.0/Din_b)[None, :] - S_ext)))
    Phi = X0*(((kb_out-k_out)*SIG_OUT)[:, None] - ((kb_in-k_in)*SIG_IN)[None, :])
    ident = np.abs(S_cont - S_box*np.exp(1j*Phi)) / np.maximum(1.0, np.abs(S_box))
    ident_max = float(np.max(ident))
    # K1/K2, K3-syn, K3 real
    okK12, rows_k12, match_qual = k1k2_check(w, dx, ref, tag)
    rows_syn = []; worst_syn = 0.0
    cov_syn = {sd: [0, 0, 0, 0] for sd in ("L", "R")}            # cobertura total por (extremo, modo)
    cov_syn_avec = {}                                              # desglose por vector a
    invalid_syn = []                                               # no finitos en K3-syn -> I-C1 falla
    rng = np.random.default_rng([SEED_K3SYN, int(round(w*1000)), N, int(round(X0*10))])
    for side in ("L", "R"):
        jb = 0 if side == "L" else info["N"]
        avecs = [("random", rng.standard_normal(4) + 1j*rng.standard_normal(4))]
        for p, s in enumerate(sols):
            a_, why = safe_solve(ref[side]["Xi"], s[1][4*jb:4*jb+4])   # Xi/yb finitos, solve ok, a finito
            if a_ is not None:
                avecs.append((f"real_p{p}", a_))
            else:                                          # entrada de K3-syn invalida -> invalidez (I-C1)
                invalid_syn.append(dict(side=side, avec=f"real_p{p}", pair=None, mode=None,
                                        reason=f"vector_a_invalido_{why}"))
        for name, a in avecs:
            tc = dict(omega=w, X0=X0, N=N, dx=dx, side=side, avec=name)
            rr, wv, cnt, inv = k3_syn(w, dx, ref, side, a, tc)
            rows_syn += rr; worst_syn = max(worst_syn, wv)          # wv finito (k3_ratios lo garantiza)
            invalid_syn += [dict(side=side, avec=name, **d) for d in inv]
            cov_syn_avec[f"{side}:{name}"] = cnt
            for m in range(4):
                cov_syn[side][m] += cnt[m]
    cov_ok = all(cov_syn[sd][m] >= THR["k3syn_min_samples"] for sd in ("L", "R") for m in range(4))
    rows_k3, invalid_k3real = k3_real(w, info, sols, ref)        # diagnostico: no entra en integridad
    cov_real = {sd: [sum(1 for r in rows_k3 if r["side"] == sd and r["mode"] == m) for m in range(4)]
                for sd in ("L", "R")}                              # K3 real: solo diagnostico
    # integridad
    # Nyquist (I-C0) calculado DIRECTAMENTE de la base congelada (no depende de las filas de K2, que pueden
    # ser filas de error sin esa clave si K2 falla)
    nyq_ok = all(abs(ref[sd]["ks"][m])*dx <= THR["nyquist"]
                 for sd in ("L", "R") for m in range(4) if ref[sd]["prop"][m])
    ic0 = dict(window=window, basis=basis, nyquist=nyq_ok, eta_back=eta,
               ok=bool(window and nyq_ok and eta <= THR["eta_back"] and
                       all(b["bitwise"] and b["labels"] and b["perm_ok"] for b in basis.values())))
    # I-C1: PRIMERO cobertura (todos los modos medidos), DESPUES el peor error. Sin cobertura -> falla.
    # I-C1 FAIL-SAFE: (1) K1/K2; (2) cobertura de los 4 modos x 2 extremos; (3) NINGUN valor no finito en
    # K3-syn; (4) solo entonces el peor cociente medido. Cualquier invalidez -> I-C1 = False (STOP).
    syn_valid = (len(invalid_syn) == 0)
    ic1 = dict(K1K2=okK12, K3syn_coverage=cov_syn, K3syn_coverage_ok=cov_ok,
               K3syn_valid=syn_valid, K3syn_invalid=invalid_syn,
               K3syn_worst_ratio=(worst_syn if (cov_ok and syn_valid) else None),
               match_quality_no_decide=match_qual,
               ok=bool(okK12 and cov_ok and syn_valid and np.isfinite(worst_syn)
                       and worst_syn <= THR["k3syn_factor"]))
    ic4 = dict(ident_max=ident_max, ok=bool(ident_max <= THR["ident"]))
    # proxy heuristico: NaN (no fiable) si el estimador de cond no convergio o es invalido
    eps_proxy = (float(est["cond"]*eta*np.max(np.abs(S_ext)))
                 if (est["converged"] and np.isfinite(est["cond"])) else np.nan)
    return dict(tag=tag, omega=w, X0=X0, N=N, dx=dx, S_ext=S_ext, S_box=S_box, S_cont=S_cont, Phi=Phi,
                recon=recon, ident=ident, eta_back=eta, resid_abs=float(max(rabs)), xnorm=xnorm, a_dec=a_dec,
                cond_est=float(est["cond"]), sigma_min=float(est["sigma_min"]), sigma_max=float(est["sigma_max"]),
                est_converged=bool(est["converged"]), maxS=float(np.max(np.abs(S_ext))), eps_proxy=eps_proxy,
                kb_in=kb_in, kb_out=kb_out, ic0=ic0, ic1=ic1, ic4=ic4,
                k3syn_cov_avec=cov_syn_avec, k3real_cov=cov_real, k3real_invalid=invalid_k3real,
                rows_k12=rows_k12, rows_syn=rows_syn, rows_k3=rows_k3, seconds=time.time()-t0)


def _cov_entry(c):
    """R2: cobertura efectiva de K3 por caso + calidad del emparejamiento K1/K2 (no decide)."""
    return dict(omega=c["omega"], tag=c["tag"], X0=c["X0"], N=c["N"],
                k3syn_total=c["ic1"]["K3syn_coverage"], k3syn_coverage_ok=c["ic1"]["K3syn_coverage_ok"],
                k3syn_valid=c["ic1"]["K3syn_valid"], k3syn_invalid=c["ic1"]["K3syn_invalid"],
                k3syn_by_avec=c["k3syn_cov_avec"], k3real_diag=c["k3real_cov"],
                k3real_invalid_diag_no_decide=c["k3real_invalid"],
                k1k2_match_quality=c["ic1"]["match_quality_no_decide"], integrity_ok=dict(
                    ic0=c["ic0"]["ok"], ic1=c["ic1"]["ok"], ic4=c["ic4"]["ok"]))


def integrity_or_stop(c):
    fails = [n for n in ("ic0", "ic1", "ic4") if not c[n]["ok"]]
    if fails:
        raise IntegrityStop(f"omega={c['omega']} X0={c['X0']} N={c['N']}: falla {fails} -> "
                            f"{ {n: {k: v for k, v in c[n].items() if k != 'basis'} for n in fails} }"
                            + (f" basis={c['ic0']['basis']}" if "ic0" in fails else ""))


# ============================================================================ metricas cientificas
def E_X(Sa, Sb):
    return float(np.max(np.abs(Sa - Sb)))


def plaquettes(S):
    """I_ij = arg(S00 Sij S0j* Si0*), i,j in {1,2}; None si algun elemento < 1e-3 ||S||_F."""
    thr = PHASE_DEF*la.norm(S)
    out = {}
    for i in (1, 2):
        for j in (1, 2):
            els = (S[0, 0], S[i, j], S[0, j], S[i, 0])
            out[(i, j)] = (float(np.angle(S[0, 0]*S[i, j]*np.conj(S[0, j])*np.conj(S[i, 0])))
                           if min(abs(e) for e in els) >= thr else None)
    return out


def E_inv(Pa, Pb):
    vals = [abs(float(_wrap(Pa[q]-Pb[q]))) for q in Pa if Pa[q] is not None and Pb[q] is not None]
    return (max(vals) if vals else None), len(vals)


def richardson(Sh, Sh2, Sh4):
    d1, d2 = np.abs(Sh-Sh2), np.abs(Sh2-Sh4)
    est = (d1 >= THR["p_floor"]) & (d2 >= THR["p_floor"])
    p = np.full((3, 3), np.nan)
    p[est] = np.log2(d1[est]/d2[est])
    n1, n2 = la.norm(Sh-Sh2), la.norm(Sh2-Sh4)
    pF = float(np.log2(n1/n2)) if (n1 >= THR["p_floor"] and n2 >= THR["p_floor"]) else np.nan
    Sref = (4.0*Sh4 - Sh2)/3.0
    return p, est, pF, Sref


# ============================================================================ evaluacion por omega
def evaluate_omega(w, cX0, cDX):
    ref = cX0[X_REF]
    # ---- C2 (invariantes) ----
    PR = plaquettes(ref["S_ext"])
    E_mod = {X0: E_X(np.abs(c["S_ext"]), np.abs(ref["S_ext"])) for X0, c in cX0.items()}
    P = {X0: plaquettes(c["S_ext"]) for X0, c in cX0.items()}
    E_I = {}; nI = {}
    for X0 in cX0:
        E_I[X0], nI[X0] = E_inv(P[X0], PR)                       # informativo (conjunto variable por X0)
    post_mod = all(E_mod[X] <= THR["C2"] for X in X0_POST)
    mono_mod = all(E_mod[a] > E_mod[b] for a, b in zip(X0_PRE[:-1], X0_PRE[1:]))   # monotonia ESTRICTA
    # Fases: conjunto FIJO de plaquetas por tramo (definidas en X_ref y en TODOS los X0 del tramo), para
    # comparar siempre las mismas cantidades. Tres estados: OK / FALLO / NO_EVALUABLE (conjunto vacio).
    def qset(Xs):
        return [q for q in PR if PR[q] is not None and all(P[X][q] is not None for X in Xs)]
    def Eq(X, Q):
        return max(abs(float(_wrap(P[X][q] - PR[q]))) for q in Q)
    Q_post, Q_pre = qset(X0_POST), qset(X0_PRE)
    E_Ipost = {X: Eq(X, Q_post) for X in X0_POST} if Q_post else {}
    E_Ipre = {X: Eq(X, Q_pre) for X in X0_PRE} if Q_pre else {}
    ph_post = ("VACIO" if not Q_post else
               ("OK" if all(E_Ipost[X] <= THR["C2"] for X in X0_POST) else "FALLO"))
    ph_pre = ("VACIO" if not Q_pre else
              ("OK" if all(E_Ipre[a] > E_Ipre[b] for a, b in zip(X0_PRE[:-1], X0_PRE[1:])) else "FALLO"))
    # Estado de fases segun P2_DESIGN_ADDENDUM_1.md A2 (mutuamente excluyentes, en este orden):
    any_defined = any(P[X][q] is not None for X in P for q in P[X])          # incluye X_ref (esta en cX0)
    if not any_defined:
        phases = "NINGUNA"                       # exactamente "Sin I_ij definidas" del §3.2(i)
    elif "FALLO" in (ph_post, ph_pre):
        phases = "FALLO"                         # prioritario aunque el otro tramo este vacio
    elif ph_post == ph_pre == "OK":
        phases = "OK"
    else:
        phases = "PARCIAL"                       # hay plaquetas definidas pero no se evaluan ambos tramos
    # Estado de C2 (addendum A3)
    if not (post_mod and mono_mod) or phases == "FALLO":
        C2_state = "FALLA"
    elif phases == "OK":
        C2_state = "PASA"
    elif phases == "NINGUNA":
        C2_state = "PASA_SOLO_MODULOS"           # §3.2(i), con anotacion obligatoria
    else:
        C2_state = "INCONCLUSO"                  # PARCIAL -> clase I (addendum A4)
    C2 = C2_state in ("PASA", "PASA_SOLO_MODULOS")
    xs_ = np.array(X0_PRE); ys_ = np.array([E_mod[X] for X in X0_PRE])
    try:                                                   # tasa SOLO informativa: nunca debe abortar el gate
        rate = (float(-np.polyfit(xs_, np.log(ys_), 1)[0])
                if (np.all(np.isfinite(ys_)) and np.all(ys_ > 0)) else np.nan)
    except (np.linalg.LinAlgError, ValueError):
        rate = np.nan
    # ---- C3 ----
    E_box = {X0: E_X(c["S_box"], ref["S_box"]) for X0, c in cX0.items()}
    E_cont = {X0: E_X(c["S_cont"], ref["S_cont"]) for X0, c in cX0.items()}
    C3 = all(E_box[X] <= THR["C3"] for X in X0_POST)
    # Ambos proxies deben ser finitos INDIVIDUALMENTE (max() con NaN depende del orden de los argumentos).
    e28, eref = cX0[28.0]["eps_proxy"], cX0[X_REF]["eps_proxy"]
    eps_pair = max(e28, eref) if (np.isfinite(e28) and np.isfinite(eref)) else np.nan
    if not np.isfinite(eps_pair):
        ratio, ratio_status = np.nan, "no_concluyente_proxy_invalido"
    elif not (np.isfinite(E_box[28.0]) and E_box[28.0] > 0):            # evita division por cero / invalida
        ratio, ratio_status = np.nan, "no_concluyente_Ebox_no_positivo"
    elif not np.isfinite(E_cont[28.0]):                                   # numerador invalido -> no concluyente
        ratio, ratio_status = np.nan, "no_concluyente_Econt_no_finito"
    elif E_box[28.0] >= THR["C3_floor_factor"]*eps_pair:
        ratio = E_cont[28.0]/E_box[28.0]
        ratio_status = "indicativo_ok" if ratio >= THR["C3_ratio"] else "indicativo_bajo"
    else:
        ratio, ratio_status = np.nan, "no_concluyente_bajo_proxy"
    # ---- C5 ----
    dx_res = {}
    for X0 in DX_SWEEP_X0:
        cs = [cDX[(X0, N)] for N in DX_SWEEP_N]
        res = {}
        for fam in ("S_box", "S_cont"):
            p, est, pF, Sref = richardson(cs[0][fam], cs[1][fam], cs[2][fam])
            res[fam] = dict(p=p, est=est, pF=pF, Sref=Sref)
        Sfine = cs[2]["S_box"]                                           # M: S_fis^box en N=2000 (§3.3)
        Mset = np.abs(Sfine) >= PHASE_DEF*la.norm(Sfine)
        pb, eb = res["S_box"]["p"], res["S_box"]["est"]
        n_mat = int(np.sum(Mset))
        mat_est = bool(n_mat > 0 and np.all(eb[Mset]))             # M vacio -> no estimable -> INCONCLUSO
        mat_inrange = bool(np.all((pb[Mset & eb] >= THR["p_lo"]) & (pb[Mset & eb] <= THR["p_hi"])))
        E_bc = E_X(res["S_box"]["Sref"], res["S_cont"]["Sref"])
        dx_res[X0] = dict(res=res, Mset=Mset, n_material=n_mat, mat_est=mat_est, mat_inrange=mat_inrange,
                          E_ref_box_cont=E_bc)
    E_12_16 = E_X(dx_res[12.0]["res"]["S_box"]["Sref"], dx_res[16.0]["res"]["S_box"]["Sref"])
    E_12_16_cont = E_X(dx_res[12.0]["res"]["S_cont"]["Sref"], dx_res[16.0]["res"]["S_cont"]["Sref"])
    any_nonest = any(not d["mat_est"] for d in dx_res.values())
    c5_fail = (any(not d["mat_inrange"] for d in dx_res.values()) or
               any(d["E_ref_box_cont"] > THR["E_ref"] for d in dx_res.values()) or E_12_16 > THR["E_ref"])
    C5 = "FALLA" if c5_fail else ("INCONCLUSO" if any_nonest else "PASA")
    # ---- clasificacion (§4 + addendum A4): C2 FALLA -> C; C2 INCONCLUSO -> I (prevalece sobre C3/C5) ----
    if C2_state == "FALLA":
        cls = "C"
    elif C2_state == "INCONCLUSO":
        cls = "I"
    elif not C3:
        cls = "B"
    elif C5 == "PASA":
        cls = "A"
    elif C5 == "FALLA":
        cls = "D"
    else:
        cls = "D*"
    if C2_state == "PASA_SOLO_MODULOS":
        cls_note = "C2=PASA_SOLO_MODULOS: las fases no aportan evidencia favorable (§3.2(i))"
    elif C2_state == "INCONCLUSO":
        cls_note = "clase I: inconcluso por fases parcialmente evaluables (addendum A4); C3/C5 solo registrados"
    else:
        cls_note = ""
    return dict(omega=w, C2=C2, C2_state=C2_state, classification_note=cls_note,
                C2_detail=dict(post_mod=post_mod, mono_mod_strict=mono_mod, phases=phases, phases_post=ph_post,
                               phases_pre_strict_monotone=ph_pre, any_plaquette_defined=any_defined, Q_post=[f"I_{i}{j}" for i, j in Q_post],
                               Q_pre=[f"I_{i}{j}" for i, j in Q_pre], E_I_post=E_Ipost, E_I_pre=E_Ipre,
                               rate_informativa=rate),
                est_converged_28_16=bool(cX0[28.0]["est_converged"] and cX0[X_REF]["est_converged"]),
                E_mod=E_mod, E_I=E_I, nI=nI, C3=bool(C3), E_box=E_box, E_cont=E_cont, eps_pair=eps_pair,
                C3_ratio=ratio, C3_ratio_status=ratio_status, C5=C5, dx_res=dx_res, E_ref_12_16=E_12_16,
                E_ref_12_16_cont=E_12_16_cont, classification=cls)


# ============================================================================ persistencia (§6)
def write_csv(path, rows):
    keys = []
    for r in rows:
        for k in r:
            if k not in keys:
                keys.append(k)
    with open(path, "w", newline="") as f:
        wtr = csv.DictWriter(f, fieldnames=keys)
        wtr.writeheader()
        for r in rows:
            wtr.writerow(r)


def rows_modal(refs):
    rows = []
    for w, ref in refs.items():
        for side in ("L", "R"):
            R = ref[side]
            for m in range(4):
                r = dict(omega=w, side=side, col=m, label=R["labels"][m], prop=R["prop"][m],
                         k_re=R["ks"][m].real, k_im=R["ks"][m].imag, min_sep=R["min_sep"],
                         near_degenerate=R["near_degenerate"], Xi_sha256=R["sha"])
                for c in range(4):
                    r[f"xi{c}_re"] = R["Xi"][c, m].real; r[f"xi{c}_im"] = R["Xi"][c, m].imag
                rows.append(r)
    return rows


def rows_X0(w, cX0):
    ref = cX0[X_REF]; rows = []
    for X0, c in sorted(cX0.items()):
        nF = la.norm(c["S_ext"])
        for o in range(3):
            for p in range(3):
                se, sb, sc = c["S_ext"][o, p], c["S_box"][o, p], c["S_cont"][o, p]
                rb, rc = ref["S_box"][o, p], ref["S_cont"][o, p]
                pdef = bool(abs(se) >= PHASE_DEF*nF) and bool(abs(ref["S_ext"][o, p]) >= PHASE_DEF*la.norm(ref["S_ext"]))
                rows.append(dict(
                    omega=w, X0=X0, N=c["N"], dx=c["dx"], out=CH[o], inp=CH[p],
                    Sext_re=se.real, Sext_im=se.imag, Sbox_re=sb.real, Sbox_im=sb.imag,
                    Scont_re=sc.real, Scont_im=sc.imag, absS=abs(se), argS_ext=float(np.angle(se)),
                    phase_defined=pdef,
                    dS_box=abs(sb-rb), dS_cont=abs(sc-rc),
                    dlogabs=float(np.log(abs(se)/abs(ref["S_ext"][o, p]))) if abs(se) > 0 and abs(ref["S_ext"][o, p]) > 0 else np.nan,
                    darg_box=float(_wrap(np.angle(sb)-np.angle(rb))) if pdef else np.nan,
                    darg_cont=float(_wrap(np.angle(sc)-np.angle(rc))) if pdef else np.nan,
                    Phi=c["Phi"][o, p], ident_err=c["ident"][o, p], recon=c["recon"],
                    eta_back=c["eta_back"], resid_abs=c["resid_abs"], cond_est=c["cond_est"],
                    est_converged=c["est_converged"], maxS=c["maxS"], eps_proxy=c["eps_proxy"],
                    xnorm=c["xnorm"], a_dec=c["a_dec"], seconds=c["seconds"]))
    return rows


def rows_X0_raw(w, cases):
    """Volcado CRUDO (sin metricas respecto a X_ref) de los casos calculados antes de un STOP/excepcion."""
    rows = []
    for c in cases:
        for o in range(3):
            for p in range(3):
                se, sb, sc = c["S_ext"][o, p], c["S_box"][o, p], c["S_cont"][o, p]
                rows.append(dict(omega=w, X0=c["X0"], N=c["N"], dx=c["dx"], tag=c["tag"], out=CH[o], inp=CH[p],
                                 Sext_re=se.real, Sext_im=se.imag, Sbox_re=sb.real, Sbox_im=sb.imag,
                                 Scont_re=sc.real, Scont_im=sc.imag, ident_err=c["ident"][o, p], recon=c["recon"],
                                 eta_back=c["eta_back"], resid_abs=c["resid_abs"], cond_est=c["cond_est"],
                                 eps_proxy=c["eps_proxy"], ic0=c["ic0"]["ok"], ic1=c["ic1"]["ok"], ic4=c["ic4"]["ok"],
                                 partial_dump=True))
    return rows


def rows_inv(w, cX0, ev):
    rows = []
    for X0, c in sorted(cX0.items()):
        P = plaquettes(c["S_ext"])
        dd = ev["C2_detail"]
        r = dict(omega=w, X0=X0, E_mod=ev["E_mod"][X0], E_inv_informativo=ev["E_I"][X0], n_inv_defined=ev["nI"][X0],
                 E_inv_Qpost=dd["E_I_post"].get(X0), E_inv_Qpre=dd["E_I_pre"].get(X0), C2_state=ev["C2_state"],
                 E_box=ev["E_box"][X0], E_cont=ev["E_cont"][X0])
        for (i, j), v in P.items():
            r[f"I_{i}{j}"] = v; r[f"I_{i}{j}_defined"] = v is not None
        for o in range(3):
            for p in range(3):
                r[f"abs_{CH[o]}_{CH[p]}"] = abs(c["S_ext"][o, p])
        rows.append(r)
    return rows


def rows_dx(w, cDX, ev):
    rows = []
    for X0 in DX_SWEEP_X0:
        d = ev["dx_res"][X0]
        for N in DX_SWEEP_N:
            c = cDX[(X0, N)]
            for o in range(3):
                for p in range(3):
                    r = dict(omega=w, X0=X0, N=N, dx=c["dx"], out=CH[o], inp=CH[p],
                             Sbox_re=c["S_box"][o, p].real, Sbox_im=c["S_box"][o, p].imag,
                             Scont_re=c["S_cont"][o, p].real, Scont_im=c["S_cont"][o, p].imag,
                             eta_back=c["eta_back"], material=bool(d["Mset"][o, p]))
                    for fam, nm in (("S_box", "box"), ("S_cont", "cont")):
                        R = d["res"][fam]
                        r[f"p_{nm}"] = R["p"][o, p] if R["est"][o, p] else "no_estimable"
                        r[f"pFrob_{nm}"] = R["pF"]
                        r[f"Sref_{nm}_re"] = R["Sref"][o, p].real; r[f"Sref_{nm}_im"] = R["Sref"][o, p].imag
                    r["E_ref_box_cont"] = d["E_ref_box_cont"]; r["E_ref_12_16_box"] = ev["E_ref_12_16"]
                    r["E_ref_12_16_cont"] = ev["E_ref_12_16_cont"]
                    rows.append(r)
    return rows


def manifest(status, exit_code, t_start, evals, refs, note, err=None, cov_log=None):
    try:
        blas = str(np.__config__.CONFIG.get("Build Dependencies", {}).get("blas", {}))
    except Exception:
        blas = "(no disponible)"
    files = {}
    for fn in ("gate_P2_phase_reference.py", "gate_B2_1_algebraic.py", "gate_B2_2_conditioning.py",
               "gate_B2_4_X0_convergence.py"):
        files[fn] = _sha_file(os.path.join(HERE, fn))
    outs = {os.path.basename(p): _sha_file(p) for p in ALL_OUT if p != F_MAN and os.path.exists(p)}
    def jsonable(v):
        if isinstance(v, dict):
            return {str(k): jsonable(x) for k, x in v.items()}
        if isinstance(v, (list, tuple)):
            return [jsonable(x) for x in v]
        if isinstance(v, (np.bool_,)):
            return bool(v)
        if isinstance(v, (np.integer,)):
            return int(v)
        if isinstance(v, (float, np.floating)):
            return None if not np.isfinite(v) else float(v)
        if isinstance(v, np.ndarray):
            return jsonable(v.tolist())
        return v
    summ = {}
    for w, ev in evals.items():
        summ[str(w)] = dict(classification=ev["classification"], classification_note=ev["classification_note"],
                            C2=ev["C2"], C2_state=ev["C2_state"], C2_detail=ev["C2_detail"],
                            est_converged_28_16=ev["est_converged_28_16"],
                            C3=ev["C3"], C3_ratio=ev["C3_ratio"], C3_ratio_status=ev["C3_ratio_status"],
                            eps_pair=ev["eps_pair"], C5=ev["C5"], E_ref_12_16_box=ev["E_ref_12_16"],
                            C5_detail={str(X0): dict(n_material=d["n_material"], material_estimable=d["mat_est"],
                                                     material_in_range=d["mat_inrange"],
                                                     E_ref_box_cont=d["E_ref_box_cont"],
                                                     excluded_nonmaterial=[f"{CH[o]}<-{CH[p]}" for o in range(3) for p in range(3)
                                                                           if not d["Mset"][o, p]])
                                       for X0, d in ev["dx_res"].items()})
    man = dict(
        gate="P2 phase reference",
        design="residual_anatomy/P2_DESIGN.md (congelado 2026-11-16) + P2_DESIGN_ADDENDUM_1.md (C2 parcial -> clase I)",
        timestamp=time.strftime("%Y-%m-%d %H:%M:%S %Z"), command=" ".join([sys.executable] + sys.argv),
        cwd=os.getcwd(), exit_code=exit_code, integrity=status, note=note, error=err,
        seconds_total=time.time()-t_start, peak_rss_bytes=_peak_rss_bytes(),
        environment=dict(python=platform.python_version(), numpy=np.__version__, scipy=scipy.__version__,
                         platform=platform.platform(), blas=blas),
        sha256=dict(scripts=files, design_P2_DESIGN_md=_sha_file(F_DESIGN),
                    design_P2_DESIGN_ADDENDUM_1_md=_sha_file(F_ADDENDUM), input_npz=_sha_file(NPZ), outputs=outs,
                    Xi_ref={str(w): {s: r[s]["sha"] for s in ("L", "R")} for w, r in refs.items()}),
        parameters=dict(omegas=OMEGAS, dx_X0sweep=DX_X0SWEEP, X0_list=X0_LIST, X_ref=X_REF, X0_pre=X0_PRE,
                        X0_post=X0_POST, dx_sweep_X0=DX_SWEEP_X0, dx_sweep_N=DX_SWEEP_N, k3_npairs=K3_NPAIRS,
                        k3syn_ncells=K3SYN_NCELLS, k3_filter=K3_FILTER, phase_def=PHASE_DEF, seed_k3syn=SEED_K3SYN,
                        sigma_in=SIG_IN, sigma_out=SIG_OUT, eta_back_repr="representacion fisica B2.2, ||M|| = sigma_max (potencia)"),
        thresholds=THR, results=summ,
        k3_coverage=cov_log or [],          # R2: cobertura efectiva K3-syn (decide) y K3 real (diagnostico) por caso
        global_classification=(None if status != "PASS" else
                               ("A" if all(ev["classification"] == "A" for ev in evals.values())
                                else "NO-A (ver por omega)")),
    )
    with open(F_MAN, "w") as f:
        json.dump(jsonable(man), f, indent=2, ensure_ascii=False)


# ============================================================================ main
def main():
    t_start = time.time()
    if any(os.path.exists(p) for p in ALL_OUT) and "--overwrite" not in sys.argv:
        print("Artefactos P2 ya existen en out/; usar --overwrite para regenerarlos. No se hace nada.")
        return 3
    os.makedirs(OUT_DIR, exist_ok=True)
    refs, evals, partial = {}, {}, {}
    cov_log = []
    rows_k = []; rows_x0 = []; rows_i = []; rows_d = []
    status, exit_code, note, err = "PASS", 0, "", None
    try:
        for w in OMEGAS:
            refs[w] = modal_reference(w)
        write_csv(F_MODAL, rows_modal(refs))
        for w in OMEGAS:
            ref = refs[w]
            cX0, cDX = {}, {}
            for X0 in X0_LIST:
                N = int(round(2*X0/DX_X0SWEEP))
                c = run_case(w, X0, N, ref, tag="X0sweep")
                rows_k += c["rows_k12"] + c["rows_syn"] + c["rows_k3"]
                cov_log.append(_cov_entry(c))
                print(f"  w={w} X0={X0:4.1f} N={N:5d}  I-C0={c['ic0']['ok']} I-C1={c['ic1']['ok']} "
                      f"I-C4={c['ic4']['ok']}  eta={c['eta_back']:.1e}  {c['seconds']:.0f}s", flush=True)
                for k in ("rows_k12", "rows_syn", "rows_k3"):
                    c[k] = None
                cX0[X0] = c; partial[w] = (cX0, cDX)
                integrity_or_stop(c)
            for X0 in DX_SWEEP_X0:
                for N in DX_SWEEP_N:
                    c = run_case(w, X0, N, ref, tag="dxsweep")
                    rows_k += c["rows_k12"] + c["rows_syn"] + c["rows_k3"]
                    cov_log.append(_cov_entry(c))
                    print(f"  w={w} X0={X0:4.1f} N={N:5d} (dx)  I-C0={c['ic0']['ok']} I-C1={c['ic1']['ok']} "
                          f"I-C4={c['ic4']['ok']}  eta={c['eta_back']:.1e}  {c['seconds']:.0f}s", flush=True)
                    for k in ("rows_k12", "rows_syn", "rows_k3"):
                        c[k] = None
                    cDX[(X0, N)] = c; partial[w] = (cX0, cDX)
                    integrity_or_stop(c)
            ev = evaluate_omega(w, cX0, cDX)
            evals[w] = ev
            rows_x0 += rows_X0(w, cX0); rows_i += rows_inv(w, cX0, ev); rows_d += rows_dx(w, cDX, ev)
            print(f"  w={w}: C2={ev['C2_state']} C3={ev['C3']} (ratio {ev['C3_ratio_status']}) C5={ev['C5']} "
                  f"-> clasificacion {ev['classification']} {ev['classification_note']}", flush=True)
    except IntegrityStop as e:
        status, exit_code, note = "STOP", 2, str(e)
        print(f"\n*** STOP POR INTEGRIDAD (sin clasificacion): {e}", flush=True)
    except Exception as e:
        status, exit_code, err = "EXCEPTION", 1, traceback.format_exc()
        print(err, flush=True)
    finally:
        # se persiste TODO lo disponible, pase o falle (§4: "todo resultado se persiste")
        for w, (pX0, pDX) in partial.items():            # omegas no evaluadas: volcado crudo sin metricas
            if w not in evals:
                rows_x0 += rows_X0_raw(w, list(pX0.values()) + list(pDX.values()))
        if rows_k: write_csv(F_KBOX, rows_k)
        if rows_x0: write_csv(F_X0, rows_x0)
        if rows_i: write_csv(F_INV, rows_i)
        if rows_d: write_csv(F_DX, rows_d)
        if status == "PASS" and not evals:
            status = "INCOMPLETE"
        manifest(status, exit_code, t_start, evals, refs, note, err, cov_log)
        print(f"\n  manifiesto -> {F_MAN}  (integridad={status}, exit={exit_code})", flush=True)
    return exit_code


if __name__ == "__main__":
    sys.exit(main())
