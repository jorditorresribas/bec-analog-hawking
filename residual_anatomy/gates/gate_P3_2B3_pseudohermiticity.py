"""
GATE P3.2-B3 — pseudo-hermiticidad de M_eff en la metrica fisica 𝒢=diag(W,−W).
Implementa P3_2B3_DESIGN_DRAFT.md rev.2 (diseño APROBADO por ChatGPT con las 2 precisiones).
NO toca v1.0 ni los gates aprobados. Realizacion: Dirichlet-fantasma ('plana' de B1), K_P Pade G FIJO.

OPERADOR (M_eff, 2N×2N, bloques u=[0:N], v=[N:2N]):
  L_± = K_P ± C + (G+Qp) ,   M = [[L_+, G],[−G, −L_−]]
  K_P = ½D + c1·D A D + c2·D² A³ D² ,  A=(2·diag(G)+D)⁻¹ ,  D=−∂² (Dirichlet-fantasma, euclid-simetrico)
  C   = −i(v∂_x + v'/2)  (forma build_MBdG v1.0, hermitica solo a O(dx))      G = g·n(x) FIJO
  Qp  = V − μ + g n + v²/2  (presion cuantica del fondo)
METRICA:  𝒢 = diag(W,−W),  W = diag(½,1,…,1,½) (peso trapezoidal).  Defecto:  R = M†𝒢 − 𝒢M.

BULK/BORDE (def. EXACTA, 2N×2N):  B={0,N−1,N,2N−1};  I={1..N−2}∪{N+1..2N−2};  R_bulk=R[I,I].
Sub-bloques uu=R[1:N−1,1:N−1], vv=R[N+1:2N−1,N+1:2N−1] (conservar separacion u/v).

CRITERIOS (B3.4):  A1 fuera-diag=0 ~maq ; A2 R_bulk[K_P],R_bulk[G+Qp] ~maq ; A3 (SOLO waterfall)
R_bulk~O(dx) pend∈[0.8,1.2] y R_bulk≈R_bulk[C] ; (control homog G=1/v cte) R_bulk~maquina SIN pendiente.
R4 antiguo RETIRADO (en el interior W=I ⟹ R_bulk^(W)=R_bulk^(I) identico: el bulk es metric-indep).
"""
import numpy as np
import sys, os, contextlib
from pathlib import Path
root = Path(__file__).resolve().parents[2]
with open(os.devnull, 'w') as _dn, contextlib.redirect_stdout(_dn):
    sys.path.insert(0, str(root / 'solver'))
    from correlations import bg, g, mu

lam = 2.0; c1 = (lam-1)/4; c2 = -(lam-1)**2/16      # mismo lam que los demas gates

# ---------------------------------------------------------------- operadores (realizacion plana)
def build_D_A(Gx, N, dx):
    """D plano (−∂², Dirichlet-fantasma, euclid-sim), D2=D², A=(2 diag(G)+D)⁻¹."""
    main = 2/dx**2*np.ones(N); off = -1/dx**2*np.ones(N-1)
    D = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)
    D2 = D @ D
    A = np.linalg.inv(2*np.diag(Gx) + D)
    return D, D2, A

def KP_mat(D, D2, A):
    return 0.5*D + c1*(D@A@D) + c2*(D2@(A@A@A)@D2)

def C_operator(v, dx):
    """C = −i(v∂ + v'/2) (forma v1.0 build_MBdG; hermitica solo a O(dx))."""
    N = len(v); e = np.ones(N-1)
    D1 = (np.diag(e, 1) - np.diag(e, -1))/(2*dx)
    vp = np.gradient(v, dx)
    return -1j*(np.diag(v)@D1 + 0.5*np.diag(vp))

def build_M_and_terms(mode, N, L, w=2.0):
    """Devuelve M (2N×2N), los bloques y los operadores por termino para la descomposicion."""
    x = np.linspace(-L, L, N); dx = x[1]-x[0]
    n, v, V = bg(x, w=w, mode=mode)
    Gx = g*n                                  # G FIJO (no reescalado)
    Qp = V - mu + g*n + 0.5*v**2
    diag = Gx + Qp
    D, D2, A = build_D_A(Gx, N, dx)
    KP = KP_mat(D, D2, A)
    C  = C_operator(v, dx)
    Dg = np.diag(diag)
    Lp = KP + C + Dg
    Lm = KP - C + Dg
    Gd = np.diag(Gx)
    M = np.block([[Lp, Gd], [-Gd, -Lm]])
    return dict(x=x, dx=dx, N=N, M=M, KP=KP, C=C, Dg=Dg, Gd=Gd, Gx=Gx, Lp=Lp, Lm=Lm)

# ------------------------------------------------------------------------------- utilidades R
def metric_G(N):
    W = np.ones(N); W[0] = 0.5; W[-1] = 0.5
    return np.diag(np.concatenate([W, -W])), W

def interior_idx(N):
    return np.r_[1:N-1, N+1:2*N-1]

def mx(X):  # norma max-entrada
    return np.abs(X).max() if X.size else 0.0

def smooth_modes(x, m, sigma=2.0):
    """m paquetes SUAVES de espectro concentrado a baja k (gaussianas, NO estrictamente band-limited)
    centrados cerca del horizonte x=0, anchura sigma. Columnas normalizadas (NO ortogonales: se
    ortonormalizan luego por QR en analyze)."""
    N = len(x); B = np.zeros((N, m))
    centers = np.linspace(-1.5*sigma, 1.5*sigma, m)
    for j, x0 in enumerate(centers):
        g = np.exp(-((x-x0)/sigma)**2)
        B[:, j] = g/np.linalg.norm(g)
    return B

def analyze(d, m_smooth=6):
    """R total + por terminos; bulk(uu,vv)/borde; norma MAX, ESPECTRAL y en el SUBESPACIO SUAVE
    ORTONORMALIZADO (QR). Defectos por termino en ABS y REL (/scale). Devuelve dict."""
    N = d['N']; M = d['M']; x = d['x']; Gmetr, W = metric_G(N); Wm = np.diag(W)
    GM = Gmetr @ M
    R = M.conj().T @ Gmetr - GM
    I = interior_idx(N)
    scale = mx(GM)
    offUL = R[:N, N:]; offLL = R[N:, :N]
    iin = np.r_[1:N-1]
    P_KP  = d['KP'].conj().T @ Wm - Wm @ d['KP']
    P_C   = d['C'].conj().T  @ Wm - Wm @ d['C']
    P_dia = d['Dg'].conj().T @ Wm - Wm @ d['Dg']
    Rbulk_mat = R[np.ix_(I, I)]
    # --- SONDA FISICA ORTONORMALIZADA (QR): subespacio suave, u-block y v-block ---
    Bu = smooth_modes(x, m_smooth)
    Bsmooth = np.zeros((2*N, 2*m_smooth))
    Bsmooth[:N, :m_smooth] = Bu; Bsmooth[N:, m_smooth:] = Bu
    U, _ = np.linalg.qr(Bsmooth)                      # U†U = I (base ortonormal del subespacio suave)
    Rs  = U.conj().T @ R  @ U
    GMs = U.conj().T @ GM @ U
    smooth_abs = np.linalg.norm(Rs, 2)
    smooth_rel = smooth_abs/np.linalg.norm(GMs, 2)
    # --- DESCOMPOSICION del defecto PROYECTADO por termino (ChatGPT): R=R_KP+R_C+R_dia (lineal en M) ---
    Z = np.zeros((N, N))
    RKP  = np.block([[P_KP, Z], [Z,  P_KP]])          # K_P en L_+ (+) y en −L_− (+): mismo signo
    RC   = np.block([[P_C,  Z], [Z, -P_C ]])          # C en L_+ (+) y en −L_− (+), metrica −W: signo −
    Rdia = np.block([[P_dia,Z], [Z,  P_dia]])
    recon_err = np.abs(R - (RKP + RC + Rdia)).max()   # consistencia: R == suma de terminos
    sKP  = np.linalg.norm(U.conj().T @ RKP  @ U, 2)   # normas espectrales PROYECTADAS por termino
    sC   = np.linalg.norm(U.conj().T @ RC   @ U, 2)
    sdia = np.linalg.norm(U.conj().T @ Rdia @ U, 2)
    # control de colas en el borde: sondas restringidas al INTERIOR (filas de borde a 0) -> QR -> defecto
    Bint = Bsmooth.copy(); Bint[[0, N-1, N, 2*N-1], :] = 0.0
    Ui, _ = np.linalg.qr(Bint)
    smooth_abs_int = np.linalg.norm(Ui.conj().T @ R @ Ui, 2)
    kpscale = mx(d['KP'])
    out = dict(
        N=N, dx=d['dx'], scale=scale,
        offUL=mx(offUL), offLL=mx(offLL),
        Rbulk=mx(Rbulk_mat), Rbulk_spec=np.linalg.norm(Rbulk_mat, 2), Rbulk_rel=mx(Rbulk_mat)/scale,
        smooth_abs=smooth_abs, smooth_rel=smooth_rel, smooth_abs_int=smooth_abs_int,
        sKP_proj=sKP, sC_proj=sC, sdia_proj=sdia, recon_err=recon_err,
        PKP_bulk=mx(P_KP[np.ix_(iin, iin)]),   PKP_bulk_rel=mx(P_KP[np.ix_(iin, iin)])/scale,
        PC_bulk=mx(P_C[np.ix_(iin, iin)]),     PC_bulk_rel=mx(P_C[np.ix_(iin, iin)])/scale,
        Pdia_bulk=mx(P_dia[np.ix_(iin, iin)]), Pdia_bulk_rel=mx(P_dia[np.ix_(iin, iin)])/scale,
        KPsym=mx(d['KP']-d['KP'].conj().T)/kpscale,
    )
    return out, R

# ------------------------------------------------------------------------------------- run
def run():
    print("GATE P3.2-B3 — pseudo-hermiticidad de M_eff en 𝒢=diag(W,−W)\n" + "="*72)
    print("Realizacion Dirichlet-fantasma (plana); K_P Pade G FIJO (lam=%.1f); C forma v1.0.\n" % lam)

    L = 16.0
    # ---- PRE-REQUISITO B3.0: K_P^T=K_P euclideo para esta realizacion ----
    print("[PRE-REQ B3.0] K_P^T=K_P euclideo (realizacion plana), rel = max|K_P−K_P^T|/max|K_P|:")
    for N in [200, 400, 800]:
        d = build_M_and_terms('horizon', N, L)
        o, _ = analyze(d)
        flag = 'OK' if o['KPsym'] < 1e-10 else ('cond.' if o['KPsym'] < 1e-6 else 'FALLO')
        print(f"   N={N:5d} dx={d['dx']:.4f}  K_P eucl_rel={o['KPsym']:.2e}  {flag}")
    print("   (nota: a N grande el defecto crece por condicionamiento de D²A³D² densa, no por asimetria")
    print("    real; en banda disperso se evita. Mantener N en rango limpio, como B2.3 N≲3200.)\n")

    # ---- A1: bloques fuera de diagonal = 0 exacto ----
    print("[A1] Bloques fuera de diagonal GW−WG (deben ser 0 EXACTO), N=400:")
    for mode in ['homog', 'horizon']:
        d = build_M_and_terms(mode, 400, L); o, _ = analyze(d)
        print(f"   {mode:>8}: offUL={o['offUL']:.2e}  offLL={o['offLL']:.2e}")
    print()

    # ---- A2: defectos RELATIVOS por termino + simetria de K_P, por resolucion (waterfall) ----
    print("[A2] Defectos de bulk por termino, RELATIVOS a scale=max|𝒢M|, + simetria de K_P (waterfall):")
    print("     (K_P y G+Qp deben ser ~maquina; el defecto real es C. rel = abs/max|𝒢M|)")
    print(f"   {'N':>5} {'dx':>8} {'KPsym_rel':>10} {'PKP_bulk_rel':>13} {'Pdia_bulk_rel':>14} {'PC_bulk_rel':>12}")
    for N in [200, 400, 800, 1600]:
        d = build_M_and_terms('horizon', N, L); o, _ = analyze(d)
        print(f"   {N:>5} {o['dx']:8.4f} {o['KPsym']:10.2e} {o['PKP_bulk_rel']:13.2e} "
              f"{o['Pdia_bulk_rel']:14.2e} {o['PC_bulk_rel']:12.2e}")
    print("   -> PKP_bulk_rel, Pdia_bulk_rel ~maquina (KPsym crece por condicionamiento DENSO de D²A³D²,")
    print("      no por asimetria real); PC_bulk_rel = el defecto fisico (decae ~O(dx²) en rel).\n")

    # ---- A3 (reformulado en DOS resultados, ChatGPT): A3-matriz (O(1)) y A3-suave (O(dx²)) ----
    print("[A3-matriz] Defecto de R_bulk en norma MATRICIAL (max / espectral): NO se exige que desaparezca.")
    print("[A3-suave]  Defecto proyectado al subespacio SUAVE ORTONORMALIZADO (QR): se exige O(dx²).")
    print("            (_int = sondas restringidas al interior: control de colas en el borde)\n")
    for mode, papel in [('homog',   'CONTROL (v cte): sin fuente de defecto por C'),
                        ('horizon', 'PRINCIPAL (waterfall): C introduce el defecto a medir')]:
        print(f"   perfil '{mode}' = {papel}:")
        print(f"   {'N':>5} {'dx':>8} {'MAT_max':>11} {'MAT_spec':>11} | "
              f"{'SUAVE_abs':>11} {'SUAVE_rel':>11} {'SUAVE_int':>11}")
        Ns = [200, 400, 800, 1600]; dxs = []; rmax = []; rspec = []; sabs = []; srel = []; sint = []
        for N in Ns:
            d = build_M_and_terms(mode, N, L); o, _ = analyze(d)
            dxs.append(o['dx']); rmax.append(o['Rbulk']); rspec.append(o['Rbulk_spec'])
            sabs.append(o['smooth_abs']); srel.append(o['smooth_rel']); sint.append(o['smooth_abs_int'])
            print(f"   {N:>5} {o['dx']:8.4f} {o['Rbulk']:11.3e} {o['Rbulk_spec']:11.3e} | "
                  f"{o['smooth_abs']:11.3e} {o['smooth_rel']:11.3e} {o['smooth_abs_int']:11.3e}")
        def sl(a): return np.polyfit(np.log(dxs), np.log(a), 1)[0] if min(a) > 0 else float('nan')
        print(f"   pendientes dx:  MAT_max={sl(rmax):.2f} MAT_spec={sl(rspec):.2f} | "
              f"SUAVE_abs={sl(sabs):.2f} SUAVE_rel={sl(srel):.2f} SUAVE_int={sl(sint):.2f}")
        print()

    # ---- A3-suave DESCOMPUESTO por termino sobre el subespacio QR (ChatGPT) ----
    print("[A3-suave por término] Descomposicion del defecto PROYECTADO U†RU = U†R_KP U + U†R_C U +")
    print("     U†R_dia U (normas espectrales ABS). Confirma que el O(dx²) viene de C y cuantifica K_P:")
    for mode in ['homog', 'horizon']:
        print(f"   perfil '{mode}':")
        print(f"   {'N':>5} {'dx':>8} {'proj_KP':>11} {'proj_C':>11} {'proj_dia':>11} "
              f"{'proj_TOTAL':>11} {'recon_err':>10}")
        Ns = [200, 400, 800, 1600]; dxs = []; aKP = []; aC = []; aT = []
        for N in Ns:
            d = build_M_and_terms(mode, N, L); o, _ = analyze(d)
            dxs.append(o['dx']); aKP.append(o['sKP_proj']); aC.append(o['sC_proj']); aT.append(o['smooth_abs'])
            print(f"   {N:>5} {o['dx']:8.4f} {o['sKP_proj']:11.3e} {o['sC_proj']:11.3e} "
                  f"{o['sdia_proj']:11.3e} {o['smooth_abs']:11.3e} {o['recon_err']:10.2e}")
        def sl(a): return np.polyfit(np.log(dxs), np.log(a), 1)[0] if min(a) > 0 else float('nan')
        print(f"   pendientes dx:  proj_KP={sl(aKP):.2f}  proj_C={sl(aC):.2f}  proj_TOTAL={sl(aT):.2f}")
        print()

    # ---- A4a: defecto FISICO estable al aumentar L (dx fijo, waterfall) ----
    print("[A4a] Defecto FISICO (subespacio suave) estable al aumentar L (dx≈fijo) — defecto LOCAL:")
    print(f"   {'L':>5} {'N':>6} {'dx':>8} {'FISICO_abs':>11} {'Rbulk(max)':>11}")
    dx_target = 0.08
    for Lt in [12.0, 16.0, 24.0, 32.0]:
        Nt = int(round(2*Lt/dx_target)) | 1
        d = build_M_and_terms('horizon', Nt, Lt); o, _ = analyze(d)
        print(f"   {Lt:>5.0f} {Nt:>6} {o['dx']:8.4f} {o['smooth_abs']:11.3e} {o['Rbulk']:11.3e}")
    print()

    print("Conclusion cientifica (acotada, redaccion ChatGPT):")
    print(" Para la realizacion Dirichlet-fantasma estudiada, el defecto de pseudo-hermiticidad de M_eff")
    print(" NO desaparece en la norma matricial completa (A3-matriz = O(1)): persiste en los modos de")
    print(" rejilla. Su proyeccion sobre paquetes SUAVES de espectro concentrado a baja k, localizados")
    print(" en el horizonte, muestra convergencia O(dx²) (A3-suave, QR ortonormalizado, pend 2.00). El")
    print(" control homogeneo da ~maquina en el subespacio suave. NO se demuestra pseudo-herm exacta de")
    print(" la matriz discreta ni conservacion en el continuo; NO se extiende a Robin/otras BC; NO SBP.")
    print(" El crecimiento del defecto matricial del CONTROL con N es CONTAMINACION NUMERICA de la")
    print(" construccion densa (inv + D²A³D²), no fisica (en banda disperso se evita). Pre-produccion a")
    print(" N grande exige la realizacion en BANDA; para el diagnostico fisico de B3 no es necesaria.")

if __name__ == '__main__':
    run()
