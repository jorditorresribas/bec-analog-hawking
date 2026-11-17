"""
BLOQUE 2 — Gate 1 ALGEBRAICO (N pequeño, SIN observables Hawking).

Verifica el ensamblado del sistema aumentado del BVP global de scattering descrito en
STABILIZED_SOLVER_BLOCK1.md (esquema de caja + matching por amplitudes aumentadas + 4 restricciones).

NO calcula espectro de Hawking, NO n_L, NO afirma pseudo-unitaridad. Solo estructura algebraica:
  (1) construccion y ordenacion de la base modal Xi_L, Xi_R por el clasificador (v_g, Omega, Im k);
  (2) invariante de conteo  #in + #grow-evan = 4;
  (3) cuadratura del sistema  (4N+12) x (4N+12);
  (4) sigma_min(Xi_b)  -> independencia modal;
  (5) rango numerico de M_aug  == 4N+12  (ker trivial => unicidad);  reportar sigma_min/sigma_max (escala);
  (6) consistencia del matching: resolver los 3 problemas (s = e1,e2,e3), residuo y verificacion de que
      las filas de matching y las restricciones se satisfacen en la solucion.

Baseline LOCAL lambda=1: usa el operador BdG de v1.0 (A_at) sobre el perfil CONGELADO (solo lectura de
out/waterfall_state.npz). NO importa ni modifica bdg_smatrix.py; reimplementa A localmente para quedar
autocontenido y no tocar v1.0.
"""
import numpy as np
import scipy.linalg as la
from scipy.interpolate import interp1d

g = 1.0
STATE = str(__import__("pathlib").Path(__file__).resolve().parents[2] / "out" / "waterfall_state.npz")
d = np.load(STATE)
xg, n0g, v0g, Vg = d["x"], d["n"], d["v"], d["V"]
mu = float(d["mu"])
dvg = np.gradient(v0g, xg)
fn = interp1d(xg, n0g, kind="cubic"); fv = interp1d(xg, v0g, kind="cubic")
fV = interp1d(xg, Vg, kind="cubic"); fdv = interp1d(xg, dvg, kind="cubic")

nL, vL, VL = n0g[:50].mean(), v0g[:50].mean(), Vg[:50].mean(); cL = np.sqrt(g*nL)
nR, vR, VR = n0g[-50:].mean(), v0g[-50:].mean(), Vg[-50:].mean(); cR = np.sqrt(g*nR)


def A_at(x, w):
    n, v, V, dv = float(fn(x)), float(fv(x)), float(fV(x)), float(fdv(x))
    a = -1j*dv + v**2 + 4*g*n + 2*V - 2*mu - 2*w
    b = 1j*dv + v**2 + 4*g*n + 2*V - 2*mu + 2*w
    return np.array([[0,1,0,0],[a,-2j*v,2*g*n,0],[0,0,0,1],[2*g*n,0,b,2j*v]], dtype=complex)


def A_homog(n, v, V, w):
    a = v**2 + 4*g*n + 2*V - 2*mu - 2*w
    b = v**2 + 4*g*n + 2*V - 2*mu + 2*w
    return np.array([[0,1,0,0],[a,-2j*v,2*g*n,0],[0,0,0,1],[2*g*n,0,b,2j*v]], dtype=complex)


def norm_current(Y, Vf):
    u, up, v, vp = Y
    return (np.imag(np.conj(u)*up) + Vf*abs(u)**2 + np.imag(np.conj(v)*vp) - Vf*abs(v)**2)


def modes(n, v, V, c, w, imtol=1e-8):
    """4 modos clasificados. Propagantes normalizados a |J|=1; evanescentes a norma unidad."""
    lam, vec = la.eig(A_homog(n, v, V, w))
    out = []
    for i in range(4):
        k = -1j*lam[i]; Y = vec[:, i].copy()
        prop = abs(k.imag) < imtol
        if prop:
            kr = k.real; Om = w - v*kr; vg = v + (c**2*kr + kr**3/2)/Om
            J = norm_current(Y, v); Y = Y/np.sqrt(abs(J))
            out.append(dict(k=kr, prop=True, Om=Om, vg=vg, norm=int(np.sign(Om)),
                            dir=int(np.sign(vg)), Y=Y))
        else:
            # guardamos Im(k) crudo; el crecimiento hacia el EXTERIOR depende del extremo
            # (|e^{ikx}|=e^{-Im(k)x}): izq x->-inf crece si Im k>0; der x->+inf crece si Im k<0.
            out.append(dict(k=k, prop=False, imk=float(k.imag), Y=Y/la.norm(Y)))
    return out


def build_Xi(side, w):
    """Devuelve Xi (4x4), etiquetas de columna y el nº de 'fijadas' (in + grow) de este extremo.

    Clasificacion de evanescentes DEPENDIENTE DEL EXTREMO (crece hacia el exterior -> excluir):
      izq (x->-inf): crece si Im k>0 ; decae si Im k<0.
      der (x->+inf): crece si Im k<0 ; decae si Im k>0.
    Orden DETERMINISTA dentro de cada grupo: por (Re k, Im k) ascendente, para que la columna
    de cada canal represente SIEMPRE el mismo modo fisico al comparar S(w) entre frecuencias.
    """
    kkey = lambda m: (float(np.real(m["k"])), float(np.imag(m["k"])))
    if side == "L":
        ms = modes(nL, vL, VL, cL, w)
        prop = [m for m in ms if m["prop"]]; evan = [m for m in ms if not m["prop"]]
        ins  = [m for m in prop if m["dir"] > 0]      # izq: v_g>0 = IN (hacia horizonte)
        outs = [m for m in prop if m["dir"] < 0]
        dec  = [m for m in evan if m["imk"] < 0]      # Im k<0 decae a x->-inf
        grow = [m for m in evan if m["imk"] > 0]      # Im k>0 crece a x->-inf -> excluir
    else:
        ms = modes(nR, vR, VR, cR, w)
        prop = [m for m in ms if m["prop"]]; evan = [m for m in ms if not m["prop"]]
        ins  = [m for m in prop if m["dir"] < 0]      # der: v_g<0 = IN (viene de +inf)
        outs = [m for m in prop if m["dir"] > 0]
        dec  = [m for m in evan if m["imk"] > 0]      # Im k>0 decae a x->+inf
        grow = [m for m in evan if m["imk"] < 0]      # Im k<0 crece a x->+inf -> excluir
    ins, outs = sorted(ins, key=kkey), sorted(outs, key=kkey)
    dec, grow = sorted(dec, key=kkey), sorted(grow, key=kkey)
    cols = ins + outs + dec + grow
    labels = (["in"]*len(ins) + ["out"]*len(outs) + ["dec"]*len(dec) + ["grow"]*len(grow))
    n_fixed = len(ins) + len(grow)
    Xi = np.column_stack([m["Y"] for m in cols])
    return Xi, labels, n_fixed, cols


def assemble(w, N, X0):
    """Ensambla M_aug (denso) y devuelve M, indices, Xi_L, Xi_R, labels."""
    xs = np.linspace(-X0, X0, N+1); dx = xs[1]-xs[0]
    nfield = 4*(N+1)
    ndof = nfield + 8                    # + a_L(4) + a_R(4)
    def iy(j): return 4*j                 # inicio de Y_j
    iaL = nfield; iaR = nfield+4
    M = np.zeros((ndof, ndof), dtype=complex)
    row = 0
    # interior: caja  (I - dx/2 A_{j+1}) Y_{j+1} - (I + dx/2 A_j) Y_j = 0
    for j in range(N):
        Aj  = A_at(xs[j],   w); Aj1 = A_at(xs[j+1], w)
        B = np.eye(4) - dx/2*Aj1
        C = np.eye(4) + dx/2*Aj
        M[row:row+4, iy(j+1):iy(j+1)+4] = B
        M[row:row+4, iy(j):iy(j)+4]     = -C
        row += 4
    XiL, labL, nfixL, _ = build_Xi("L", w)
    XiR, labR, nfixR, _ = build_Xi("R", w)
    # matching izq: Y_0 - Xi_L a_L = 0
    M[row:row+4, iy(0):iy(0)+4] = np.eye(4); M[row:row+4, iaL:iaL+4] = -XiL; row += 4
    # matching der: Y_N - Xi_R a_R = 0
    M[row:row+4, iy(N):iy(N)+4] = np.eye(4); M[row:row+4, iaR:iaR+4] = -XiR; row += 4
    # restricciones: construir selectores segun labels (in -> prescrita; grow -> =0)
    def sel(labels, tag):
        return [i for i, L in enumerate(labels) if L == tag]
    # izq: 1 in + 1 grow (en ventana)
    constr = []   # (col_base, idx_local, tag) ; tag='in' se prescribe, 'grow' =0
    for i in sel(labL, "grow"): constr.append((iaL, i, "grow0"))
    for i in sel(labL, "in"):   constr.append((iaL, i, "inL"))
    for i in sel(labR, "grow"): constr.append((iaR, i, "grow0"))
    for i in sel(labR, "in"):   constr.append((iaR, i, "inR"))
    # añadir filas de restriccion
    rhs_tags = []
    for (base, idx, tag) in constr:
        M[row, base+idx] = 1.0; rhs_tags.append((row, tag)); row += 1
    info = dict(ndof=ndof, N=N, X0=X0, dx=dx, iaL=iaL, iaR=iaR, iyN=iy(N),
                XiL=XiL, XiR=XiR, labL=labL, labR=labR, nfixL=nfixL, nfixR=nfixR,
                rhs_tags=rhs_tags, nrows=row, constr=constr)
    return M, info


def rhs_for(info, s):
    """s = (s_L, s_R1, s_R2) amplitudes entrantes prescritas. Devuelve b."""
    b = np.zeros(info["ndof"], dtype=complex)
    # asignar a cada fila de restriccion su valor
    inL_seen = 0; inR_seen = 0
    for (r, tag) in info["rhs_tags"]:
        if tag == "grow0":
            b[r] = 0.0
        elif tag == "inL":
            b[r] = s[0]; inL_seen += 1
        elif tag == "inR":
            b[r] = s[1] if inR_seen == 0 else s[2]; inR_seen += 1
    return b


def run(w=0.10, N=60, X0=12.0, verbose=True):
    M, info = assemble(w, N, X0)
    ndof = info["ndof"]
    ok = True; report = {}
    # (3) cuadratura
    square = (M.shape[0] == M.shape[1] == ndof == 4*N+12) and (info["nrows"] == ndof)
    report["square"] = (M.shape, 4*N+12, info["nrows"], square); ok &= square
    # (2) invariante de conteo
    inv = info["nfixL"] + info["nfixR"]
    report["count_invariant"] = (info["nfixL"], info["nfixR"], inv, inv == 4); ok &= (inv == 4)
    # (2b) GUARD de ventana 3x3: el solve de 3 problemas (e1,e2,e3) SOLO es valido si la
    # ESTRUCTURA COMPLETA de canales es la prevista (no basta el nº de entrantes: otra
    # distribucion prop/evan podria dar los mismos #in). Exigimos:
    #   izq: 1 in, 1 out, 1 dec, 1 grow  |  der: 2 in, 2 out, 0 dec, 0 grow.
    cL_ = info["labL"]; cR_ = info["labR"]
    n_in_L = cL_.count("in"); n_in_R = cR_.count("in")
    window_3x3 = (cL_.count("in") == 1 and cL_.count("out") == 1 and
                  cL_.count("dec") == 1 and cL_.count("grow") == 1 and
                  cR_.count("in") == 2 and cR_.count("out") == 2 and
                  cR_.count("dec") == 0 and cR_.count("grow") == 0)
    report["window_3x3"] = (n_in_L, n_in_R, window_3x3); ok &= window_3x3
    # (4) sigma_min de las bases modales
    sL = la.svdvals(info["XiL"]); sR = la.svdvals(info["XiR"])
    report["sigma_min_Xi"] = (sL.min(), sL.max(), sR.min(), sR.max())
    # (5) rango numerico de M_aug
    sv = la.svdvals(M)
    tol = sv.max() * ndof * np.finfo(float).eps
    rank = int((sv > tol).sum())
    report["rank"] = (rank, ndof, rank == ndof, sv.min(), sv.max(), sv.min()/sv.max())
    ok &= (rank == ndof)
    # (6) consistencia del matching: resolver 3 problemas e1,e2,e3  (SOLO en ventana 3x3)
    if not window_3x3:
        report["solve_residual"] = np.nan; report["matching_err"] = np.nan
        report["grow_zero_err"] = np.nan; report["in_amp_err"] = np.nan
        report["PASS"] = ok   # ok ya es False por el guard (2b)
        if verbose:
            print("="*72)
            print(f"BLOQUE 2 — Gate 1 ALGEBRAICO   (w={w}, N={N}, X0={X0}, dx={info['dx']:.4f})")
            print("="*72)
            print(f"  Clasificacion izq: {info['labL']}  der: {info['labR']}")
            print(f"  [2b] canales IN: izq={n_in_L}, der={n_in_R}  -> ventana 3x3? {window_3x3}")
            print(f"  FUERA DE VENTANA 3x3: no se ejecuta el solve de scattering.  RESULTADO: FAIL(guard)")
        return report
    iLin = info["labL"].index("in")                                   # columna 'in' izquierda
    iRin = [i for i, L in enumerate(info["labR"]) if L == "in"]       # columnas 'in' derechas
    res_list = []; match_err = []; constr_err = []; in_err = []
    for s in (np.array([1,0,0]), np.array([0,1,0]), np.array([0,0,1])):
        b = rhs_for(info, s)
        x = la.solve(M, b)
        res_list.append(la.norm(M@x - b))
        # verificar Y_0 = Xi_L a_L y Y_N = Xi_R a_R
        Y0 = x[0:4]; YN = x[info["iyN"]:info["iyN"]+4]
        aL = x[info["iaL"]:info["iaL"]+4]; aR = x[info["iaR"]:info["iaR"]+4]
        match_err.append(max(la.norm(Y0 - info["XiL"]@aL), la.norm(YN - info["XiR"]@aR)))
        # verificar restricciones de crecimiento nulo (grow=0)
        ce = 0.0
        for (r, tag) in info["rhs_tags"]:
            col = int(np.argmax(np.abs(M[r])))
            if tag == "grow0": ce = max(ce, abs(x[col]))
        constr_err.append(ce)
        # verificacion EXPLICITA de las amplitudes entrantes prescritas (distingue fallo de
        # matching de fallo de asignacion de fuentes)
        ie = abs(aL[iLin] - s[0])
        if len(iRin) >= 1: ie = max(ie, abs(aR[iRin[0]] - s[1]))
        if len(iRin) >= 2: ie = max(ie, abs(aR[iRin[1]] - s[2]))
        in_err.append(ie)
    report["solve_residual"] = max(res_list)
    report["matching_err"] = max(match_err)
    report["grow_zero_err"] = max(constr_err)
    report["in_amp_err"] = max(in_err)
    ok &= (max(res_list) < 1e-8) and (max(match_err) < 1e-8) \
          and (max(constr_err) < 1e-8) and (max(in_err) < 1e-8)
    report["PASS"] = ok

    if verbose:
        print("="*72)
        print(f"BLOQUE 2 — Gate 1 ALGEBRAICO   (w={w}, N={N}, X0={X0}, dx={info['dx']:.4f})")
        print("="*72)
        print(f"  Clasificacion izq: {info['labL']}  | fijadas(in+grow)={info['nfixL']}")
        print(f"  Clasificacion der: {info['labR']}  | fijadas(in+grow)={info['nfixR']}")
        print(f"  [2] invariante #in+#grow = {inv}  (==4 ? {inv==4})")
        wn = report["window_3x3"]
        print(f"  [2b] canales IN: izq={wn[0]}, der={wn[1]}  -> ventana 3x3? {wn[2]} (gatea PASS)")
        sh, exp, nr, sq = report["square"]
        print(f"  [3] cuadratura: M={sh}, esperado {exp}x{exp}, filas usadas={nr}  -> {sq}")
        print(f"  [4] sigma(Xi_b) [INFORMATIVO, no gatea PASS]: "
              f"Xi_L min={sL.min():.3e} max={sL.max():.3e} | "
              f"Xi_R min={sR.min():.3e} max={sR.max():.3e}")
        rk = report["rank"]
        print(f"  [5] rango(M_aug)={rk[0]} / {rk[1]}  ker trivial? {rk[2]}  | "
              f"sigma_min={rk[3]:.3e} sigma_max={rk[4]:.3e} ratio={rk[5]:.3e}")
        print(f"  [6] residuo solve   = {report['solve_residual']:.3e}")
        print(f"      matching err    = {report['matching_err']:.3e}  (Y_b == Xi_b a_b)")
        print(f"      grow==0 err     = {report['grow_zero_err']:.3e}")
        print(f"      in-amp err      = {report['in_amp_err']:.3e}  (a_in == s prescrita)")
        print("-"*72)
        print(f"  RESULTADO: {'PASS' if ok else 'FAIL'}")
    return report


if __name__ == "__main__":
    # test principal
    r = run(w=0.10, N=60, X0=12.0)
    # invariancia del rango/ker con N (no es convergencia fisica; solo estructura)
    print("\n[control] ker trivial e invariante de conteo a varios N (w=0.10, X0=12):")
    for N in (20, 40, 80, 120):
        rr = run(w=0.10, N=N, X0=12.0, verbose=False)
        rk = rr["rank"]
        print(f"   N={N:4d}  rango={rk[0]}/{rk[1]}  ker_trivial={rk[2]}  "
              f"inv={rr['count_invariant'][2]}  ratio_sv={rk[5]:.2e}  PASS={rr['PASS']}")
    # control: otra frecuencia en ventana
    print("\n[control] otra frecuencia en ventana (w=0.05, N=60):")
    run(w=0.05, N=60, X0=12.0)

    # control: CLASIFICACION con evanescentes DERECHOS (frecuencias donde aparecen) ->
    # verifica que la correccion del extremo no es cosmetica. Se muestra Im(k) y la etiqueta:
    # en el extremo derecho el evanescente a EXCLUIR es Im(k)<0 (crece a +inf).
    print("\n[control] clasificacion de evanescentes DERECHOS (barrido de w):")
    print("   (der: grow=excluir <=> Im k<0 ; dec=mantener <=> Im k>0)")
    for w in (0.30, 0.50, 0.80, 1.20):
        XiL, labL, nfixL, _    = build_Xi("L", w)
        XiR, labR, nfixR, colsR = build_Xi("R", w)
        evanR = [(f"Im k={m['imk']:+.3f}", lab)
                 for m, lab in zip(colsR, labR) if not m["prop"]]
        print(f"   w={w:.2f}  izq={labL}(fix={nfixL})  der={labR}(fix={nfixR})  "
              f"inv={nfixL+nfixR}  evan_der={evanR}")

    # control: GUARD de ventana 3x3 -> run() fuera de ventana debe dar FAIL(guard), no PASS
    print("\n[control] guard de ventana 3x3 (w=0.80, fuera de ventana):")
    rg = run(w=0.80, N=60, X0=12.0)
    print(f"   -> PASS={rg['PASS']} (debe ser False; solo 1+1 canales IN, S no es 3x3)")
