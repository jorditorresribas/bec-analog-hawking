"""
BLOQUE 2 — B2.4b: REANALISIS de FASE y SENSIBILIDAD de S(X0).

Objetivo (aprobado por ChatGPT con el signo corregido): comprobar si la variacion de S con X0
observada en B2.4 es una ROTACION DE FASE por la referencia espacial de las amplitudes de borde,
NO una no-convergencia fisica; y separar la sensibilidad fisica real de la inestabilidad del
calculo adjunto.

Tres bloques:
  A. Convergencia de MODULOS |S_ext,ij(X0)| (invariantes de fase).
  B. Reconstruccion S_fis = D_out^{-1} S_ext D_in  (fases EXACTAS de k, bordes y convencion de modes;
     NO se ajustan para minimizar). Control algebraico: D_out S_fis D_in^{-1} == S_ext.
     Convergencia de S_fis (compleja) entre X0 (Frobenius + elemento a elemento).
  C. Perturbacion controlada: delta s ~ y^dag (delta b - delta M x). Se compara la prediccion adjunta
     con la diferencia finita (resolviendo el sistema perturbado), para varias magnitudes.
     Si delta_adj es grande, esta prueba revela si el adjunto es fiable o no.

Convencion de fase (ChatGPT): S_fis[o,p] = exp(-i k_out(o) x_b,out(o)) S_ext[o,p] exp(+i k_in(p) x_b,in(p)),
canales izquierdos en x_b=-X0, derechos en x_b=+X0.

NO modifica B2.1/B2.2/B2.3/v1.0. NO toca BLOCK2_LOG.md.
"""
import numpy as np
import scipy.linalg as la
import gate_B2_1_algebraic as g1
import gate_B2_2_conditioning as g2

RNG = np.random.default_rng(1)


def channel_ks(w):
    """k reales y signo del borde (x_b=+-X0) de los 3 canales OUT y los 3 IN, en el orden que usan
       out_indices (out: left,right1,right2) y rhs_for (in: left, right1, right2)."""
    _, labL, _, colsL = g1.build_Xi("L", w)
    _, labR, _, colsR = g1.build_Xi("R", w)
    oR = [i for i, L in enumerate(labR) if L == "out"]
    iR = [i for i, L in enumerate(labR) if L == "in"]
    k_out = np.array([colsL[labL.index("out")]["k"], colsR[oR[0]]["k"], colsR[oR[1]]["k"]], float)
    k_in  = np.array([colsL[labL.index("in")]["k"],  colsR[iR[0]]["k"], colsR[iR[1]]["k"]], float)
    xb_out = np.array([-1.0, +1.0, +1.0])   # left out en -X0; right out en +X0
    xb_in  = np.array([-1.0, +1.0, +1.0])
    return k_out, xb_out, k_in, xb_in


def S_ext_and_fis(w, N, X0):
    """S extraida (del solver) y S_fis reconstruida con las fases EXACTAS. Devuelve ambas + control."""
    M, info = g1.assemble(w, N, X0)
    dr, dc = g2.physical_scaling(info)
    S_ext, _, _, _ = g2.smatrix_and_solves(info, M, dr, dc)
    k_out, xb_out, k_in, xb_in = channel_ks(w)
    Dout = np.exp(1j * k_out * xb_out * X0)           # diagonales (vectores)
    Din  = np.exp(1j * k_in  * xb_in  * X0)
    # S_fis[o,p] = exp(-i k_out x_b,out X0) S_ext[o,p] exp(+i k_in x_b,in X0)
    S_fis = (1.0/Dout)[:, None] * S_ext * Din[None, :]
    # control algebraico: D_out S_fis D_in^{-1} debe reproducir S_ext
    S_back = Dout[:, None] * S_fis * (1.0/Din)[None, :]
    recon_err = float(np.max(np.abs(S_back - S_ext)))
    return S_ext, S_fis, recon_err


# ------------------- A + B: convergencia modulos y S_fis ----------------------
def convergence(w, X0_list, dx=0.02, ref_X0=16):
    print("="*96)
    print(f"[B2.4b A+B] w={w} dx={dx}   modulos |S_ext| y reconstruccion S_fis (fase exacta)")
    print("="*96)
    data = {}
    for X0 in X0_list:
        N = int(round(2*X0/dx))
        S_ext, S_fis, recon = S_ext_and_fis(w, N, X0)
        data[X0] = dict(S_ext=S_ext, S_fis=S_fis, recon=recon)
        print(f"  X0={X0:4.1f} N={N:5d}  control inversion D_out S_fis D_in^-1 == S_ext: max|dif|={recon:.1e}")
    ref = data[ref_X0]
    print(f"\n  Convergencia vs X0={ref_X0}:  (A) |S_ext| moduli   (B) S_fis compleja")
    print(f"    {'X0':>5} {'d|S|_Frob':>11} {'d|S|_max':>10} | {'dS_fis_Frob':>12} {'dS_fis_max':>11}")
    for X0 in X0_list:
        dmodF = la.norm(np.abs(data[X0]["S_ext"]) - np.abs(ref["S_ext"]))
        dmodM = float(np.max(np.abs(np.abs(data[X0]["S_ext"]) - np.abs(ref["S_ext"]))))
        dfisF = la.norm(data[X0]["S_fis"] - ref["S_fis"])
        dfisM = float(np.max(np.abs(data[X0]["S_fis"] - ref["S_fis"])))
        print(f"    {X0:5.1f} {dmodF:11.3e} {dmodM:10.3e} | {dfisF:12.3e} {dfisM:11.3e}")
    print("  (A moduli pequenos => variacion es de fase. B S_fis pequenos => convergencia fisica en la ventana.)")
    return data


# ------------------- C: perturbacion controlada vs adjunto --------------------
def perturbation_test(w, X0, dx=0.02, out_sel=0, in_prob=0, eps_list=(1e-6, 1e-8, 1e-10),
                      do_dM=True):
    """Compara delta_s medido (diferencia finita, sistema perturbado) con prediccion adjunta
       y^dag(delta b - delta M x). out_sel: 0=left-out,1=right-out1,2=right-out2. in_prob: 0/1/2 (e_p+1)."""
    print("\n" + "="*96)
    print(f"[B2.4b C] w={w} X0={X0}  perturbacion controlada (out={out_sel}, in_prob={in_prob})")
    print("="*96)
    N = int(round(2*X0/dx))
    M, info = g1.assemble(w, N, X0)
    dr, dc = g2.physical_scaling(info)
    Mrep = dr[:, None] * M * dc[None, :]
    ndof = info["ndof"]
    sin = np.zeros(3); sin[in_prob] = 1.0
    brep = dr * g1.rhs_for(info, sin)
    z = la.solve(Mrep, brep)
    # selector de salida en vars escaladas: c_rep = conj(dc)*e_idx
    outs = g2.out_indices(info)
    c = np.zeros(ndof, dtype=complex); c[outs[out_sel]] = 1.0
    crep = np.conj(dc) * c
    y = la.solve(Mrep.conj().T, crep)
    s0 = np.vdot(crep, z)
    print(f"  s0 = c^dag z = {s0:.6e}   ||y||={la.norm(y):.3e}  (y^dag brep = {np.vdot(y, brep):.6e})")
    # --- perturbacion de b (delta M = 0) ---
    print("  (1) delta b  (delta M=0):  pred = y^dag db ;  mide = c^dag(z_pert - z)")
    u = RNG.standard_normal(ndof) + 1j*RNG.standard_normal(ndof); u /= la.norm(u)
    for eps in eps_list:
        db = eps * u
        z1 = la.solve(Mrep, brep + db)
        ds_meas = np.vdot(crep, z1) - s0
        ds_pred = np.vdot(y, db)
        rel = abs(ds_meas - ds_pred)/max(abs(ds_meas), 1e-300)
        print(f"      eps={eps:.0e}  ds_meas={ds_meas:+.3e}  ds_pred={ds_pred:+.3e}  err_1er_orden_rel={rel:.3e}")
    # --- perturbacion de M (delta b = 0), rango-1 para no duplicar memoria a N grande ---
    if do_dM:
        print("  (2) delta M (delta b=0, rango-1 eps*||M||*a b^H):  pred = -y^dag (dM x) ;  mide = c^dag(z_pert-z)")
        Mnorm = la.norm(Mrep, 2) if ndof <= 6000 else float(np.max(np.abs(Mrep)))*np.sqrt(ndof)
        a = RNG.standard_normal(ndof) + 1j*RNG.standard_normal(ndof); a /= la.norm(a)
        bb = RNG.standard_normal(ndof) + 1j*RNG.standard_normal(ndof); bb /= la.norm(bb)
        for eps in eps_list:
            scale = eps * Mnorm
            dMz = scale * a * np.vdot(bb, z)             # (dM) z = scale a (b^H z)
            z1 = la.solve(Mrep + scale*np.outer(a, bb.conj()), brep)
            ds_meas = np.vdot(crep, z1) - s0
            ds_pred = -np.vdot(y, dMz)
            rel = abs(ds_meas - ds_pred)/max(abs(ds_meas), 1e-300)
            print(f"      eps={eps:.0e}  ds_meas={ds_meas:+.3e}  ds_pred={ds_pred:+.3e}  err_1er_orden_rel={rel:.3e}")
    print("  (Si el adjunto es fiable, err_1er_orden_rel -> 0 al reducir eps. Si delta_adj esta roto,")
    print("   la prediccion NO coincidira con la diferencia finita aunque ||y|| sea enorme.)")


if __name__ == "__main__":
    # A+B: ventana fisica 8-16 + control de sensibilidad extrema 20-28
    for w in (0.10, 0.05):
        convergence(w, X0_list=(8,10,12,14,16,20,24,28), dx=0.02, ref_X0=16)
    # C: perturbacion en X0 moderado (12, adjunto fiable) y extremo (28, delta_adj roto)
    perturbation_test(0.10, X0=12, dx=0.02, do_dM=True)
    perturbation_test(0.10, X0=28, dx=0.02, do_dM=False)   # dM off a N grande por memoria
    print("\n>>> B2.4b COMPLETO")
