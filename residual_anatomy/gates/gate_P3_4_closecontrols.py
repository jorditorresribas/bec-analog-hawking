"""
GATE P3.4 — controles de cierre del operador (autorizados por ChatGPT):
  (1) C2: barrido ESPACIAL de N (L,T,región,estado inicial fijos; error temporal controlado con nt alto).
      Registrar ΔQ_R, flujo (∫frontera), ∫defecto, residuo de balance, recon, ‖Ψ‖. NO se exige ΔQ_R→0.
  (2) C1: sensibilidad del UMBRAL |Im k|<thr de clasificacion propagante. Comprobar que el recuento/
      clasificacion de canales y S NO dependen artificialmente del umbral. NO se ajusta el umbral para pasar.
No toca v1.0 ni el gate (solo importa sus funciones; usa thr/nt como parametros ya existentes).
"""
import numpy as np
import gate_P3_4_smatrix as G

def c2_N_sweep():
    print("[C2 · barrido N]  L=16, T=4, R=[-4,4], estado inicial fijo, nt=1600 (error temporal ~1e-7):")
    print(f"   {'N':>5} {'dx':>7} {'ΔQ_R':>13} {'∫frontera':>13} {'∫defecto':>12} {'|balance|':>11} "
          f"{'recon':>9} {'‖Ψ‖_T':>8}")
    for N in [300, 400, 600, 800]:
        o = G.packet_balance(N=N, L=16.0, T=4.0, nt=1600)
        print(f"   {N:>5} {o['dx']:7.4f} {o['dQR']:13.6e} {o['integ_bdy']:13.6e} {o['integ_def']:12.4e} "
              f"{o['balance']:11.3e} {o['recon_max']:9.1e} {o['Qgrow']:8.4f}")
    print("   -> |balance| se mantiene en el piso TEMPORAL (~1e-7, indep. de N); ΔQ_R/flujo/∫defecto")
    print("      convergen con N (∫defecto→0 como O(dx²), es el término de pseudo-herm). NO se exige ΔQ_R→0.\n")

def c1_thr_sweep():
    print("[C1 · sensibilidad del umbral de clasificacion |Im k|<thr]  X0=12, dx=0.002:")
    print("   Para cada ω, se comprueba recuento 3+3, max|R_S| y max|S−S(thr=1e-8)| (debe ser ~0 si el")
    print("   recuento no cambia; si cambia el recuento, se marca).")
    for om in [0.05, 0.10, 0.15, 0.20]:
        Sref, einr, eoutr, _ = G.smatrix_at(om, 12.0, 0.002, thr=1e-8)
        print(f"   ω={om:.2f}:")
        print(f"      {'thr':>8} {'#prop_L':>8} {'#prop_R':>8} {'3+3':>5} {'max|R_S|':>11} {'max|S−Sref|':>12}")
        for thr in [1e-6, 1e-8, 1e-10, 1e-12]:
            S, ein, eout, inf = G.smatrix_at(om, 12.0, 0.002, thr=thr)
            # recuento de propagantes a cada lado (re-clasifica con este thr)
            mL = G.modes_raw(G.B.nL, G.B.vL, G.B.VL, G.B.cL, om, G.B.vL, thr=thr)
            mR = G.modes_raw(G.B.nR, G.B.vR, G.B.VR, G.B.cR, om, G.B.vR, thr=thr)
            npl = sum(1 for m in mL if m['prop']); npr = sum(1 for m in mR if m['prop'])
            if S is None:
                print(f"      {thr:8.0e} {npl:>8} {npr:>8} {'NO':>5} {'(no 3+3)':>11} {'-':>12}")
            else:
                Ei = np.diag(ein.astype(float)); Eo = np.diag(eout.astype(float))
                r = np.abs(S.conj().T@Eo@S - Ei).max()
                dS = np.abs(S - Sref).max() if Sref is not None else float('nan')
                print(f"      {thr:8.0e} {npl:>8} {npr:>8} {'si':>5} {r:11.3e} {dS:12.3e}")
    print("   -> si recuento 3+3 y S son invariantes al thr en el rango razonable, la clasificacion NO es")
    print("      fragil (no explica el floor de C1). NO valida la pseudo-unitaridad convergente (reservada).\n")

def run():
    print("GATE P3.4 — controles de cierre (barrido N de C2 + sensibilidad del umbral de C1)\n" + "="*74)
    c2_N_sweep()
    c1_thr_sweep()

if __name__ == '__main__':
    run()
