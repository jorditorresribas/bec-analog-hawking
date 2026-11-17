"""
GATE P3.4 (baseline λ=1) — conservación de la matriz S y balance regional.
Implementa P3_4_SMATRIX_DESIGN_DRAFT.md (diseño APROBADO por ChatGPT para el baseline λ=1).
Implementacion NUEVA que USA v1.0 (bdg_smatrix.py, fd_bdg_operator.py) como REFERENCIA; NO los modifica.

Controles:
  C1  pseudo-unitaridad MATRICIAL de S:  max|S†η_out S − η_in|,  ≥4 resoluciones (pendiente en dx) + meseta
      en X0 + recuento de canales estable en la ventana de ω. Registro de canales (k,vg,Ω,sign(Ω),J crudo).
  C2  balance de flujo regional con PAQUETE (independiente), identidad GENERAL con TERMINO DE DEFECTO:
      dQ_R/dt = i⟨Ψ,[M,P_R]Ψ⟩_𝒢 + iΨ†(M†𝒢−𝒢M)P_R Ψ   (el 2º termino NO se asume 0; se mide y reporta).
  A4  ancla: la S de v1.0 reproduce el benchmark (NO sustituye a C1/C2).
Resultados de pseudo-unitaridad global (C1) y balance regional (C2) se mantienen SEPARADOS.
"""
import numpy as np
import sys, os, contextlib
from pathlib import Path
import scipy.linalg as la
root = Path(__file__).resolve().parents[2]
os.chdir(root)                                    # bdg_smatrix usa rutas relativas (out/waterfall_state.npz)
with open(os.devnull, 'w') as _dn, contextlib.redirect_stdout(_dn):
    sys.path.insert(0, str(root / 'solver'))
    import bdg_smatrix as B                       # solver de transferencia λ=1 de v1.0 (REFERENCIA)
    from fd_bdg_operator import build_MBdG, g
    from correlations import bg as corr_bg, mu as corr_mu

# ------------------------------------------------------- clasificacion de canales (J crudo para la tabla)
def modes_raw(n, v, V, c, omega, Vf, thr=1e-8):
    """Como B.modes pero conserva J ANTES de normalizar (Jraw) para el registro de canales.
    thr = umbral de clasificacion propagante |Im k|<thr (parametrizable: control de sensibilidad)."""
    lam, vec = la.eig(B.A_homog(n, v, V, omega))
    out = []
    for i in range(4):
        k = -1j*lam[i]; Y = vec[:, i]; prop = abs(k.imag) < thr
        if prop:
            kr = k.real; Om = omega - v*kr; vg = v + (c**2*kr + kr**3/2)/Om
            Jraw = B.norm_current(Y, Vf)
            out.append(dict(k=kr, prop=True, Om=Om, vg=vg, Jraw=Jraw,
                            norm=int(np.sign(Om)), Y=Y/np.sqrt(abs(Jraw))))
        else:
            out.append(dict(k=k, prop=False, decay=(k.imag < 0), Y=Y))
    return out

def smatrix_at(omega, X0, dx, ev_scale=1.0, thr=1e-8):
    """Reconstruye S 3×3 y (η_in,η_out) a resolucion dx (replica la logica de B.smatrix, param. en dx).
    ev_scale reescala los autovectores EVANESCENTES (control de invariancia: S no debe depender de su
    normalizacion arbitraria). thr = umbral de clasificacion propagante."""
    Prop = B.transfer(omega, X0, dx)
    mL = modes_raw(B.nL, B.vL, B.VL, B.cL, omega, B.vL, thr=thr)
    mR = modes_raw(B.nR, B.vR, B.VR, B.cR, omega, B.vR, thr=thr)
    Lin  = next(m for m in mL if m['prop'] and m['vg'] > 0)
    Lout = next(m for m in mL if m['prop'] and m['vg'] < 0)
    Lev  = [m for m in mL if not m['prop']]
    Lev_d = next(m for m in Lev if m['decay']); Lev_g = next(m for m in Lev if not m['decay'])
    Rins  = [m for m in mR if m['prop'] and m['vg'] < 0]
    Routs = [m for m in mR if m['prop'] and m['vg'] > 0]
    if len(Rins) != 2 or len(Routs) != 2:
        return None, None, None, (mL, mR)              # recuento de canales fuera de la ventana 3+3
    PhiL = np.column_stack([Lin['Y'], Lout['Y'], ev_scale*Lev_d['Y'], ev_scale*Lev_g['Y']])
    PhiR = np.column_stack([Rins[0]['Y'], Rins[1]['Y'], Routs[0]['Y'], Routs[1]['Y']])
    T = la.solve(PhiR, Prop@PhiL)
    S = np.zeros((3, 3), dtype=complex)
    for col, src in enumerate(['Lin', 'Rin1', 'Rin2']):
        aLin = 1.0 if src == 'Lin' else 0.0
        aRin1 = 1.0 if src == 'Rin1' else 0.0
        aRin2 = 1.0 if src == 'Rin2' else 0.0
        Mmat = np.array([[T[0, 1], T[0, 2]], [T[1, 1], T[1, 2]]])
        rhs = np.array([aRin1 - T[0, 0]*aLin, aRin2 - T[1, 0]*aLin])
        aLout, aLev_d = la.solve(Mmat, rhs)
        aL = np.array([aLin, aLout, aLev_d, 0.0])
        S[:, col] = [aLout, T[2, :]@aL, T[3, :]@aL]
    eta_in = np.array([Lin['norm'], Rins[0]['norm'], Rins[1]['norm']])
    eta_out = np.array([Lout['norm'], Routs[0]['norm'], Routs[1]['norm']])
    return S, eta_in, eta_out, (Lin, Lout, Rins, Routs)

def R_S(omega, X0, dx):
    S, ein, eout, _ = smatrix_at(omega, X0, dx)
    if S is None: return None
    Ein = np.diag(ein.astype(float)); Eout = np.diag(eout.astype(float))
    return np.abs(S.conj().T @ Eout @ S - Ein).max()

# -------------------------------------------------------------------- C2: balance regional con paquete
def packet_balance(N=400, L=16.0, T=4.0, nt=400, a=-4.0, b=4.0, x0=-7.0, k0=0.8, sig=2.0):
    """Evoluciona i∂_tΨ=MΨ (M=build_MBdG λ=1) con PROPAGADOR EXACTO expm(−iM dt) (estable; RK4
    explícito es inestable por rigidez |λ|~2/dx²). ALCANCE (ChatGPT): el estado inicial es un ESTADO DE
    PRUEBA LOCALIZADO (gaussiana en la 1ª componente, 2ª=0), NO un paquete construido de un modo BdG
    propagante (no satisface la relación U/V de una rama de dispersión); sirve para VERIFICAR la identidad
    regional, no para afirmar física de Hawking. Verifica ΔQ_R = ∫(flujo frontera + término de defecto) dt
    (trapezoidal). Reporta la magnitud del defecto. (‖Ψ‖ final NO se interpreta como amplificación física.)"""
    x = np.linspace(-L, L, N); dx = x[1]-x[0]
    n, v, V = corr_bg(x, w=2.0, mode='horizon')
    M = build_MBdG(x, n, v, V, corr_mu).toarray()
    W = np.ones(N); W[0] = 0.5; W[-1] = 0.5
    Gm = np.diag(np.concatenate([W, -W]))
    chi = ((x >= a) & (x <= b)).astype(float)
    PR = np.diag(np.concatenate([chi, chi]))
    defectM = M.conj().T @ Gm - Gm @ M                 # M†𝒢 − 𝒢M (NO cero; B3)
    comm_op = Gm @ (M @ PR - PR @ M)                   # 𝒢[M,P_R]
    dir_op  = M.conj().T @ (PR @ Gm) - (PR @ Gm) @ M   # para dQ_R/dt exacto
    dt = T/nt
    U = la.expm(-1j*M*dt)                              # propagador exacto (estable)
    env = np.exp(-((x-x0)/sig)**2) * np.exp(1j*k0*x)   # estado de PRUEBA localizado (NO modo BdG propagante)
    Psi = np.concatenate([env, np.zeros(N)]).astype(complex)
    Psi /= np.sqrt(np.vdot(Psi, Psi).real)
    def QR(P): return np.real(np.vdot(P, (PR @ Gm) @ P))
    def rates(P):
        commut = (1j*np.vdot(P, comm_op @ P)).real
        defect = (1j*np.vdot(P, defectM @ (PR @ P))).real
        direct = (1j*np.vdot(P, dir_op @ P)).real
        return commut, defect, direct
    QR0 = QR(Psi)
    integ_bdy = 0.0; integ_def = 0.0; recon_max = 0.0
    c_prev, d_prev, dr0 = rates(Psi); recon_max = max(recon_max, abs(dr0-(c_prev+d_prev)))
    for _ in range(nt):
        Psi = U @ Psi
        c_cur, d_cur, dr = rates(Psi); recon_max = max(recon_max, abs(dr-(c_cur+d_cur)))
        integ_bdy += 0.5*(c_prev+c_cur)*dt             # trapezoidal
        integ_def += 0.5*(d_prev+d_cur)*dt
        c_prev, d_prev = c_cur, d_cur
    dQR = QR(Psi) - QR0; integ_tot = integ_bdy + integ_def
    return dict(dx=dx, dQR=dQR, integ_tot=integ_tot, integ_bdy=integ_bdy, integ_def=integ_def,
                balance=abs(dQR - integ_tot), recon_max=recon_max,
                rel=abs(dQR - integ_tot)/max(abs(dQR), 1e-30),
                Qgrow=np.sqrt(np.vdot(Psi, Psi).real))   # ||Psi|| final (detecta inestabilidad dinamica)

# ---------------------------------------------------------------------------------------- run
def run():
    print("GATE P3.4 (baseline λ=1) — conservación de S y balance regional\n" + "="*72)
    print("Implementacion NUEVA; usa bdg_smatrix.py (v1.0) como referencia, NO lo modifica.\n")

    # --- Registro de canales (convencion fija) ---
    om0, X0, dx0 = 0.1, 12.0, 0.002
    S, ein, eout, info = smatrix_at(om0, X0, dx0)
    Lin, Lout, Rins, Routs = info
    print(f"[REGISTRO DE CANALES]  ω={om0}  (|J| crudo ANTES de normalizar):")
    print(f"   {'canal':>6} {'lado':>4} {'dir':>4} {'k':>9} {'vg':>8} {'Ω':>9} {'norma':>6} {'Jraw':>11}")
    def row(nm, lado, dirn, m): print(f"   {nm:>6} {lado:>4} {dirn:>4} {m['k']:9.4f} {m['vg']:8.4f} "
                                      f"{m['Om']:9.4f} {m['norm']:>6} {m['Jraw']:11.3e}")
    row('Lin', 'izq', 'IN', Lin); row('Rin1', 'der', 'IN', Rins[0]); row('Rin2', 'der', 'IN', Rins[1])
    row('Lout', 'izq', 'OUT', Lout); row('Rout1', 'der', 'OUT', Routs[0]); row('Rout2', 'der', 'OUT', Routs[1])
    print(f"   η_in={list(ein)}  η_out={list(eout)}  (SEPARADAS)")
    print(f"   |S| =\n{np.abs(S).round(4)}\n")

    # --- C1: pseudo-unitaridad, convergencia en dx (≥4 resoluciones) + pendiente ---
    print("[C1] max|S†η_out S − η_in|  vs dx  (ω=0.1, X0=12):")
    print(f"   {'dx':>8} {'max|R_S|':>12}")
    dxs = [0.008, 0.004, 0.002, 0.001]; rs = []
    for dxv in dxs:
        r = R_S(om0, X0, dxv); rs.append(r); print(f"   {dxv:8.4f} {r:12.3e}")
    slope = np.polyfit(np.log(dxs), np.log(rs), 1)[0]
    print(f"   pendiente en dx = {slope:.2f}   (converge: {'si' if rs[-1] < rs[0] else 'NO'})\n")

    # --- C1: meseta en X0 ---
    print("[C1] meseta en X0 (ω=0.1, dx=0.002):")
    print(f"   {'X0':>6} {'max|R_S|':>12}")
    for X0v in [8.0, 10.0, 12.0, 14.0, 16.0, 20.0]:
        print(f"   {X0v:6.1f} {R_S(om0, X0v, 0.002):12.3e}")
    print("   OJO: el método de TRANSFERENCIA de v1.0 es estable solo en una VENTANA de X0: para X0 grande")
    print("   acumula el modo evanescente CRECIENTE (~e^{+κX0}) y R_S explota. NO hay meseta monótona; hay")
    print("   un mínimo (~X0=12) y luego blow-up. Limitación del baseline (un solver estabilizado/scattering")
    print("   directo lo evita). Es un HALLAZGO, no un pase de A3-meseta.\n")

    # --- C1: recuento de canales estable en la ventana de ω ---
    print("[C1] recuento de canales propagantes en la ventana de ω (dx=0.002, X0=12):")
    print(f"   {'ω':>8} {'#IN_der':>8} {'#OUT_der':>9} {'max|R_S|':>12}")
    for omv in [0.05, 0.1, 0.15, 0.2]:
        S2, ei, eo, inf = smatrix_at(omv, 12.0, 0.002)
        if S2 is None:
            mL, mR = inf; nin = sum(1 for m in mR if m['prop'] and m['vg'] < 0)
            nout = sum(1 for m in mR if m['prop'] and m['vg'] > 0)
            print(f"   {omv:8.3f} {nin:>8} {nout:>9} {'(no 3+3)':>12}")
        else:
            Ei = np.diag(ei.astype(float)); Eo = np.diag(eo.astype(float))
            r = np.abs(S2.conj().T@Eo@S2 - Ei).max()
            print(f"   {omv:8.3f} {2:>8} {2:>9} {r:12.3e}")
    print()

    # --- C2: balance regional con paquete, con termino de defecto (SEPARADO de C1) ---
    print("[C2] balance regional con PAQUETE (independiente de C1): ΔQ_R vs ∫(frontera+defecto) dt")
    print("     identidad: dQ_R/dt = i⟨Ψ,[M,P_R]Ψ⟩_𝒢 + iΨ†(M†𝒢−𝒢M)P_R Ψ  (defecto NO asumido 0)")
    print(f"   {'N':>5} {'dx':>7} {'ΔQ_R':>12} {'∫total':>12} {'|balance|':>11} {'rel':>10} "
          f"{'∫defecto':>11} {'recon':>9} {'‖Ψ‖_T':>8}")
    for N in [300, 400, 600]:
        o = packet_balance(N=N)
        print(f"   {N:>5} {o['dx']:7.4f} {o['dQR']:12.4e} {o['integ_tot']:12.4e} {o['balance']:11.2e} "
              f"{o['rel']:10.2e} {o['integ_def']:11.2e} {o['recon_max']:9.1e} {o['Qgrow']:8.4f}")
    print("   -> |balance|/rel → 0 al refinar = la identidad regional cierra DURANTE la evolucion;")
    print("      ∫defecto = magnitud del termino de pseudo-herm (NO despreciable, se reporta); recon =")
    print("      max|dQ_R/dt_exacto − (frontera+defecto)| por paso (consistencia de la identidad).\n")

    print("Lectura (alcance baseline λ=1, redacción ChatGPT):")
    print(" - C1 y C2 son INDEPENDIENTES y SEPARADOS. η_in≠η_out; registro de canales (k,vg,Ω,norma,Jraw).")
    print(" - C1: pseudo-unitaridad numéricamente pequeña (~5e−5<1e−3), CONVERGENCIA con dx NO DEMOSTRADA")
    print("   (floor), VENTANA de estabilidad en X0 (no meseta) → el método de TRANSFERENCIA es INSUFICIENTE")
    print("   para satisfacer A1/A3. (El residuo <1e−3 NO sustituye al criterio de convergencia.) Controles")
    print("   A/B (gate_P3_4_controls.py): réplica fiel a v1.0 (0.0) + invariante a escala evanescente →")
    print("   el floor apunta a la PROPAGACION por transferencia, no a la construcción de S.")
    print(" - C2: balance regional discreto cierra CON el término de defecto (se mide, no se asume 0).")
    print("   El estado inicial es un ESTADO DE PRUEBA localizado (no modo BdG); recon~máquina = identidad")
    print("   ALGEBRAICA, no comprobación independiente de física de flujo. Falta barrido de dt (orden temporal).")
    print(" - Baseline λ=1 sobre v1.0 como REFERENCIA; NO se modifica v1.0; NO se extiende a λ≠1 (Padé).")

if __name__ == '__main__':
    run()
