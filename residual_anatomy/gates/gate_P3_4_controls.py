"""
GATE P3.4 — controles de C1 pedidos por ChatGPT (los dos más baratos, ANTES de plantear solver nuevo):
  (A) Comparacion directa S_gate vs B.smatrix (v1.0): la matriz S COMPLETA debe coincidir a ~maquina
      (misma convencion/orden) -> valida la REPLICA del matching (no solo el residuo de pseudo-unit).
  (B) Reescalado de los autovectores EVANESCENTES (1e-6, 1, 1e6): S debe ser INVARIANTE (su normalizacion
      es arbitraria) -> si cambia, el matching esta mal condicionado o mal planteado.
No modifica v1.0 ni el gate. Importa las funciones del gate (no ejecuta su run()).
"""
import numpy as np
import gate_P3_4_smatrix as G
B = G.B

def run():
    print("GATE P3.4 — controles de C1 (fidelidad de la réplica + invariancia evanescente)\n" + "="*72)
    X0, dx = 12.0, 0.002
    oms = [0.05, 0.1, 0.15, 0.2]

    # --- (A) S_gate vs B.smatrix (v1.0) ---
    print("[A] Comparacion S_gate vs S_v1.0 (B.smatrix), X0=12, dx=0.002 — debe ~maquina:")
    print(f"   {'ω':>8} {'max|Sg−Sv1|':>13} {'η_in iguales':>13} {'η_out iguales':>14}")
    for om in oms:
        Sg, ein_g, eout_g, _ = G.smatrix_at(om, X0, dx)
        Sv, ein_v, eout_v, _ = B.smatrix(om, X0)
        dS = np.abs(Sg - Sv).max()
        eqi = bool(np.array_equal(ein_g, ein_v)); eqo = bool(np.array_equal(eout_g, eout_v))
        print(f"   {om:8.3f} {dS:13.2e} {str(eqi):>13} {str(eqo):>14}")
    print()

    # --- (B) invariancia bajo reescalado de los evanescentes ---
    print("[B] Invariancia de S al reescalar los autovectores EVANESCENTES (ω=0.1, X0=12, dx=0.002):")
    print("    S no debe depender de la normalizacion arbitraria de los evanescentes.")
    Sref, _, _, _ = G.smatrix_at(0.1, X0, dx, ev_scale=1.0)
    print(f"   {'ev_scale':>10} {'max|S−S(ev=1)|':>16} {'max|R_S|':>12}")
    for s in [1e-6, 1e-3, 1.0, 1e3, 1e6]:
        S, ein, eout, _ = G.smatrix_at(0.1, X0, dx, ev_scale=s)
        dS = np.abs(S - Sref).max()
        r = np.abs(S.conj().T @ np.diag(eout.astype(float)) @ S - np.diag(ein.astype(float))).max()
        print(f"   {s:10.0e} {dS:16.2e} {r:12.3e}")
    print()
    print("Lectura:")
    print(" - Si [A] ~maquina: la réplica del matching es FIEL a v1.0 (el floor NO viene de la réplica).")
    print(" - Si [B] ~maquina e invariante: el matching NO depende de la escala evanescente arbitraria.")
    print(" - Ambos OK -> el floor/blow-up apuntan a la PROPAGACION por transferencia (modo creciente),")
    print("   no a la construccion de S. (No prueba la causa exacta por si solo; la acota.)")

if __name__ == '__main__':
    run()
