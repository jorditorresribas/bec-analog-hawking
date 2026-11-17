"""
GATE P3.4 — C2: barrido TEMPORAL (nt), N y T fijos. Control pedido por ChatGPT (el más barato primero).
Verifica que el error de balance |ΔQ_R − ∫(frontera+defecto)dt| baja como O(dt²) (trapecio temporal).
NOTA: con el propagador exacto U=expm(−iM·dt), aplicar U nt veces = expm(−iM·T) EXACTO → Ψ(T) y por
tanto ΔQ_R son INDEPENDIENTES de nt (hasta la precisión de expm). El único efecto de nt es el error de
la CUADRATURA temporal del ∫ → el balance debe escalar O(dt²). No toca v1.0 ni el gate.
"""
import numpy as np
import gate_P3_4_smatrix as G

def run():
    print("GATE P3.4 — C2 barrido temporal nt (N=400, T=4 fijos)\n" + "="*68)
    N, T = 400, 4.0
    print(f"   {'nt':>6} {'dt':>8} {'ΔQ_R':>12} {'∫total':>12} {'|balance|':>11} {'∫defecto':>11} {'recon':>9}")
    nts = [100, 200, 400, 800, 1600]; dts = []; bal = []
    for nt in nts:
        o = G.packet_balance(N=N, T=T, nt=nt)
        dt = T/nt; dts.append(dt); bal.append(o['balance'])
        print(f"   {nt:>6} {dt:8.4f} {o['dQR']:12.6e} {o['integ_tot']:12.6e} {o['balance']:11.3e} "
              f"{o['integ_def']:11.3e} {o['recon_max']:9.1e}")
    slope = np.polyfit(np.log(dts), np.log(bal), 1)[0]
    print(f"\n   pendiente de |balance| vs dt = {slope:.2f}  (esperado ≈2 = O(dt²) del trapecio temporal)")
    print("   ΔQ_R y ∫defecto deben ser ~constantes en nt (Ψ(T)=expm(−iM·T) indep. de nt); solo cambia")
    print("   la cuadratura temporal del ∫ → el balance decae O(dt²). recon ~máquina (identidad algebraica).")

if __name__ == '__main__':
    run()
