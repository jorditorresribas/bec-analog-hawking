"""
H3 — AUDITORIA DE CONVERGENCIA NUMERICA del mode-matching (N_omega validado).
NO asumir que L insuficiente es la causa del mal pseudo-unit en w=8 baja omega.
Separar L (tamano finito) de N/dx (discretizacion). Registrar pseudo-unit Y N_omega
Y su variacion entre resoluciones (estimador directo del error, no solo pseudo-unit).

Experimentos: A) L con dx fijo; B) N con L fijo; D) barrido baja-omega en el mejor
config; control w=2. Caso problematico: w=8 baja omega (pseudo-unit ~7.5e-2).
NO modifica scripts originales. Usa mode_matching.solve_inject.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import csv
import mode_matching as MM

g = MM.g


def Nw_pu(w, om, L, N):
    MM.w = w; MM.omega = om
    x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
    mR = MM.modes_full(n[-1], v[-1], V[-1])
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
    out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=incR)
    Nw = sum(a for (nrm, s, a) in out if nrm > 0)
    pu = abs(sum(nrm*a for (nrm, s, a) in out) - (-1.0))
    return Nw, pu


rows = []
print("="*74)
print("H3 — convergencia numerica (w=8 baja-omega problematico, w=2 control)")
print("="*74)

# --- Experimento A: L con dx fijo (=0.01). w=8, omega baja ---
print("\n[A] Convergencia en L (dx=0.01 fijo).  w=8, omega=0.008  (T_grav~0.0108)")
print(f"    {'L':>5} {'N':>7} {'N_omega':>10} {'pseudo-unit':>12} {'dN/N vs prev':>12}")
dx = 0.01; prevN = None
for L in (34, 50, 70, 90):
    N = int(2*L/dx)
    Nw, pu = Nw_pu(8.0, 0.008, L, N)
    dv = abs(Nw-prevN)/Nw if prevN else np.nan
    print(f"    {L:5.0f} {N:7d} {Nw:10.5f} {pu:12.2e} {dv:12.2e}")
    rows.append([8.0, 0.008, 'A_varyL', L, N, dx, Nw, pu])
    prevN = Nw

# --- Experimento B: N con L fijo (=60). w=8, omega baja ---
print("\n[B] Convergencia en N (L=60 fijo -> dx decreciente).  w=8, omega=0.008")
print(f"    {'N':>7} {'dx':>8} {'N_omega':>10} {'pseudo-unit':>12} {'dN/N vs prev':>12}")
prevN = None
for N in (6000, 12000, 24000):
    dx = 120.0/N
    Nw, pu = Nw_pu(8.0, 0.008, 60.0, N)
    dv = abs(Nw-prevN)/Nw if prevN else np.nan
    print(f"    {N:7d} {dx:8.4f} {Nw:10.5f} {pu:12.2e} {dv:12.2e}")
    rows.append([8.0, 0.008, 'B_varyN', 60.0, N, dx, Nw, pu])
    prevN = Nw

# --- Experimento D: barrido baja-omega w=8 en config amplio (L=80, dx=0.01) ---
print("\n[D] Barrido baja-omega w=8 (L=80, dx=0.01):")
print(f"    {'omega':>7} {'omega/Tg':>9} {'N_omega':>10} {'pseudo-unit':>12}")
Tg8 = 0.01075
for om in (0.004, 0.006, 0.008, 0.012, 0.016, 0.022):
    Nw, pu = Nw_pu(8.0, om, 80.0, 16000)
    print(f"    {om:7.4f} {om/Tg8:9.2f} {Nw:10.5f} {pu:12.2e}")
    rows.append([8.0, om, 'D_lowom', 80.0, 16000, 0.01, Nw, pu])

# --- control w=2 (debe estar ya convergido) ---
print("\n[control] w=2, omega=0.026 (~0.6 Tg): convergencia en L")
print(f"    {'L':>5} {'N':>7} {'N_omega':>10} {'pseudo-unit':>12}")
for L in (24, 40):
    N = int(2*L/0.01)
    Nw, pu = Nw_pu(2.0, 0.026, L, N)
    print(f"    {L:5.0f} {N:7d} {Nw:10.5f} {pu:12.2e}")
    rows.append([2.0, 0.026, 'ctrl_w2', L, N, 0.01, Nw, pu])

with open('out/h3_convergence.csv', 'w', newline='') as f:
    w = csv.writer(f); w.writerow(['w', 'omega', 'exp', 'L', 'N', 'dx', 'N_omega', 'pseudo_unit'])
    w.writerows(rows)
print("\nCSV: out/h3_convergence.csv")
print("\n(interpretacion: mejora con L=tamano finito; con N=discretizacion;")
print(" si N_omega converge aunque pseudo-unit sea alto -> revisar extraccion.)")
