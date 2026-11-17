"""
H4.1 — AUDITORIA DE ROBUSTEZ del limite IR (sin tocar fisica/solvers/kappa).
Usa los datos guardados de H4 (out/h4_ir.csv) + convergencia nueva en N.
NO recalcular para obtener un valor mas favorable. NO elegir ventana por cercania
a T_grav. Si algo falla, documentar y reducir el alcance.

A. Reproducibilidad: reproducir los 6 ratios desde el CSV.
B. Robustez: leave-one-out + variar extremo superior de la ventana IR.
C. Convergencia: N_omega (no solo omega*N) vs N para puntos IR.
D. Presupuesto de incertidumbre (sin sumar como independientes sin justificar).
E. Conclusion.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import csv
import mode_matching as MM

# --- A. cargar datos H4 ---
data = {}
with open('out/h4_ir.csv') as f:
    for row in csv.DictReader(f):
        w = float(row['w']); data.setdefault(w, []).append(
            (float(row['omega']), float(row['N_omega']), float(row['omegaN']),
             float(row['pseudo_u']), float(row['Tgrav']), float(row['om_over_Tg'])))

print("="*78)
print("H4.1 — auditoria de robustez del limite IR")
print("="*78)


def extrap(oms, omN, deg):
    c = np.polyfit(oms, omN, deg); return c[-1]  # intercepto = T_IR


def planck(oms, Nw):
    y = np.log(1+1/Nw); c = np.polyfit(oms, y, 1); return 1/c[0]


for w in sorted(data):
    d = sorted(data[w]); oms = np.array([r[0] for r in d])
    Nw = np.array([r[1] for r in d]); omN = np.array([r[2] for r in d])
    pu = np.array([r[3] for r in d]); Tg = d[0][4]; xr = np.array([r[5] for r in d])
    print(f"\n{'='*40}\n w={w:.0f}  T_grav={Tg:.5f}  ({len(oms)} puntos, om/Tg in [{xr.min():.2f},{xr.max():.2f}])")
    # A. reproducir
    selP = xr <= 0.6
    print(f"  [A] Reproducir: Planck(om/Tg<0.6)={planck(oms[selP],Nw[selP])/Tg:.3f}  "
          f"IR_lin={extrap(oms,omN,1)/Tg:.3f}  IR_quad={extrap(oms,omN,2)/Tg:.3f}")

    # B. leave-one-out en la extrapolacion
    loo_lin = []; loo_quad = []
    for i in range(len(oms)):
        m = np.ones(len(oms), bool); m[i] = False
        loo_lin.append(extrap(oms[m], omN[m], 1)/Tg)
        loo_quad.append(extrap(oms[m], omN[m], 2)/Tg)
    print(f"  [B] Leave-one-out: IR_lin ratio in [{min(loo_lin):.3f},{max(loo_lin):.3f}] "
          f"(spread {max(loo_lin)-min(loo_lin):.3f}); IR_quad in [{min(loo_quad):.3f},{max(loo_quad):.3f}]")
    # variar extremo superior de la ventana
    print(f"  [B] Variar extremo superior (om/Tg max):")
    for xmax in (0.5, 0.65, 0.82, 1.0):
        m = xr <= xmax
        if m.sum() < 3: continue
        print(f"      <= {xmax:.2f} ({m.sum()} pts): IR_lin={extrap(oms[m],omN[m],1)/Tg:.3f} "
              f"IR_quad={extrap(oms[m],omN[m],2)/Tg:.3f} Planck={planck(oms[m],Nw[m])/Tg:.3f}")

# C. convergencia en N de puntos IR (nuevo calculo, w=2 y w=8 en 2 omegas bajas)
print(f"\n{'='*40}\n  [C] Convergencia en N de N_omega (puntos IR, L fijo):")
def Nw_pu(w, om, L, N):
    MM.w = w; MM.omega = om
    x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
    mR = MM.modes_full(n[-1], v[-1], V[-1]); Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
    out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=incR)
    return sum(a for (nrm, s, a) in out if nrm > 0), abs(sum(nrm*a for (nrm, s, a) in out)+1.0)
for w, L, om in ((2.0, 60.0, 0.006), (8.0, 100.0, 0.0015)):
    print(f"    w={w:.0f} omega={om}:")
    for N in (int(2*L/0.02), int(2*L/0.01), int(2*L/0.005)):
        nw, p = Nw_pu(w, om, L, N)
        print(f"      N={N:6d} (dx={2*L/N:.4f}): N_omega={nw:.5f}  omega*N={om*nw:.6f}  pu={p:.1e}")

print(f"\n{'='*40}\n  [D] Presupuesto de incertidumbre (ratio T/T_grav):")
print("    - modelo extrapolacion (lin vs quad): ~1.7-1.8%")
print("    - leave-one-out (seleccion de puntos): ver [B] arriba")
print("    - variar ventana superior: ver [B] arriba")
print("    - convergencia N: ver [C] arriba")
print("    - tamano finito L: N_omega estable 1e-4 (H3, L convergido)")
print("  (combinar la mayor dominante, NO sumar como independientes sin justificar)")
