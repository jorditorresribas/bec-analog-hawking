"""
H4 — AUDITORIA DEL LIMITE INFRARROJO omega->0. DOS estimadores independientes:
  (1) ajuste Planck en ventanas desplazadas a baja omega.
  (2) limite IR: T_IR = lim_{omega->0} omega*N_omega  (exacto para Planck).
NO imponer T_IR=T_grav. NO diseñar para que converja. Dejar que los datos digan
si NO converge. Solo extrapolar puntos con convergencia numerica demostrada.

Para Planck ideal N=1/(e^{om/T}-1): omega*N_omega -> T. Si el espectro no es
exactamente termico, los dos estimadores pueden converger a valores distintos.
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


def Tgrav(w):
    xs = np.linspace(-8*w-12, 8*w+12, 20000); MM.w = w
    t = np.tanh(xs/w); s = 0.5*(1+t); nd = (0.5/2.0)**(2/3); n = 1+(nd-1)*s
    vv = 0.5/n; c = np.sqrt(n); f = c-vv; iH = np.argmin(np.abs(f))
    return abs(np.gradient(c, xs)[iH]-np.gradient(vv, xs)[iH])/(2*np.pi)


rows = []
print("="*80)
print("H4 — limite infrarrojo omega->0 (dos estimadores, sin imponer T_grav)")
print("="*80)

for w, L, N in ((2.0, 60.0, 12000), (8.0, 100.0, 20000)):
    Tg = Tgrav(w)
    xr = np.array([0.14, 0.20, 0.28, 0.38, 0.50, 0.65, 0.82, 1.0])   # omega/Tg
    oms = xr*Tg
    data = [Nw_pu(w, om, L, N) for om in oms]
    Nw = np.array([d[0] for d in data]); pu = np.array([d[1] for d in data])
    omN = oms*Nw                                     # omega*N_omega -> T_IR
    print(f"\n--- w={w:.0f} (L={L:.0f}, N={N})  T_grav={Tg:.5f} ---")
    print(f"  {'om/Tg':>7} {'N_omega':>10} {'omega*N':>10} {'(omN)/Tg':>9} {'pseudo-u':>9}")
    for i in range(len(oms)):
        print(f"  {xr[i]:7.2f} {Nw[i]:10.4f} {omN[i]:10.5f} {omN[i]/Tg:9.4f} {pu[i]:9.1e}")
        rows.append([w, Tg, xr[i], oms[i], Nw[i], omN[i], pu[i]])

    # convergencia numerica del punto MAS BAJO (variar L)
    om_low = oms[0]
    print(f"  convergencia del punto mas bajo (om/Tg={xr[0]:.2f}) variando L:")
    for Lc in (L, L*1.5):
        nwc, puc = Nw_pu(w, om_low, Lc, int(Lc*N/L))
        print(f"    L={Lc:.0f}: omega*N={om_low*nwc:.5f}  pseudo-u={puc:.1e}")

    # estimador 1: ajuste Planck (ventana baja om/Tg in [0.14,0.6])
    sel = xr <= 0.6
    y = np.log(1+1/Nw[sel]); cf = np.polyfit(oms[sel], y, 1); T_pl = 1/cf[0]
    # estimador 2: extrapolacion de omega*N a om->0 (lineal y cuadratica)
    cl = np.polyfit(oms, omN, 1); T_lin = cl[-1]
    cq = np.polyfit(oms, omN, 2); T_quad = cq[-1]
    print(f"  ESTIMADOR 1 (Planck, om/Tg<0.6): T={T_pl:.5f}  ratio={T_pl/Tg:.3f}")
    print(f"  ESTIMADOR 2 (omega*N -> 0, lineal):    T_IR={T_lin:.5f}  ratio={T_lin/Tg:.3f}")
    print(f"  ESTIMADOR 2 (omega*N -> 0, cuadratica): T_IR={T_quad:.5f}  ratio={T_quad/Tg:.3f}")
    print(f"  dependencia del modelo (lin vs quad): {abs(T_lin-T_quad)/T_quad:.1%}")

with open('out/h4_ir.csv', 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(['w', 'Tgrav', 'om_over_Tg', 'omega', 'N_omega', 'omegaN', 'pseudo_u'])
    wr.writerows(rows)
print("\nCSV: out/h4_ir.csv")
print("\n  Diagnostico: si ambos estimadores -> T_grav = escenario A (era ventana).")
print("  Si -> otra T comun = escenario B (correccion fisica). Si discrepan = C.")
