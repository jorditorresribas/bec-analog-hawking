"""
H1 — SENSIBILIDAD DE T_fit A LA VENTANA DE AJUSTE (caracterizar el ~10%).
NO modifica fisica ni scripts originales. N_omega = el VALIDADO por mode-matching
(mode_matching.solve_inject), NO los scripts viejos del transfer.

Pregunta: T_fit = T_fit(omega_min, omega_max)? Si T_fit/T_grav cambia mucho con la
ventana -> parte del ~10% es sesgo del ajuste (H1). Si es estable -> H1 debil.
Ajuste Planck: log(1+1/N_omega) = omega/T  (misma funcion que el analisis original).
Guarda CSV + grafica. NO redefine kappa ni introduce factores.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import csv
import mode_matching as MM
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g = MM.g; n_u = 1.0; mach_u = 0.5; md = 2.0
vLb = mach_u*np.sqrt(g*n_u); jc = n_u*vLb; mu = g*n_u + 0.5*vLb**2
n_d = (jc/(md*np.sqrt(g)))**(2/3)


def T_grav_geom(w):
    xs = np.linspace(-8*w-12, 8*w+12, 20000); MM.w = w
    t = np.tanh(xs/w); s = 0.5*(1+t); n = n_u+(n_d-n_u)*s; v = jc/n; c = np.sqrt(g*n)
    f = c - v; iH = np.argmin(np.abs(f))
    kappa = abs(np.gradient(c, xs)[iH] - np.gradient(v, xs)[iH])  # |d(c-v)/dx|_H (geom)
    return kappa/(2*np.pi)


def N_omega(w, om, L, N):
    MM.w = w; MM.omega = om
    x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
    mR = MM.modes_full(n[-1], v[-1], V[-1])
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
    out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=incR)
    Nw = sum(a for (nrm, s, a) in out if nrm > 0)
    pu = abs(sum(nrm*a for (nrm, s, a) in out) - (-1.0))
    return Nw, pu


def planck_fit(oms, Nw):
    y = np.log(1 + 1/Nw); coef = np.polyfit(oms, y, 1)
    T = 1/coef[0]; resid = np.sqrt(np.mean((y - np.polyval(coef, oms))**2))
    return T, coef[1], resid


print("="*78)
print("H1 — sensibilidad de T_fit a la ventana (N_omega mode-matching validado)")
print("="*78)

rows_csv = []
for w, L, N in ((2.0, 24.0, 4800), (8.0, 34.0, 6800)):
    Tg = T_grav_geom(w)
    # malla densa de N_omega cubriendo ~0.4 a 4.5 * Tg
    om_dense = np.linspace(0.4*Tg, 4.5*Tg, 24)
    Nw = np.array([N_omega(w, om, L, N)[0] for om in om_dense])
    pus = np.array([N_omega(w, om, L, N)[1] for om in om_dense])
    print(f"\n--- w={w:.0f}  T_grav_geom={Tg:.5f}  (malla N_omega: {om_dense[0]:.4f}..{om_dense[-1]:.4f}) ---")
    print(f"  pseudo-unit en la malla: max {pus.max():.1e}  (baja omega: {pus[0]:.1e})")

    # BASELINE: ventana usada en hydro_convergence [0.6 Tg, 3 Tg] -> reproducir
    base = (om_dense >= 0.6*Tg) & (om_dense <= 3.0*Tg)
    Tb, ib, rb = planck_fit(om_dense[base], Nw[base])
    print(f"  BASELINE [0.6Tg,3Tg]: T_fit={Tb:.5f} ratio={Tb/Tg:.3f} (reproduce hydro_convergence)")

    # (a) variar omega_min (omega_max fijo = 3.2 Tg)
    print(f"  variar omega_min  (omega_max=3.2 Tg):")
    print(f"    {'wmin/Tg':>8} {'npts':>5} {'T_fit':>9} {'ratio':>7} {'intc':>8} {'rmsRes':>8}")
    omax = 3.2*Tg
    for fmin in (0.4, 0.6, 0.8, 1.0, 1.3, 1.6):
        sel = (om_dense >= fmin*Tg) & (om_dense <= omax)
        if sel.sum() < 4: continue
        T, ic, r = planck_fit(om_dense[sel], Nw[sel])
        print(f"    {fmin:8.1f} {sel.sum():5d} {T:9.5f} {T/Tg:7.3f} {ic:8.4f} {r:8.1e}")
        rows_csv.append([w, Tg, 'vary_min', fmin, 3.2, sel.sum(), T, T/Tg, ic, r])
    # (b) variar omega_max (omega_min fijo = 0.6 Tg)
    print(f"  variar omega_max  (omega_min=0.6 Tg):")
    print(f"    {'wmax/Tg':>8} {'npts':>5} {'T_fit':>9} {'ratio':>7} {'intc':>8} {'rmsRes':>8}")
    omin = 0.6*Tg
    for fmax in (1.8, 2.4, 3.0, 3.6, 4.2):
        sel = (om_dense >= omin) & (om_dense <= fmax*Tg)
        if sel.sum() < 4: continue
        T, ic, r = planck_fit(om_dense[sel], Nw[sel])
        print(f"    {fmax:8.1f} {sel.sum():5d} {T:9.5f} {T/Tg:7.3f} {ic:8.4f} {r:8.1e}")
        rows_csv.append([w, Tg, 'vary_max', 0.6, fmax, sel.sum(), T, T/Tg, ic, r])

with open('out/temperature_window_audit.csv', 'w', newline='') as fcsv:
    wtr = csv.writer(fcsv)
    wtr.writerow(['w', 'Tgrav', 'tipo', 'fmin', 'fmax', 'npts', 'Tfit', 'ratio', 'intercept', 'rms_resid'])
    wtr.writerows(rows_csv)
print("\nCSV: out/temperature_window_audit.csv")

# grafica: ratio vs ventana
fig, ax = plt.subplots(1, 2, figsize=(12, 4.5))
for w in (2.0, 8.0):
    d = [r for r in rows_csv if r[0] == w and r[2] == 'vary_min']
    ax[0].plot([r[3] for r in d], [r[7] for r in d], 'o-', label=f'w={w:.0f}')
    d2 = [r for r in rows_csv if r[0] == w and r[2] == 'vary_max']
    ax[1].plot([r[4] for r in d2], [r[7] for r in d2], 's-', label=f'w={w:.0f}')
ax[0].axhline(1, color='gray', ls=':'); ax[0].set_xlabel('omega_min / T_grav')
ax[0].set_ylabel('T_fit / T_grav_geom'); ax[0].set_title('sensibilidad a omega_min'); ax[0].legend()
ax[1].axhline(1, color='gray', ls=':'); ax[1].set_xlabel('omega_max / T_grav')
ax[1].set_ylabel('T_fit / T_grav_geom'); ax[1].set_title('sensibilidad a omega_max'); ax[1].legend()
plt.tight_layout(); plt.savefig('out/temperature_window_audit.png', dpi=120)
print("Figura: out/temperature_window_audit.png")
