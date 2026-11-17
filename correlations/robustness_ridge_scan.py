"""
RUMBO #1 (b) — DIAGNOSTICO: la metrica limpia (envolvente) usaba pendiente FIJA
-1.26, pero la pendiente del compañero depende de v_g, que cambia con lam. Si la
cresta se desplaza, la envolvente a pendiente fija subestima -> ARTEFACTO.

Aqui, para cada lam, ESCANEAMOS la pendiente del ridge: para un abanico de
pendientes integramos la envolvente de Hilbert a lo largo de esa recta; la
pendiente que MAXIMIZA la envolvente integrada ES la cresta real. Reportamos:
  - slope_opt(lam)         (pendiente real de la cresta de Hawking)
  - env_int_opt(lam)       (fuerza en la cresta real: limpia Y sin sesgo de pendiente)
  - slope predicho por v_g de los 2 compañeros a cada lam (prediccion independiente)
Asi separamos FISICA (la cresta se mueve y/o debilita) de ARTEFACTO de muestreo.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
from scipy.signal import hilbert
from scipy.interpolate import RegularGridInterpolator
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import correlations as C

L, N, w = 30.0, 3000, 2.0
ncoarse = 360
lams = [0.25, 0.5, 1.0, 2.0, 4.0]
g = C.g

# fondo asintotico (para v_g de los compañeros)
n_u = C.n_u; vu = C.vLb; cu = np.sqrt(g*n_u)
n_d = (C.jcur/(2.0*np.sqrt(g)))**(2/3); vd = C.jcur/n_d; cd = np.sqrt(g*n_d)

def omega_max(lam):
    lo, hi = 0.05, 1.3
    for _ in range(26):
        mid = 0.5*(lo+hi)
        r = C.scatter_field(L=24, N=2000, omega=mid, channel='d2', lam=lam,
                            mode='horizon', op='lambda')
        lo, hi = (mid, hi) if (r['Rin_neg'] is not None and r['nchan'] == 2) else (lo, mid)
    return lo

def vg_slope_pred(lam, omega=0.05):
    """Pendiente del compañero = v_g(upstream Hawking out) / v_g(downstream partner),
    medidos con la dispersion lam a baja omega (limite IR de la correlacion).
    Devuelve x_d/x_u = v_g(down)/v_g(up) con signos (upstream vg<0 sale a -inf,
    downstream partner vg<0 sale a -inf; en el ridge x'=slope x)."""
    # usamos la clasificacion de canales del propio solver a esta omega
    x = np.linspace(-L, L, N); n, v, V = C.bg(x, w, 'horizon')
    mL = C._classify(n[0], v[0], V[0], cu, omega, lam, 'lambda')
    mR = C._classify(n[-1], v[-1], V[-1], cd, omega, lam, 'lambda')
    # Hawking out upstream: propagante vg<0 (sale a x->-inf) en upstream
    up = [m for m in mL if m['prop'] and m['vg'] < 0]
    dn = [m for m in mR if m['prop'] and m['vg'] < 0]   # compañeros que salen a -inf? no
    # compañero lento downstream: el de norma - (d2) OUT. Tomamos vg de cada
    vgu = up[0]['vg'] if up else np.nan
    # partner: canal downstream norma negativa saliente (vg<0 en down = hacia horizonte? )
    vgd = min((m['vg'] for m in mR if m['prop']), key=abs) if mR else np.nan
    return vgu, vgd

def ridge_env_int(xc, Cmat, slope, smax=13.0, ns=400):
    interp = RegularGridInterpolator((xc, xc), Cmat, bounds_error=False, fill_value=0.0)
    s = np.linspace(0.5, smax, ns)
    g_ = interp(np.column_stack([-s, abs(slope)*s]))
    E = np.abs(hilbert(g_))
    return np.trapezoid(E, s), E, s

def compute_Cd2(lam):
    omax = omega_max(lam)
    omegas = np.linspace(0.02, 0.9*omax, 90)
    xc, Cd2, _ = C.density_correlation(L, N, omegas, w=w, lam=lam, mode='horizon',
                                       ncoarse=ncoarse, terms='d2', op='lambda')
    return omax, xc, Cd2

slopes = np.linspace(0.8, 2.4, 49)   # magnitudes de pendiente a escanear
print("Escaneo de pendiente de la cresta (envolvente integrada) por lam:")
print(f"{'lam':>5} {'omax':>6} | {'slope_opt':>9} {'env@opt':>9} {'env@-1.26':>10} "
      f"{'opt/fix':>7}")
R = {}
for lam in lams:
    omax, xc, Cd2 = compute_Cd2(lam)
    vals = np.array([ridge_env_int(xc, Cd2, -sl)[0] for sl in slopes])
    jopt = int(np.argmax(vals)); slope_opt = -slopes[jopt]; env_opt = vals[jopt]
    env_fix = ridge_env_int(xc, Cd2, -1.26)[0]
    R[lam] = dict(omax=omax, xc=xc, Cd2=Cd2, slope_opt=slope_opt,
                  env_opt=env_opt, env_fix=env_fix, slopes=-slopes, vals=vals)
    print(f"{lam:>5} {omax:>6.3f} | {slope_opt:>9.2f} {env_opt:>9.3e} {env_fix:>10.3e} "
          f"{env_opt/env_fix:>7.2f}")

print("\n--- ROBUSTEZ con metrica limpia EN LA CRESTA REAL (sin sesgo de pendiente) ---")
base = R[1.0]['env_opt']
print(f"{'lam':>5} | {'slope_opt':>9} | {'env_opt/lam1':>12}")
for lam in lams:
    print(f"{lam:>5} | {R[lam]['slope_opt']:>9.2f} | {R[lam]['env_opt']/base:>12.3f}")

np.savez("out/robustness_ridge_scan.npz", lams=np.array(lams),
         slope_opt=np.array([R[l]['slope_opt'] for l in lams]),
         env_opt=np.array([R[l]['env_opt'] for l in lams]),
         env_fix=np.array([R[l]['env_fix'] for l in lams]))

fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
cmap = plt.cm.viridis(np.linspace(0, 1, len(lams)))
for col, lam in zip(cmap, lams):
    ax[0].plot(R[lam]['slopes'], R[lam]['vals'], color=col, lw=1.6, label=f"lam={lam}")
    ax[0].axvline(R[lam]['slope_opt'], color=col, lw=0.6, ls=':')
ax[0].axvline(-1.26, color='k', lw=1.0, ls='--', label='pred -1.26 (lam=1)')
ax[0].set_xlabel("pendiente de la recta de muestreo")
ax[0].set_ylabel("envolvente integrada")
ax[0].set_title("Escaneo de pendiente: la cresta real maximiza la envolvente")
ax[0].legend(fontsize=8)
ax[1].plot(lams, [R[l]['env_opt']/base for l in lams], 'o-', label='env en cresta real')
ax[1].plot(lams, [R[l]['env_fix']/R[1.0]['env_fix'] for l in lams], 's--',
           label='env a pendiente fija -1.26', alpha=0.6)
ax[1].axhline(1, color='k', lw=0.5, alpha=0.4); ax[1].set_xscale('log', base=2)
ax[1].set_xlabel("lambda"); ax[1].set_ylabel("fuerza / lam=1")
ax[1].set_title("Robustez: cresta real vs pendiente fija")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
fig.savefig("out/robustness_ridge_scan.png", dpi=120)
print("\nFig -> out/robustness_ridge_scan.png")
