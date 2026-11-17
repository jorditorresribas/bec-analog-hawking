"""
CONVERGENCIA HIDRODINAMICA con MODE-MATCHING (estable, alcanza w grande).
N_omega validado por 2 metodos. Aqui: para cada anchura w, fit termico de
N_omega(omega) y comparacion T_fit/T_grav. Tendencia -> 1 = dispersion.

Protocolo (ChatGPT): reconstruir fondo por w, recalcular kappa=|d(c_s-v)/dx|_H
(correcto, SIN 1/2), T_grav=kappa/2pi; N_omega en ventana comun lejos de umbrales;
fit termico documentado (ventana, puntos, residuo); pseudo-unit por frecuencia;
convergencia; NO cambiar el procedimiento de extraccion entre anchuras.
Cautela: w grande NO garantiza regimen hidrodinamico; T_fit/T_grav->1 es evidencia
para ESTE analogo y regimen, no Hawking gravitacional.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import mode_matching as MM

# ventana de ajuste ESCALADA a T_grav de cada anchura (regimen termico omega~T,
# no la cola dispersiva). Se construye en fit_T a partir de Tg.
def fit_window(Tg):
    return np.linspace(max(0.012, 0.6*Tg), 3.0*Tg, 10)


def kappa_Tgrav(w):
    xs = np.linspace(-8*w-12, 8*w+12, 8000); MM.w = w
    n, v, V = MM.bg(xs, 'horizon'); c = np.sqrt(MM.g*n)
    f = c - v; iH = np.argmin(np.abs(f))
    kappa = abs(np.gradient(f, xs)[iH])      # SIN 1/2
    return kappa, kappa/(2*np.pi)


def N_omega(w, om, L, N):
    MM.w = w; MM.omega = om
    x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
    mR = MM.modes_full(n[-1], v[-1], V[-1])
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
    out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=incR)
    Nw = sum(a for (nrm, s, a) in out if nrm > 0)         # emision norma +
    pu = abs(sum(nrm*a for (nrm, s, a) in out) - (-1.0))  # eta_in=-1 (neg-norm)
    return Nw, pu


def fit_T(w, L, N, Tg):
    oms = fit_window(Tg)
    Nw = []; pus = []
    for om in oms:
        nw, pu = N_omega(w, om, L, N); Nw.append(nw); pus.append(pu)
    Nw = np.array(Nw)
    y = np.log(1 + 1/Nw)                     # = omega/T si termico
    coef, cov = np.polyfit(oms, y, 1, cov=True)
    T_fit = 1/coef[0]; dT = T_fit**2*np.sqrt(cov[0, 0])
    resid = np.max(np.abs(y - np.polyval(coef, oms)))
    return T_fit, dT, resid, coef[1], np.max(pus)


print("="*86)
print("CONVERGENCIA HIDRODINAMICA (mode-matching, N_omega validado)")
print("  ventana de ajuste escalada a T_grav de cada w (omega ~ 0.6-3 T); fit log(1+1/N)=omega/T")
print("="*86)
print(f"{'w':>4} {'T_fit':>9} {'+-dT':>8} {'T_grav':>9} {'ratio':>7} {'resid':>8} {'eU':>8} {'intercept':>9}")
print("-"*86)
rows = []
for w in (2.0, 3.0, 4.0, 6.0, 8.0):
    L = max(22.0, 3.2*w + 8); N = int(L*200)
    kappa, Tg = kappa_Tgrav(w)
    Tf, dT, resid, intc, pu = fit_T(w, L, N, Tg)
    rows.append((w, Tf, Tg, Tf/Tg))
    print(f"{w:4.0f} {Tf:9.4f} {dT:8.1e} {Tg:9.4f} {Tf/Tg:7.3f} {resid:8.1e} {pu:8.1e} {intc:9.4f}")
print("-"*86)
print("  tendencia de T_fit/T_grav con w:")
for w, Tf, Tg, r in rows:
    print(f"    w={w:.0f}: ratio={r:.3f}")
r0, rN = rows[0][3], rows[-1][3]
trend = ("CONVERGE hacia 1 (desviacion era dispersiva)" if abs(rN-1) < abs(r0-1)*0.6
         else "no converge claramente a 1 (investigar)")
print(f"  -> w=2: {r0:.3f}  ->  w=8: {rN:.3f}   [{trend}]")

# convergencia en N a w=4 (control)
print("\n  CONVERGENCIA EN N (T_fit a w=4):")
for N in (3000, 5000):
    Tf, *_ = fit_T(4.0, 22.0, N, kappa_Tgrav(4.0)[1]); print(f"    N={N}: T_fit={Tf:.4f}")
print("\n  (T_fit/T_grav->1 = evidencia de convergencia al espectro termico de Hawking")
print("   para ESTE analogo y regimen; NO Hawking gravitacional.)")
