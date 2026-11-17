"""
RUMBO #1 (b) — OPCION A (recomendada por ChatGPT): DEPENDENCIA CON w.

Pregunta fisica: la senal de Hawking (la lengua) se vuelve MENOS sensible a la
dispersion UV (lambda) cuando el horizonte se SUAVIZA (w grande -> mas
hidrodinamico)? Si la variacion-lambda de la lengua DECRECE con w, es evidencia
de que la robustez trans-planckiana es FISICA (Unruh-Schuetzhold/Jacobson), no una
particularidad numerica de w=2.

DISCIPLINA (ChatGPT): se cambia SOLO w. Observable, criterios de validacion,
metrica y controles IDENTICOS a los pasos previos. Metrica = integral de la parte
negativa de C_d2 en el cuadrante cruzado (la que convergio limpio en N y L).
Dentro de cada w se normaliza a su propio lambda=1 -> se compara la DISPERSION del
cociente entre w's.

w = {2, 4, 6}.  L, N, ncoarse fijos.  lambda = {0.5, 1, 2} (rango x4 convergido).
Se reporta kappa(w)=|d(c_s-v)/dx| en el horizonte (gravedad superficial, ~1/w).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import correlations as C

L, N, ncoarse = 30.0, 3000, 360
ws = [2.0, 4.0, 6.0]
lams = [0.5, 1.0, 2.0]
g = C.g

def omega_max(lam):
    """Borde de la ventana de Hawking: propiedad ASINTOTICA (Mach_d=2 fijo),
    independiente de w. Se calcula una vez por lam."""
    lo, hi = 0.05, 1.3
    for _ in range(26):
        mid = 0.5*(lo+hi)
        r = C.scatter_field(L=24, N=2000, omega=mid, channel='d2', lam=lam,
                            mode='horizon', op='lambda')
        lo, hi = (mid, hi) if (r['Rin_neg'] is not None and r['nchan'] == 2) else (lo, mid)
    return lo

def kappa(w):
    """Gravedad superficial acustica geometrica kappa=|d(c_s - v)/dx| en el horizonte
    (donde v = c_s). Derivada validada en kappa_derive.py (sin 1/2)."""
    x = np.linspace(-L, L, N); n, v, _ = C.bg(x, w, 'horizon')
    cs = np.sqrt(g*n)
    f = cs - v
    df = np.gradient(f, x)
    ih = np.argmin(np.abs(cs - v))          # horizonte: c_s = v
    return abs(df[ih]), x[ih]

def tongue_strength(w, lam, nomega=90):
    omax = OMAX[lam]
    omegas = np.linspace(0.02, 0.9*omax, nomega)
    xc, Cd2, _ = C.density_correlation(L, N, omegas, w=w, lam=lam, mode='horizon',
                                       ncoarse=ncoarse, terms='d2', op='lambda')
    dx = xc[1]-xc[0]
    iu = (xc < -0.5) & (xc > -14); idn = (xc > 0.5) & (xc < 14)
    blk = Cd2[np.ix_(iu, idn)]
    return blk[blk < 0].sum()*dx**2

# sanity: canal Hawking presente a w grande
print("Sanity (lam=1, omega=0.1): nchan y flujo por w")
for w in ws:
    r = C.scatter_field(L, N, 0.1, 'd2', w=w, lam=1.0, mode='horizon', op='lambda')
    print(f"  w={w}: nchan={r['nchan']}  flux_out={r['flux_out']:+.4f}  eta_in={r['eta_in']}")

OMAX = {lam: omega_max(lam) for lam in lams}
print(f"\nomega_max(lam) (asintotico, indep. de w): " +
      "  ".join(f"lam={l}:{OMAX[l]:.3f}" for l in lams))
print("\nkappa(w) (gravedad superficial, ~1/w):")
for w in ws:
    k, xh = kappa(w)
    print(f"  w={w}: kappa={k:.4f}  (x_H={xh:+.3f})  T_H=kappa/2pi={k/(2*np.pi):.4f}")

print("\n" + "="*64)
print("FUERZA DE LA LENGUA str(w,lam) y cociente str/str(lam=1) por w")
print("="*64)
S = {}
for w in ws:
    S[w] = {lam: tongue_strength(w, lam) for lam in lams}
print(f"{'w':>4} | " + " ".join(f"str(lam={l})" for l in lams) +
      " | " + " ".join(f"r(lam={l})" for l in lams) + " | spread")
spreads = {}
for w in ws:
    base = S[w][1.0]
    r = {lam: S[w][lam]/base for lam in lams}
    spread = max(abs(r[0.5]-1), abs(r[2.0]-1))     # desviacion maxima vs lam=1
    spreads[w] = spread
    print(f"{w:>4} | " + " ".join(f"{S[w][l]:+.3e}" for l in lams) +
          " | " + " ".join(f"{r[l]:.3f}" for l in lams) + f" | {spread:.3f}")

print("\n" + "="*64)
print("HIPOTESIS: la DISPERSION del cociente (sensibilidad a lambda) debe DECRECER")
print("al crecer w (horizonte mas hidrodinamico -> mas robusto frente al UV k^4).")
print("="*64)
for w in ws:
    print(f"  w={w}: spread (max|r-1|) = {spreads[w]:.3f}  "
          f"= sensibilidad a lambda sobre x4")
trend = "DECRECE (confirma robustez fisica)" if spreads[6.0] < spreads[2.0] \
        else "NO decrece (revisar)"
print(f"\n  Tendencia w=2 -> w=6: {spreads[2.0]:.3f} -> {spreads[6.0]:.3f}  => {trend}")

np.savez("out/w_dependence.npz", ws=np.array(ws), lams=np.array(lams),
         kappa=np.array([kappa(w)[0] for w in ws]),
         S=np.array([[S[w][l] for l in lams] for w in ws]),
         spreads=np.array([spreads[w] for w in ws]))

fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
cmap = plt.cm.plasma(np.linspace(0.1, 0.8, len(ws)))
for col, w in zip(cmap, ws):
    base = S[w][1.0]
    ax[0].plot(lams, [S[w][l]/base for l in lams], 'o-', color=col,
               label=f"w={w} (kappa={kappa(w)[0]:.3f})")
ax[0].axhline(1, color='k', lw=0.5, alpha=0.4); ax[0].set_xscale('log', base=2)
ax[0].set_xlabel("lambda (dispersion UV k^4)")
ax[0].set_ylabel("fuerza lengua / lam=1 (dentro de cada w)")
ax[0].set_title("Sensibilidad a lambda por anchura de horizonte w")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
ax[1].plot(ws, [spreads[w] for w in ws], 'o-', color='crimson')
ax[1].set_xlabel("w (anchura del horizonte)")
ax[1].set_ylabel("spread = max|r-1| sobre lambda x4")
ax[1].set_title("Sensibilidad UV vs w  (hipotesis: DECRECE)")
ax[1].grid(alpha=0.3)
fig.savefig("out/w_dependence.png", dpi=120)
print("\nFig -> out/w_dependence.png   Datos -> out/w_dependence.npz")
