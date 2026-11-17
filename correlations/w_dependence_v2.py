"""
RUMBO #1 (b) — OPCION A corregida: DEPENDENCIA CON w con BANDA ESCALADA a T_H(w).

El v1 (w_dependence.py) uso banda FIJA [0.02, 0.9 omega_max] para todos los w. Pero
T_H ~ kappa/2pi cae con w (0.043 w2 -> 0.014 w6): la banda fija cae en la COLA
DISPERSIVA (omega >> T_H) para w grande, justo donde el k^4 muerde mas -> inflaba
artificialmente la sensibilidad a lambda. ES EL MISMO BUG QUE H1 de la auditoria
termica (ventana fija con T_grav que encoge).

CORRECCION: banda en unidades de T_H(w), IGUAL para todos los w: omega in
[a,b]*T_H(w). Asi la lengua se construye con los MISMOS modos termicos relativos en
cada w -> comparacion apples-to-apples de la sensibilidad a lambda.
Ademas: L crece con w (leccion H3: baja omega + horizonte ancho necesita caja grande)
y se CONTROLA la pseudo-unitaridad por frecuencia (gate: solo fiarse si <1e-2).

Se cambia SOLO w (y los controles numericos L,N que w exige). Observable y metrica
identicos. lambda = {0.5,1,2}.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import correlations as C

ws = [2.0, 4.0, 6.0]
lams = [0.5, 1.0, 2.0]
g = C.g
# ventana termica comun (en unidades de T_H). Baja omega domina la lengua; cota
# inferior no demasiado baja para no exigir L enorme, superior en la zona termica.
A_TH, B_TH = 1.0, 4.0
NOMEGA = 60

def omega_max(lam):
    lo, hi = 0.05, 1.3
    for _ in range(26):
        mid = 0.5*(lo+hi)
        r = C.scatter_field(L=24, N=2000, omega=mid, channel='d2', lam=lam,
                            mode='horizon', op='lambda')
        lo, hi = (mid, hi) if (r['Rin_neg'] is not None and r['nchan'] == 2) else (lo, mid)
    return lo

def kappa_TH(w):
    x = np.linspace(-40, 40, 4000); n, v, _ = C.bg(x, w, 'horizon')
    cs = np.sqrt(g*n); df = np.gradient(cs - v, x)
    ih = np.argmin(np.abs(cs - v))
    k = abs(df[ih]); return k, k/(2*np.pi)

def LN_for(w):
    """Caja y malla que w exige (H3): L crece con w; dx fino ~0.02 fijo."""
    L = max(40.0, 8*w + 16)       # w2->40, w4->48, w6->64
    N = int(round(2*L/0.02))
    return L, N

def pseudo_unit_defect(L, N, w, lam, omega):
    """max|flux_out - eta_in| sobre los 3 canales entrantes a esta omega (gate)."""
    d = 0.0
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L, N, omega, ch, w, lam, 'horizon', 'lambda')
        if r['Rin_neg'] is None or r['nchan'] != 2:
            return None
        d = max(d, abs(r['flux_out'] - r['eta_in']))
    return d

def tongue_strength(w, lam):
    L, N = LN_for(w)
    _, TH = kappa_TH(w)
    omax = OMAX[lam]
    olo, ohi = A_TH*TH, min(B_TH*TH, 0.9*omax)
    omegas = np.linspace(olo, ohi, NOMEGA)
    # gate de pseudo-unitaridad en 3 puntos de la banda
    pud = [pseudo_unit_defect(L, N, w, lam, om) for om in
           (omegas[0], omegas[len(omegas)//2], omegas[-1])]
    pud = [p for p in pud if p is not None]
    pu = max(pud) if pud else np.nan
    xc, Cd2, _ = C.density_correlation(L, N, omegas, w=w, lam=lam, mode='horizon',
                                       ncoarse=360, terms='d2', op='lambda')
    dx = xc[1]-xc[0]
    # ventana espacial escalada a w (horizonte mas ancho -> lengua mas extendida)
    xhi = 14.0 + 2.0*w
    iu = (xc < -0.5) & (xc > -xhi); idn = (xc > 0.5) & (xc < xhi)
    blk = Cd2[np.ix_(iu, idn)]
    return blk[blk < 0].sum()*dx**2, pu, (olo, ohi), (L, N)

OMAX = {lam: omega_max(lam) for lam in lams}
print("Banda termica comun: omega in [%.1f, %.1f]*T_H(w)" % (A_TH, B_TH))
print("w |  kappa   T_H   |  L    N   | banda omega(lam=1)   | pseudo-unit(lam=1)")
INFO = {}
for w in ws:
    k, TH = kappa_TH(w); L, N = LN_for(w)
    INFO[w] = (k, TH, L, N)

print("\n" + "="*72)
print("FUERZA DE LA LENGUA con banda escalada a T_H(w)  (metrica integral, d2)")
print("="*72)
S = {}; PU = {}; BAND = {}
for w in ws:
    S[w] = {}; PU[w] = {}
    for lam in lams:
        s, pu, band, (L, N) = tongue_strength(w, lam)
        S[w][lam] = s; PU[w][lam] = pu
        if lam == 1.0: BAND[w] = (band, L, N)
    k, TH, L, N = INFO[w]
    print(f"w={w}: kappa={k:.3f} T_H={TH:.4f}  L={L:.0f} N={N}  "
          f"banda(l1)=[{BAND[w][0][0]:.3f},{BAND[w][0][1]:.3f}]  "
          f"pu(l1)={PU[w][1.0]:.1e}")

print(f"\n{'w':>4} | " + " ".join(f"str(l={l})" for l in lams) +
      " | " + " ".join(f"r(l={l})" for l in lams) + " | spread | max pseudo-unit")
spreads = {}
for w in ws:
    base = S[w][1.0]; r = {l: S[w][l]/base for l in lams}
    spread = max(abs(r[0.5]-1), abs(r[2.0]-1)); spreads[w] = spread
    pumax = max(PU[w][l] for l in lams)
    print(f"{w:>4} | " + " ".join(f"{S[w][l]:+.3e}" for l in lams) +
          " | " + " ".join(f"{r[l]:.3f}" for l in lams) +
          f" | {spread:.3f}  | {pumax:.1e}")

print("\n" + "="*72)
print("HIPOTESIS (ChatGPT/A): spread (sensibilidad a lambda) DECRECE al crecer w.")
print("="*72)
for w in ws:
    flag = "" if max(PU[w][l] for l in lams) < 1e-2 else "  <-- PU alto, NO fiable"
    print(f"  w={w}: spread={spreads[w]:.3f}{flag}")
reliable = [w for w in ws if max(PU[w][l] for l in lams) < 1e-2]
if len(reliable) >= 2:
    trend = "DECRECE (confirma)" if spreads[reliable[-1]] < spreads[reliable[0]] \
            else ("CRECE (contra la hipotesis)" if spreads[reliable[-1]] > spreads[reliable[0]]
                  else "plano")
    print(f"\n  Tendencia (w fiables {reliable[0]}->{reliable[-1]}): "
          f"{spreads[reliable[0]]:.3f} -> {spreads[reliable[-1]]:.3f}  => {trend}")

np.savez("out/w_dependence_v2.npz", ws=np.array(ws), lams=np.array(lams),
         S=np.array([[S[w][l] for l in lams] for w in ws]),
         PU=np.array([[PU[w][l] for l in lams] for w in ws]),
         spreads=np.array([spreads[w] for w in ws]),
         kappa=np.array([INFO[w][0] for w in ws]),
         TH=np.array([INFO[w][1] for w in ws]))

fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
cmap = plt.cm.plasma(np.linspace(0.1, 0.8, len(ws)))
for col, w in zip(cmap, ws):
    base = S[w][1.0]
    ax[0].plot(lams, [S[w][l]/base for l in lams], 'o-', color=col,
               label=f"w={w} (T_H={INFO[w][1]:.3f})")
ax[0].axhline(1, color='k', lw=0.5, alpha=0.4); ax[0].set_xscale('log', base=2)
ax[0].set_xlabel("lambda"); ax[0].set_ylabel("fuerza lengua / lam=1")
ax[0].set_title("Sensibilidad a lambda por w (banda escalada a T_H)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
ax[1].plot(ws, [spreads[w] for w in ws], 'o-', color='crimson')
ax[1].set_xlabel("w"); ax[1].set_ylabel("spread = max|r-1|")
ax[1].set_title("Sensibilidad UV vs w  (banda ~T_H, pseudo-unit controlada)")
ax[1].grid(alpha=0.3)
fig.savefig("out/w_dependence_v2.png", dpi=120)
print("\nFig -> out/w_dependence_v2.png   Datos -> out/w_dependence_v2.npz")
