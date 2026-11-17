"""
RUMBO #1 (b) — OPCION A, version de CIERRE (satisface el criterio de ChatGPT).

Dependencia con w de la sensibilidad a lambda de la lengua de Hawking, con:
  (1) banda termica comun omega in [a,b]*T_H(w), MISMA regla de seleccion para todo
      (w,lambda); se verifica que NO hay truncado diferencial contra omega_max.
  (2) pseudo-unitaridad POR FRECUENCIA eps(omega)=max_canal|flux_out-eta_in|, de los
      MISMOS solves que construyen C (coherencia total); se cuenta n_valid y se mira
      si los puntos descartados se concentran en algun extremo de lambda.
  (3) spread(w) = max_lambda R_lambda - min_lambda R_lambda (def. de ChatGPT),
      publicando los TRES cocientes R_lambda (no solo el spread).
  (4) tabla completa: w, kappa, T_H, banda, n_valid/n_tot, eps mediana/max, R_lambda.

Se cambia SOLO w (y L,N que w exige). Observable y metrica identicos a §7.
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
A_TH, B_TH = 1.0, 4.0        # banda termica comun (unidades de T_H)
NOMEGA = 60
PU_TOL = 1e-2                 # umbral de frecuencia fiable

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
    ih = np.argmin(np.abs(cs - v)); k = abs(df[ih])
    return k, k/(2*np.pi)

def LN_for(w):
    L = max(40.0, 8*w + 16); N = int(round(2*L/0.02))
    return L, N

def correlation_with_pu(L, N, omegas, w, lam, ncoarse=360):
    """Construye C_d2 Y registra eps(omega) de los MISMOS solves (3 canales/freq).
    Devuelve xc, C_d2, y arrays por frecuencia: eps, valido(bool)."""
    xfull = np.linspace(-L, L, N)
    idx = np.linspace(0, N-1, ncoarse).astype(int); xc = xfull[idx]
    nc = C.bg(xfull, w, 'horizon')[0][idx]
    acc = np.zeros((ncoarse, ncoarse), complex)
    dω = omegas[1]-omegas[0]
    eps = np.full(len(omegas), np.nan); valid = np.zeros(len(omegas), bool)
    win = 0.5*(1-np.cos(np.pi*np.minimum(np.arange(len(omegas)),
               len(omegas)-1-np.arange(len(omegas)))/(0.15*len(omegas))))
    win = np.clip(win, 0, 1)
    for iw, om in enumerate(omegas):
        modes = {}; e = 0.0; ok = True
        for ch in ('u', 'd1', 'd2'):
            r = C.scatter_field(L, N, om, ch, w, lam, 'horizon', 'lambda')
            if r['Rin_neg'] is None or r['nchan'] != 2:
                ok = False; break
            e = max(e, abs(r['flux_out'] - r['eta_in']))
            modes[ch] = dict(u=r['u'][idx], r=(r['u']+r['v_mode'])[idx])
        if not ok:
            continue
        eps[iw] = e; valid[iw] = (e < PU_TOL)
        br = np.outer(np.conj(modes['d2']['u']), modes['d2']['r'])
        br = br + np.conj(br.T)
        acc += br*dω*win[iw]
    pref = np.sqrt(np.outer(nc, nc))/(4*np.pi)
    return xc, np.real(pref*acc), eps, valid

def measure(w, lam):
    L, N = LN_for(w); _, TH = kappa_TH(w)
    olo, ohi_req = A_TH*TH, B_TH*TH
    ohi_cap = 0.9*OMAX[lam]
    truncated = ohi_req > ohi_cap
    ohi = min(ohi_req, ohi_cap)
    omegas = np.linspace(olo, ohi, NOMEGA)
    xc, Cd2, eps, valid = correlation_with_pu(L, N, omegas, w, lam)
    dx = xc[1]-xc[0]; xhi = 14.0 + 2.0*w
    iu = (xc < -0.5) & (xc > -xhi); idn = (xc > 0.5) & (xc < xhi)
    blk = Cd2[np.ix_(iu, idn)]
    return dict(str=blk[blk < 0].sum()*dx**2, L=L, N=N, TH=TH, band=(olo, ohi),
                truncated=truncated, eps=eps, valid=valid,
                n_valid=int(valid.sum()), n_tot=len(omegas),
                eps_med=np.nanmedian(eps), eps_max=np.nanmax(eps))

OMAX = {lam: omega_max(lam) for lam in lams}
print("Banda termica comun: omega in [%.1f,%.1f]*T_H(w).  PU_TOL=%.0e" % (A_TH, B_TH, PU_TOL))
print("omega_max(lam):  " + "  ".join(f"l={l}:{OMAX[l]:.3f}" for l in lams))

M = {}
for w in ws:
    M[w] = {lam: measure(w, lam) for lam in lams}

print("\n" + "="*86)
print("TABLA DE CIERRE (opcion A)")
print("="*86)
hdr = f"{'w':>3} {'kappa':>6} {'T_H':>6} {'L':>4} {'N':>5} {'lam':>4} " \
      f"{'banda':>15} {'trunc':>5} {'str':>11} {'R_lam':>6} {'n_val':>6} " \
      f"{'eps_med':>8} {'eps_max':>8}"
print(hdr)
spreads = {}; Rtab = {}
for w in ws:
    base = M[w][1.0]['str']; R = {}
    for lam in lams:
        m = M[w][lam]; R[lam] = m['str']/base
        k = kappa_TH(w)[0]
        print(f"{w:>3} {k:>6.3f} {m['TH']:>6.4f} {m['L']:>4.0f} {m['N']:>5} {lam:>4} "
              f"[{m['band'][0]:.3f},{m['band'][1]:.3f}] {str(m['truncated']):>5} "
              f"{m['str']:>+11.3e} {R[lam]:>6.3f} {m['n_valid']:>3}/{m['n_tot']:<2} "
              f"{m['eps_med']:>8.1e} {m['eps_max']:>8.1e}")
    Rtab[w] = R
    spreads[w] = max(R.values()) - min(R.values())    # def. ChatGPT: max-min
    print("-"*86)

print("\nRESUMEN  (spread = max_lam R_lam - min_lam R_lam):")
print(f"{'w':>4} | {'R(0.5)':>7} {'R(1)':>6} {'R(2)':>6} | {'spread':>7} | "
      f"{'eps_max banda':>13} | fiable?")
for w in ws:
    emax = max(M[w][l]['eps_max'] for l in lams)
    rel = all(M[w][l]['eps_max'] < PU_TOL for l in lams) and \
          all(M[w][l]['n_valid'] == M[w][l]['n_tot'] for l in lams)
    print(f"{w:>4} | {Rtab[w][0.5]:>7.3f} {Rtab[w][1.0]:>6.3f} {Rtab[w][2.0]:>6.3f} | "
          f"{spreads[w]:>7.3f} | {emax:>13.1e} | {'SI' if rel else 'NO'}")

# puntos descartados: se concentran en algun extremo de lambda?
print("\nFrecuencias descartadas (eps>%.0e o canal ausente) por (w,lam):" % PU_TOL)
for w in ws:
    drops = {l: M[w][l]['n_tot']-M[w][l]['n_valid'] for l in lams}
    print(f"  w={w}: " + "  ".join(f"lam={l}:{drops[l]}" for l in lams))

print("\n" + "="*86)
rel_ws = [w for w in ws if all(M[w][l]['eps_max'] < PU_TOL and
          M[w][l]['n_valid'] == M[w][l]['n_tot'] for l in lams)]
print(f"w con controles limpios: {rel_ws}")
if len(rel_ws) >= 2:
    s0, s1 = spreads[rel_ws[0]], spreads[rel_ws[-1]]
    verdict = ("DISMINUYE con w -> compatible con menor sensibilidad UV" if s1 < s0*0.8
               else "AUMENTA con w -> hipotesis simple NO respaldada" if s1 > s0*1.2
               else "aprox CONSTANTE -> sin mejora clara de robustez con w")
    print(f"spread: w={rel_ws[0]} ({s0:.3f}) -> w={rel_ws[-1]} ({s1:.3f})  =>  {verdict}")
else:
    print("Menos de 2 w con controles limpios -> NO concluyente (contaminacion numerica).")

np.savez("out/w_dependence_v3.npz", ws=np.array(ws), lams=np.array(lams),
         R=np.array([[Rtab[w][l] for l in lams] for w in ws]),
         spreads=np.array([spreads[w] for w in ws]),
         kappa=np.array([kappa_TH(w)[0] for w in ws]),
         TH=np.array([kappa_TH(w)[1] for w in ws]),
         eps_max=np.array([[M[w][l]['eps_max'] for l in lams] for w in ws]),
         n_valid=np.array([[M[w][l]['n_valid'] for l in lams] for w in ws]))

fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
cmap = plt.cm.plasma(np.linspace(0.1, 0.8, len(ws)))
for col, w in zip(cmap, ws):
    ax[0].plot(lams, [Rtab[w][l] for l in lams], 'o-', color=col,
               label=f"w={w} (T_H={kappa_TH(w)[1]:.3f})")
ax[0].axhline(1, color='k', lw=0.5, alpha=0.4); ax[0].set_xscale('log', base=2)
ax[0].set_xlabel("lambda"); ax[0].set_ylabel("R_lam = str/str(lam=1)")
ax[0].set_title("Sensibilidad a lambda por w (banda ~T_H, PU por freq)")
ax[0].legend(fontsize=8); ax[0].grid(alpha=0.3)
ax[1].plot(ws, [spreads[w] for w in ws], 'o-', color='crimson')
ax[1].set_xlabel("w"); ax[1].set_ylabel("spread = max R - min R")
ax[1].set_title("Sensibilidad UV vs w")
ax[1].grid(alpha=0.3)
fig.savefig("out/w_dependence_v3.png", dpi=120)
print("\nFig -> out/w_dependence_v3.png   Datos -> out/w_dependence_v3.npz")
