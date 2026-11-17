"""
Experimento RUMBO #1 (b): ROBUSTEZ de la lengua de Hawking frente al termino k^4
(primer paso trans-planckiano). Variamos lam en (omega-vk)^2 = c^2 k^2 + lam k^4/4
con c,v (metrica, kappa) FIJOS. op='lambda' para TODOS (lam=1 reproduce validado).

Para cada lam: banda omega adaptada a omega_max(lam); correlacion d2 (Hawking puro)
y total; medimos fuerza de la anti-correlacion cruzada (lengua) y la pendiente.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import correlations as C

L, N, w = 30.0, 3000, 2.0
ncoarse = 360
lams = [0.5, 1.0, 2.0]
idx = np.linspace(0, N-1, ncoarse).astype(int)
xc = np.linspace(-L, L, N)[idx]

def omega_max(lam):
    lo, hi = 0.05, 1.3
    for _ in range(26):
        mid = 0.5*(lo+hi)
        r = C.scatter_field(L=24, N=2000, omega=mid, channel='d2', lam=lam,
                            mode='horizon', op='lambda')
        if r['Rin_neg'] is not None and r['nchan'] == 2:
            lo = mid
        else:
            hi = mid
    return lo

def probe(Cmat, slope, xs):
    return np.array([Cmat[np.argmin(np.abs(xc-x)), np.argmin(np.abs(xc-slope*x))]
                     for x in xs])

res = {}
for lam in lams:
    omax = omega_max(lam)
    omegas = np.linspace(0.02, 0.9*omax, 90)
    print(f"lam={lam}: omega_max={omax:.3f}, banda 0.02..{0.9*omax:.3f} ({len(omegas)} pts)")
    xc_, Cd2, nc = C.density_correlation(L, N, omegas, w=w, lam=lam, mode='horizon',
                                         ncoarse=ncoarse, terms='d2', op='lambda')
    _, Call, _ = C.density_correlation(L, N, omegas, w=w, lam=lam, mode='horizon',
                                       ncoarse=ncoarse, terms='all', op='lambda')
    # fuerza de la lengua: integral de la parte NEGATIVA de d2 en cuadrante cruzado
    iu = (xc < -0.5) & (xc > -14); idn = (xc > 0.5) & (xc < 14)
    blk = Cd2[np.ix_(iu, idn)]
    neg_int = blk[blk < 0].sum()*(xc[1]-xc[0])**2
    xs = np.linspace(-13, -0.5, 26)
    pr = probe(Cd2, -1.26, xs)
    res[lam] = dict(Cd2=Cd2, Call=Call, omax=omax, neg_int=neg_int,
                    tongue_min=blk.min(), probe_mean=pr.mean(), probe_min=pr.min())
    print(f"   lengua d2: integral_neg={neg_int:.3e}  min={blk.min():.3e}  "
          f"recta-1.26 mean={pr.mean():+.3e}")

np.savez("out/lambda_sweep.npz", xc=xc, lams=lams,
         **{f"Cd2_{l}": res[l]['Cd2'] for l in lams},
         **{f"Call_{l}": res[l]['Call'] for l in lams})

print("\n--- ROBUSTEZ de la lengua de Hawking vs lam ---")
print(" lam   omega_max   integral_neg(d2)   min(d2)    mean recta-1.26")
for l in lams:
    r = res[l]
    print(f" {l:4}   {r['omax']:.3f}      {r['neg_int']:+.3e}     {r['tongue_min']:+.3e}  "
          f"{r['probe_mean']:+.3e}")
base = res[1.0]['neg_int']
print("\n normalizado a lam=1:")
for l in lams:
    print(f"   lam={l}: fuerza lengua / lam1 = {res[l]['neg_int']/base:.2f}")

# figura: d2 (Hawking) para cada lam
win15 = np.abs(xc) < 15
vmax = np.percentile(np.abs(res[1.0]['Cd2'][np.ix_(win15, win15)]), 99)
fig, ax = plt.subplots(1, len(lams), figsize=(5.4*len(lams), 5.4), constrained_layout=True)
for a, lam in zip(ax, lams):
    im = a.pcolormesh(xc, xc, res[lam]['Cd2'], cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      shading='auto', rasterized=True)
    xx = np.linspace(-15, 0, 10)
    a.plot(xx, -1.26*xx, 'g--', lw=1.0, alpha=0.7); a.plot(-1.26*xx, xx, 'g--', lw=1.0, alpha=0.7)
    a.axhline(0, color='k', lw=0.4, alpha=0.3); a.axvline(0, color='k', lw=0.4, alpha=0.3)
    a.set_xlim(-15, 15); a.set_ylim(-15, 15); a.set_aspect('equal')
    a.set_xlabel("x"); a.set_ylabel("x'")
    a.set_title(f"lambda={lam}  (omega_max={res[lam]['omax']:.2f})")
    fig.colorbar(im, ax=a, shrink=0.75)
fig.suptitle("Robustez de la lengua de Hawking (canal d2) frente al termino k^4 "
             "[metrica c,v FIJA]", fontsize=13)
fig.savefig("out/lambda_sweep.png", dpi=120)
print("\nFig -> out/lambda_sweep.png   Datos -> out/lambda_sweep.npz")
