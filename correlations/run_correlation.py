"""
Experimento RUMBO #1 (a): la LENGUA de correlacion de Hawking (refinado).
C(x,x') = <delta n(x) delta n(x')> en vacio-in, Ec.32 Recati-Pavloff-Carusotto.

Senal esperada: banda NEGATIVA (anti-correlacion Hawking-companero) en el cuadrante
upstream(x<0) x downstream(x'>0) a lo largo de la recta x' = slope * x con
slope_lento = (v_d-c_d)/(v_u-c_u) = -1.26 (companero lento) y slope_rapido = -3.78.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import correlations as C

L, N, w = 30.0, 3600, 2.0
omegas = np.linspace(0.02, 0.45, 140)
print(f"Banda omega: {omegas[0]:.3f}..{omegas[-1]:.3f} ({len(omegas)} pts)  L={L} N={N} w={w}")

print("HORIZONTE, integrando correlacion...")
xc, Ch, nc = C.density_correlation(L, N, omegas, w=w, lam=1.0, mode='horizon',
                                   ncoarse=420)
np.savez("out/correlation_lambda1.npz", xc=xc, Chor=Ch, nc=nc, omegas=omegas,
         L=L, N=N, w=w)

# --- chequeo de convergencia en omega (mitad de puntos) ---
xc2, Ch2, _ = C.density_correlation(L, N, omegas[::2], w=w, lam=1.0, mode='horizon',
                                    ncoarse=420)
# comparar en ventana interior |x|,|x'|<15
win = (np.abs(xc) < 15)
rel = np.abs(Ch[np.ix_(win, win)]-Ch2[np.ix_(win, win)]).max() / np.abs(Ch[np.ix_(win, win)]).max()
print(f"convergencia omega (140 vs 70 pts), ventana |x|<15: dif rel max = {rel:.2e}")

# --- sondeo a lo largo de las rectas predichas ---
def probe(slope, xs):
    out = []
    for x in xs:
        xp = slope*x
        i = np.argmin(np.abs(xc-x)); j = np.argmin(np.abs(xc-xp))
        out.append(Ch[i, j])
    return np.array(out)

xs = np.linspace(-14, -0.5, 28)      # x upstream (<0) -> xp = slope*x >0 downstream
for slope, name in [(-1.26, 'companero lento (v-c)'), (-3.78, 'companero rapido (v+c)')]:
    p = probe(slope, xs)
    print(f"\nrecta x'={slope}x [{name}]: C medio={p.mean():+.3e}  min={p.min():+.3e}  "
          f"frac negativa={np.mean(p<0):.2f}")

# integral de la banda negativa cruzada (fuerza total de la lengua)
iu = (xc < -0.5) & (xc > -15); idn = (xc > 0.5) & (xc < 15)
blk = Ch[np.ix_(iu, idn)]
print(f"\nBloque cruzado upstream x downstream (|x|<15):")
print(f"  min={blk.min():.3e}  max={blk.max():.3e}  mean={blk.mean():+.3e}")
print(f"  integral de la parte negativa = {blk[blk<0].sum()*(xc[1]-xc[0])**2:.3e}")

# --- figura: mapa + zoom con rectas predichas ---
vmax = np.percentile(np.abs(Ch[np.ix_(win, win)]), 99)
fig, ax = plt.subplots(1, 2, figsize=(13.5, 6), constrained_layout=True)
for a, (xlim) in zip(ax, [38, 15]):
    im = a.pcolormesh(xc, xc, Ch, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      shading='auto', rasterized=True)
    a.axhline(0, color='k', lw=0.5, alpha=0.3); a.axvline(0, color='k', lw=0.5, alpha=0.3)
    xx = np.linspace(-xlim, 0, 10)
    a.plot(xx, -1.26*xx, 'g--', lw=1.1, alpha=0.8, label="x'=-1.26x (companero lento)")
    a.plot(xx, -3.78*xx, 'm--', lw=1.1, alpha=0.8, label="x'=-3.78x (rapido)")
    a.set_xlabel("x"); a.set_ylabel("x'"); a.set_aspect('equal')
    a.set_xlim(-xlim, xlim); a.set_ylim(-xlim, xlim)
    fig.colorbar(im, ax=a, shrink=0.8, label=r"$\langle\delta n(x)\delta n(x')\rangle$")
ax[0].set_title("Mapa completo"); ax[1].set_title("Zoom al horizonte + rectas predichas")
ax[1].legend(loc='upper right', fontsize=8)
fig.suptitle("Lengua de correlacion de Hawking (BEC acoustic BH, vacio-in T=0, lambda=1)",
             fontsize=13)
fig.savefig("out/correlation_tongue.png", dpi=130)
print("\nFig -> out/correlation_tongue.png   Datos -> out/correlation_lambda1.npz")
