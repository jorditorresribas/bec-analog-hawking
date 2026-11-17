"""
Descomposicion de la correlacion de Hawking (Ec.32) en:
  - d2  : canal ANOMALO (norma negativa) = produccion de pares = lengua Hawking PURA
  - ud1 : canales ordinarios (norma +)
  - all : suma (Ec.32 completa)
Una sola pasada (reusa cada scatter_field una vez). taper Hann en bordes de banda.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt
import correlations as C

L, N, w = 30.0, 3600, 2.0
omegas = np.linspace(0.02, 0.45, 140)
ncoarse = 420
idx = np.linspace(0, N-1, ncoarse).astype(int)
xfull = np.linspace(-L, L, N); xc = xfull[idx]
nc = C.bg(xfull, w, 'horizon')[0][idx]
dω = omegas[1]-omegas[0]
nn = len(omegas)
win = np.clip(0.5*(1-np.cos(np.pi*np.minimum(np.arange(nn), nn-1-np.arange(nn))
                            /(0.15*nn))), 0, 1)

acc = {'d2': np.zeros((ncoarse, ncoarse), complex),
       'ud1': np.zeros((ncoarse, ncoarse), complex)}
print(f"Descomposicion, {nn} omega-pts, N={N}...")
for iw, om in enumerate(omegas):
    M = {}
    ok = True
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L, N, om, ch, w, 1.0, 'horizon')
        if r['Rin_neg'] is None or r['nchan'] != 2:
            ok = False; break
        M[ch] = dict(u=r['u'][idx], v=r['v_mode'][idx], r=r['u'][idx]+r['v_mode'][idx])
    if not ok:
        continue
    bud = np.outer(np.conj(M['u']['v']), M['u']['r']) + \
          np.outer(np.conj(M['d1']['v']), M['d1']['r'])
    bd2 = np.outer(np.conj(M['d2']['u']), M['d2']['r'])
    acc['ud1'] += (bud+np.conj(bud.T))*dω*win[iw]
    acc['d2'] += (bd2+np.conj(bd2.T))*dω*win[iw]

pref = np.sqrt(np.outer(nc, nc))/(4*np.pi)
Cd2 = np.real(pref*acc['d2'])
Cud1 = np.real(pref*acc['ud1'])
Call = Cd2+Cud1
np.savez("out/correlation_decomp.npz", xc=xc, Cd2=Cd2, Cud1=Cud1, Call=Call, nc=nc)

# --- sondeo de la lengua (solo d2) por la recta del companero lento ---
def probe(Cmat, slope, xs):
    return np.array([Cmat[np.argmin(np.abs(xc-x)), np.argmin(np.abs(xc-slope*x))]
                     for x in xs])
xs = np.linspace(-13, -0.5, 26)
for nm, Cm in [('d2 (Hawking puro)', Cd2), ('all', Call)]:
    p = probe(Cm, -1.26, xs)
    print(f"  {nm:18s} recta -1.26x: mean={p.mean():+.3e} min={p.min():+.3e} "
          f"frac_neg={np.mean(p<0):.2f}")

# --- figura 3 paneles ---
win15 = np.abs(xc) < 15
vmax = np.percentile(np.abs(Call[np.ix_(win15, win15)]), 99)
fig, ax = plt.subplots(1, 3, figsize=(17, 5.6), constrained_layout=True)
for a, (Cm, ti) in zip(ax, [(Call, 'TOTAL (Ec.32)'),
                            (Cd2, 'd2: canal anomalo = Hawking'),
                            (Cud1, 'u,d1: ordinarios')]):
    im = a.pcolormesh(xc, xc, Cm, cmap='RdBu_r', vmin=-vmax, vmax=vmax,
                      shading='auto', rasterized=True)
    xx = np.linspace(-15, 0, 10)
    a.plot(xx, -1.26*xx, 'g--', lw=1.0, alpha=0.7)
    a.plot(-1.26*xx, xx, 'g--', lw=1.0, alpha=0.7)
    a.axhline(0, color='k', lw=0.4, alpha=0.3); a.axvline(0, color='k', lw=0.4, alpha=0.3)
    a.set_xlim(-15, 15); a.set_ylim(-15, 15); a.set_aspect('equal')
    a.set_xlabel("x"); a.set_ylabel("x'"); a.set_title(ti)
    fig.colorbar(im, ax=a, shrink=0.75)
fig.suptitle("Correlacion Hawking: la lengua vive en el canal anomalo d2 "
             "(norma negativa)", fontsize=13)
fig.savefig("out/correlation_decomp.png", dpi=125)
print("Fig -> out/correlation_decomp.png")
