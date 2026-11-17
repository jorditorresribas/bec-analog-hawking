"""
RUMBO #1 (b) — PROFUNDIZAR ROBUSTEZ, paso 3: METRICA LIMPIA de la lengua
+ paso 2: EXTENDER lam a {0.25, 0.5, 1, 2, 4}.

Problema (caveat 1 de CORRELATIONS.md): la metrica "integral de la parte negativa
en el cuadrante cruzado" mezcla la cresta genuina de Hawking con las FRANJAS de
oscilacion, cuya longitud de onda CAMBIA con lam (k del compañero cambia) -> parte
del +-10% podria ser artefacto trivial del espaciado de franjas.

Metrica limpia (independiente de la fase de las franjas): tomar la correlacion
C_d2 a lo largo de la RECTA del compañero (ridge) x' = slope*x, calcular su
ENVOLVENTE de Hilbert (elimina la oscilacion) y:
  - integrar la envolvente a lo largo de la recta        -> fuerza limpia
  - su valor extrapolado al IR (s grande = separacion grande) -> señal Hawking IR
La envolvente no depende de donde caigan los nodos de las franjas -> robusta a lam.

slope = (v_d - c_d)/(v_u - c_u) (velocidad de grupo del par lento), medido -1.17,
predicho -1.26; usamos el PREDICHO geometrico (no ajustado) para no sesgar.
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

# pendiente geometrica predicha del compañero
n_u, mach_u = C.n_u, 0.5
vu = C.vLb; cu = np.sqrt(C.g*n_u)
n_d = (C.jcur/(2.0*np.sqrt(C.g)))**(2/3); vd = C.jcur/n_d; cd = np.sqrt(C.g*n_d)
slope = (vd - cd)/(vu - cu)
print(f"pendiente predicha compañero = {slope:.3f}  (v_d-c_d={vd-cd:.3f}, v_u-c_u={vu-cu:.3f})")

def omega_max(lam):
    lo, hi = 0.05, 1.3
    for _ in range(26):
        mid = 0.5*(lo+hi)
        r = C.scatter_field(L=24, N=2000, omega=mid, channel='d2', lam=lam,
                            mode='horizon', op='lambda')
        lo, hi = (mid, hi) if (r['Rin_neg'] is not None and r['nchan'] == 2) else (lo, mid)
    return lo

def ridge_profile(xc, Cmat, smax=13.0, ns=400):
    """C_d2 a lo largo de la recta del compañero: x_u=-s (upstream, x<0),
    x_d=+|slope|*s (downstream, x>0). Devuelve s, g(s), envolvente E(s)."""
    interp = RegularGridInterpolator((xc, xc), Cmat, bounds_error=False, fill_value=0.0)
    s = np.linspace(0.5, smax, ns)
    xu = -s; xd = abs(slope)*s
    g = interp(np.column_stack([xu, xd]))
    # envolvente de Hilbert (elimina oscilacion de franjas; conserva la cresta suave)
    E = np.abs(hilbert(g))
    return s, g, E

def metrics(lam, ncoarse=ncoarse, N=N, L=L, nomega=90):
    omax = omega_max(lam)
    omegas = np.linspace(0.02, 0.9*omax, nomega)
    xc, Cd2, _ = C.density_correlation(L, N, omegas, w=w, lam=lam, mode='horizon',
                                       ncoarse=ncoarse, terms='d2', op='lambda')
    s, g, E = ridge_profile(xc, Cd2)
    env_int = np.trapz(E, s)                 # fuerza limpia (envolvente integrada)
    # señal IR: media de la envolvente en la cola lejana (s grande) y ajuste exp
    tail = E[s > 8.0]; ir_tail = tail.mean()
    # metrica integral antigua (para comparar)
    dx = xc[1]-xc[0]
    iu = (xc < -0.5) & (xc > -14); idn = (xc > 0.5) & (xc < 14)
    blk = Cd2[np.ix_(iu, idn)]
    neg_int = blk[blk < 0].sum()*dx**2
    return dict(omax=omax, xc=xc, Cd2=Cd2, s=s, g=g, E=E,
                env_int=env_int, ir_tail=ir_tail, env_peak=E.max(), neg_int=neg_int)

print("\nExtendiendo lam a {0.25, 0.5, 1, 2, 4} con metrica LIMPIA (envolvente ridge)")
print(f"{'lam':>5} {'omax':>6} | {'env_int':>10} {'env_peak':>9} {'ir_tail':>9} | "
      f"{'neg_int(viejo)':>13}")
R = {}
for lam in lams:
    m = metrics(lam); R[lam] = m
    print(f"{lam:>5} {m['omax']:>6.3f} | {m['env_int']:>10.3e} {m['env_peak']:>9.2e} "
          f"{m['ir_tail']:>9.2e} | {m['neg_int']:>13.3e}")

base_e = R[1.0]['env_int']; base_o = R[1.0]['neg_int']; base_ir = R[1.0]['ir_tail']
print("\n--- ROBUSTEZ: cocientes normalizados a lam=1 ---")
print(f"{'lam':>5} | {'env_int/1':>9} {'ir_tail/1':>9} | {'neg_int/1 (viejo)':>16}")
for lam in lams:
    print(f"{lam:>5} | {R[lam]['env_int']/base_e:>9.3f} {R[lam]['ir_tail']/base_ir:>9.3f} "
          f"| {R[lam]['neg_int']/base_o:>16.3f}")

np.savez("out/robustness_clean.npz", lams=np.array(lams), slope=slope,
         **{f"s_{l}": R[l]['s'] for l in lams},
         **{f"g_{l}": R[l]['g'] for l in lams},
         **{f"E_{l}": R[l]['E'] for l in lams},
         env_int=np.array([R[l]['env_int'] for l in lams]),
         ir_tail=np.array([R[l]['ir_tail'] for l in lams]),
         neg_int=np.array([R[l]['neg_int'] for l in lams]))

# figura: perfil de ridge + envolvente para cada lam
fig, ax = plt.subplots(1, 2, figsize=(13, 5), constrained_layout=True)
cmap = plt.cm.viridis(np.linspace(0, 1, len(lams)))
for col, lam in zip(cmap, lams):
    ax[0].plot(R[lam]['s'], R[lam]['g'], color=col, lw=0.8, alpha=0.5)
    ax[0].plot(R[lam]['s'], R[lam]['E'], color=col, lw=1.8, label=f"lam={lam}")
ax[0].set_xlabel("s  (separacion a lo largo de la recta del compañero)")
ax[0].set_ylabel("C_d2 en la recta  (fino) + envolvente (grueso)")
ax[0].set_title("Ridge de Hawking y su envolvente de Hilbert")
ax[0].legend(fontsize=8); ax[0].axhline(0, color='k', lw=0.4, alpha=0.3)
# robustez: cocientes vs lam
le = [R[l]['env_int']/base_e for l in lams]
li = [R[l]['ir_tail']/base_ir for l in lams]
lo = [R[l]['neg_int']/base_o for l in lams]
ax[1].plot(lams, le, 'o-', label='envolvente integrada (limpia)')
ax[1].plot(lams, li, 's-', label='cola IR (envolvente, s>8)')
ax[1].plot(lams, lo, '^--', label='integral neg (vieja, mezcla franjas)', alpha=0.6)
ax[1].axhline(1, color='k', lw=0.5, alpha=0.4); ax[1].set_xscale('log', base=2)
ax[1].set_xlabel("lambda (coeficiente dispersivo k^4)")
ax[1].set_ylabel("fuerza de la lengua / lam=1")
ax[1].set_title("Robustez trans-planckiana (metrica limpia vs vieja)")
ax[1].legend(fontsize=8); ax[1].grid(alpha=0.3)
fig.suptitle("Robustez de la lengua de Hawking: metrica limpia (envolvente) sobre "
             "lam x16 (0.25->4)", fontsize=12)
fig.savefig("out/robustness_clean.png", dpi=120)
print("\nFig -> out/robustness_clean.png   Datos -> out/robustness_clean.npz")
