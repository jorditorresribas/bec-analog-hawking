"""
FASE 3A — Canales asintoticos de Bogoliubov (conteo de canales).
SOLO dispersion y clasificacion. NO matriz S todavia.

Relacion de dispersion (hbar=m=1) en cada region asintotica homogenea:
    (omega - v k)^2 = c_s^2 k^2 + k^4/4
  <=> k^4/4 + (c_s^2 - v^2) k^2 + 2 omega v k - omega^2 = 0   (cuartico en k)

Para cada raiz real (canal propagante):
    frecuencia comovil   Omega = omega - v k
    signo de la NORMA    sign(N) = sign(Omega)
    velocidad de grupo   v_g = v + (c_s^2 k + k^3/2)/Omega
    direccion de propagacion = sign(v_g)   (NO sign(k))
Raices complejas -> evanescentes (no son canales asintoticos).

El mecanismo de Hawking: en la region SUPERSONICA aparecen canales de NORMA
NEGATIVA que no existen en la subsonica. Esa asimetria es la radiacion de Hawking.
"""
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------- Parametros asintoticos del waterfall ----------------
d = np.load("out/waterfall_state.npz")
x, n0, v0 = d["x"], d["n"], d["v"]
g, j = float(d["g"]), float(d["j"])
# region upstream (subsonica, x->-inf) y downstream (supersonica, x->+inf)
nL, vL = n0[:50].mean(), v0[:50].mean();  cL = np.sqrt(g*nL)
nR, vR = n0[-50:].mean(), v0[-50:].mean(); cR = np.sqrt(g*nR)
print("="*68)
print("FASE 3A — canales asintoticos de Bogoliubov")
print("="*68)
print(f"  Upstream  (subsonico) : n={nL:.4f} v={vL:.4f} c_s={cL:.4f} Mach={vL/cL:.3f}")
print(f"  Downstream(supersonic): n={nR:.4f} v={vR:.4f} c_s={cR:.4f} Mach={vR/cR:.3f}")


def roots_k(omega, v, c):
    """Raices del cuartico k^4/4 + (c^2-v^2)k^2 + 2 w v k - w^2 = 0."""
    return np.roots([0.25, 0.0, c**2 - v**2, 2*omega*v, -omega**2])


def classify(omega, v, c, imtol=1e-7):
    """Devuelve lista de dicts por raiz."""
    out = []
    for k in roots_k(omega, v, c):
        prop = abs(k.imag) < imtol*(1+abs(k.real))
        kr = k.real
        Omega = omega - v*kr
        if prop and abs(Omega) > 1e-12:
            vg = v + (c**2*kr + kr**3/2)/Omega
            out.append(dict(k=kr, prop=True, Omega=Omega,
                            norm=np.sign(Omega), vg=vg, dir=np.sign(vg),
                            kc=k))
        else:
            out.append(dict(k=k, prop=False, Omega=np.nan, norm=0,
                            vg=np.nan, dir=0, kc=k))
    return out


def channel_table(omega):
    print(f"\n--- omega = {omega:.3f} ---")
    for name, v, c in (("UPSTREAM  (sub)", vL, cL), ("DOWNSTREAM(sup)", vR, cR)):
        rs = classify(omega, v, c)
        prop = [r for r in rs if r['prop']]
        pos = [r for r in prop if r['norm'] > 0]
        neg = [r for r in prop if r['norm'] < 0]
        print(f"  {name}: {len(prop)} propagantes "
              f"({len(pos)} norma+, {len(neg)} norma-), "
              f"{len(rs)-len(prop)} evanescentes")
        for r in prop:
            io = "IN " if ((v < 0 or name[0]=='U') and False) else ""
            # in/out: upstream(x<0) -> vg>0 va al horizonte=IN; downstream(x>0) -> vg>0 se aleja=OUT
            if name[0] == 'U':
                ch = "IN " if r['vg'] > 0 else "OUT"
            else:
                ch = "OUT" if r['vg'] > 0 else "IN "
            sgn = '+' if r['norm'] > 0 else '-'
            print(f"      k={r['k']:+.4f}  vg={r['vg']:+.4f}  "
                  f"Omega={r['Omega']:+.4f}  norma[{sgn}]  canal {ch}")


# ---------------- Tabla de canales a varias frecuencias ----------------------
for om in (0.05, 0.1, 0.2, 0.4):
    channel_table(om)

# --------- Barrido fino: numero de canales de norma +/- vs omega -------------
omegas = np.linspace(1e-3, 1.0, 600)
nposL, nnegL, nposR, nnegR = [], [], [], []
for om in omegas:
    rL = [r for r in classify(om, vL, cL) if r['prop']]
    rR = [r for r in classify(om, vR, cR) if r['prop']]
    nposL.append(sum(r['norm']>0 for r in rL)); nnegL.append(sum(r['norm']<0 for r in rL))
    nposR.append(sum(r['norm']>0 for r in rR)); nnegR.append(sum(r['norm']<0 for r in rR))
nposL,nnegL,nposR,nnegR = map(np.array,(nposL,nnegL,nposR,nnegR))

# Umbral de Hawking: omega max donde la region supersonica tiene canal norma-
has_neg = nnegR > 0
om_max = omegas[has_neg].max() if has_neg.any() else 0.0
print("\n" + "="*68)
print(f"  VENTANA DE HAWKING: la region supersonica tiene canal de NORMA NEGATIVA")
print(f"  para 0 < omega < omega_max = {om_max:.4f}")
print(f"  (la subsonica nunca tiene norma negativa: {nnegL.max()} canales neg.)")
print("  -> asimetria de norma al cruzar el horizonte = base de la emision Hawking")
print("="*68)

# -------------------- VALIDACION: residuos y v_g numerica --------------------
print("\nVALIDACION:")
om = 0.1
maxres = 0.0; maxvgerr = 0.0
for v, c in ((vL, cL), (vR, cR)):
    for r in classify(om, v, c):
        k = r['kc']
        res = abs(0.25*k**4 + (c**2-v**2)*k**2 + 2*om*v*k - om**2)
        maxres = max(maxres, res)
        if r['prop']:
            # v_g numerica: dω/dk con ω(k)=vk+sign(Ω)sqrt(c^2k^2+k^4/4)
            kk = r['k']; s = np.sign(r['Omega']); h = 1e-6
            wf = lambda q: v*q + s*np.sqrt(c**2*q**2 + q**4/4)
            vg_num = (wf(kk+h)-wf(kk-h))/(2*h)
            maxvgerr = max(maxvgerr, abs(vg_num - r['vg']))
print(f"  max residuo del cuartico en las raices   = {maxres:.2e}")
print(f"  max |v_g analitica - v_g numerica|        = {maxvgerr:.2e}")

# ------------------------------- Figuras -------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
kk = np.linspace(-3, 3, 800)
for a, (name, v, c) in zip(ax[:2], (("UPSTREAM subsonico", vL, cL),
                                     ("DOWNSTREAM supersonico", vR, cR))):
    for s in (+1, -1):
        W = v*kk + s*np.sqrt(c**2*kk**2 + kk**4/4)
        a.plot(kk, W, color='navy' if s>0 else 'crimson',
               label=r'$\Omega>0$ (norma +)' if s>0 else r'$\Omega<0$ (norma -)')
    for om in (0.1, 0.2):
        a.axhline(om, color='gray', ls=':', lw=0.8)
    a.axhline(0, color='k', lw=0.5); a.set_title(name)
    a.set_xlabel('k'); a.set_ylabel(r'$\omega$'); a.set_ylim(-1, 1); a.legend(fontsize=8)
    a.grid(alpha=0.3)

ax[2].plot(omegas, nposR, 'b-', label='norma + (downstream)')
ax[2].plot(omegas, nnegR, 'r-', lw=2, label='norma - (downstream) = Hawking')
ax[2].plot(omegas, nposL, 'b--', label='norma + (upstream)')
ax[2].plot(omegas, nnegL, 'r--', label='norma - (upstream)')
ax[2].axvline(om_max, color='k', ls=':', label=f'$\\omega_{{max}}$={om_max:.3f}')
ax[2].set_xlabel(r'$\omega$'); ax[2].set_ylabel('nº canales propagantes')
ax[2].set_title('Conteo de canales vs frecuencia'); ax[2].legend(fontsize=8)
ax[2].grid(alpha=0.3)
plt.tight_layout(); plt.savefig("out/bdg_channels.png", dpi=120); plt.close()
print("\nFigura: out/bdg_channels.png")
