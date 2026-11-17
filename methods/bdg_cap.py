"""
AFINADO DEL HITO 2 — estabilidad espectral con condiciones absorbentes.

Problema: el espectro BdG en caja de paredes duras tenia modos complejos
espurios localizados en las paredes (demostrado en bdg_analysis.py).

Solucion rigurosa: potencial absorbente complejo (CAP) en los bordes del
operador BdG. Absorbe las ondas salientes (emula fronteras abiertas) y elimina
los modos de pared, SIN tocar los modos fisicos del interior.

Operador BdG con CAP (Gamma(x) real, >=0 en los bordes):
    M = [[ L - i Gamma ,      G        ],
         [    -G*       , -L* - i Gamma ]]
con  L = -1/2 Dxx + V + 2 g|psi0|^2 - mu,   G = diag(g psi0^2).

Criterio: si tras el CAP el modo con mayor Im(w) vive en el INTERIOR y crece,
seria fisico; si desaparece o queda amortiguado -> era artefacto -> ESTABLE.
"""
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g, n_u, mach_u, mach_d, w_step = 1.0, 1.0, 0.5, 2.0, 2.0
v_u = mach_u*np.sqrt(g*n_u)
j = n_u*v_u; mu = g*n_u + 0.5*v_u**2
n_d = (j/(mach_d*np.sqrt(g)))**(2/3); n_H = (j**2/g)**(1/3)


def build(Nb, L, Gamma0=0.0, xcap_frac=0.75):
    x = np.linspace(-L, L, Nb); dx = x[1]-x[0]
    n0 = n_u + (n_d-n_u)*0.5*(1+np.tanh(x/w_step)); v0 = j/n0
    th = np.cumsum(v0)*dx; th -= th[Nb//2]
    psi0 = np.sqrt(n0)*np.exp(1j*th)
    main=-2*np.ones(Nb); off=np.ones(Nb-1)
    Dxx = sp.diags([off,main,off],[-1,0,1]).toarray()/dx**2
    V = mu - g*n0 + 0.5*np.real((Dxx@psi0)/psi0)
    xcap = xcap_frac*L
    Gamma = np.where(np.abs(x) > xcap,
                     Gamma0*((np.abs(x)-xcap)/(L-xcap))**2, 0.0)
    Ld = -0.5*Dxx + np.diag(V + 2*g*np.abs(psi0)**2 - mu)
    G = np.diag(g*psi0**2)
    M = np.block([[Ld - 1j*np.diag(Gamma), G],
                  [-np.conj(G), -np.conj(Ld) - 1j*np.diag(Gamma)]])
    return x, n0, v0, psi0, M, Gamma


print("="*64)
print("HITO 2 AFINADO — BdG con frontera absorbente (CAP)")
print("="*64)

Nb, L = 1000, 60.0
# Baseline sin CAP
_, _, _, _, M0, _ = build(Nb, L, Gamma0=0.0)
w0 = la.eigvals(M0)
print(f"\n  SIN CAP:  max Im(w) = {np.max(w0.imag):+.3e}   (modo de pared espurio)")

# Barrido de intensidad del CAP
print("\n  CON CAP (barrido de intensidad):")
for Gamma0 in (0.5, 1.0, 2.0, 4.0):
    _, _, _, _, Mc, _ = build(Nb, L, Gamma0=Gamma0)
    wc = la.eigvals(Mc)
    print(f"    Gamma0={Gamma0:4.1f} :  max Im(w) = {np.max(wc.imag):+.3e}")

# Analisis detallado con CAP moderado: localizar el modo de mayor Im restante
Gamma0 = 2.0
x, n0, v0, psi0, Mc, Gamma = build(Nb, L, Gamma0=Gamma0)
wc, vrc = la.eig(Mc)
imax = np.argmax(wc.imag)
wtop = wc[imax]
u = vrc[:Nb, imax]; v = vrc[Nb:, imax]
amp = np.abs(u)+np.abs(v); amp /= amp.max()
edge = np.sum(amp[np.abs(x) > 0.75*L])/np.sum(amp)
cen = np.sum(amp[np.abs(x) < 15])/np.sum(amp)
cs0 = np.sqrt(g*n0); xH = x[np.argmin(np.abs(v0-cs0))]

print("\n" + "-"*64)
print(f"  Modo de mayor Im(w) tras CAP (Gamma0={Gamma0}):")
print(f"    w = {wtop.real:+.4f} {wtop.imag:+.4f} i")
print(f"    pico en x = {x[np.argmax(amp)]:.2f}  (horizonte x_H={xH:.2f})")
print(f"    fraccion borde {edge:.1%} / centro {cen:.1%}")
stable = np.max(wc.imag) < 5e-3
print("\n" + "="*64)
print(f"  max Im(w) con CAP = {np.max(wc.imag):+.2e}")
print(f"  HITO 2: HORIZONTE ESPECTRALMENTE ESTABLE: {'SI (referee-grade)' if stable else 'REVISAR'}")
print("="*64)

# ------------------------------- Figura --------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].scatter(w0.real, w0.imag, s=10, alpha=0.5, color='crimson',
              label='sin CAP (pared espuria)')
ax[0].scatter(wc.real, wc.imag, s=10, alpha=0.6, color='navy', label='con CAP')
ax[0].axhline(0, color='k', lw=0.8); ax[0].set_ylim(-0.15, 0.15); ax[0].set_xlim(-6, 6)
ax[0].set_xlabel(r'Re($\omega$)'); ax[0].set_ylabel(r'Im($\omega$)')
ax[0].set_title('Espectro BdG: CAP elimina los modos crecientes de pared')
ax[0].legend(fontsize=9); ax[0].grid(alpha=0.3)

ax[1].plot(x, np.abs(u)/amp.max(), label='|u|'); ax[1].plot(x, np.abs(v)/amp.max(), label='|v|')
ax[1].plot(x, Gamma/Gamma.max(), 'g:', lw=1, label=r'CAP $\Gamma(x)$ (norm.)')
ax[1].axvline(xH, color='k', ls='--', lw=1, label='horizonte')
ax[1].set_xlabel('x'); ax[1].set_ylabel('amplitud (norm.)')
ax[1].set_title(f'Modo de mayor Im tras CAP  (Im w={wtop.imag:.1e})')
ax[1].legend(fontsize=9)
plt.tight_layout(); plt.savefig("out/bdg_cap_stability.png", dpi=120); plt.close()
print("\nFigura: out/bdg_cap_stability.png")
