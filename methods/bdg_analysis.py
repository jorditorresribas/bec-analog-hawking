"""
Diagnostico de los modos BdG complejos: distinguir INESTABILIDAD FISICA de
ARTEFACTO de frontera/discretizacion.

Dos pruebas (las que faltaban):
  1) Localizacion espacial del modo mas inestable: si vive en los BORDES de la
     caja -> espurio; si vive en el HORIZONTE -> candidato fisico.
  2) Dependencia con el TAMANO DE CAJA L: si max|Im w| cambia con L -> artefacto
     de frontera (un sistema abierto real no depende de donde pongas la pared).
"""
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g, n_u, mach_u, mach_d, w_step = 1.0, 1.0, 0.5, 2.0, 2.0
c_u = np.sqrt(g*n_u); v_u = mach_u*c_u
j = n_u*v_u; mu = g*n_u + 0.5*v_u**2
n_d = (j/(mach_d*np.sqrt(g)))**(2/3); n_H = (j**2/g)**(1/3)


def bdg_matrix(Nb, L):
    x = np.linspace(-L, L, Nb); dx = x[1]-x[0]
    n0 = n_u + (n_d-n_u)*0.5*(1+np.tanh(x/w_step))
    v0 = j/n0
    th = np.cumsum(v0)*dx; th -= th[Nb//2]
    psi0 = np.sqrt(n0)*np.exp(1j*th)
    main=-2*np.ones(Nb); off=np.ones(Nb-1)
    Dxx = sp.diags([off,main,off],[-1,0,1]).toarray()/dx**2
    V = mu - g*n0 + 0.5*np.real((Dxx@psi0)/psi0)
    Ld = -0.5*Dxx + np.diag(V + 2*g*np.abs(psi0)**2 - mu)
    G = np.diag(g*psi0**2)
    M = np.block([[Ld, G],[-np.conj(G), -np.conj(Ld)]])
    return x, n0, v0, psi0, M


# ---------- PRUEBA 1: localizacion espacial del modo mas inestable -----------
print("="*64)
print("PRUEBA 1 — localizacion del modo mas inestable (Nb=1000, L=60)")
print("="*64)
Nb, L = 1000, 60.0
x, n0, v0, psi0, M = bdg_matrix(Nb, L)
wv, vr = la.eig(M)
imax = np.argmax(wv.imag)
w_unst = wv[imax]
u = vr[:Nb, imax]; v = vr[Nb:, imax]
amp = np.abs(u) + np.abs(v)
amp /= amp.max()
xpeak = x[np.argmax(amp)]
# fraccion de la "masa" del modo en los bordes (|x|>0.8L) vs centro (|x|<15)
edge = np.sum(amp[np.abs(x) > 0.8*L]) / np.sum(amp)
cen = np.sum(amp[np.abs(x) < 15]) / np.sum(amp)
cs0 = np.sqrt(g*n0); xH = x[np.argmin(np.abs(v0-cs0))]
print(f"  omega mas inestable = {w_unst.real:.4f} + {w_unst.imag:.4f} i")
print(f"  pico del modo en x  = {xpeak:.2f}   (horizonte en x_H={xH:.2f})")
print(f"  fraccion en bordes |x|>0.8L = {edge:.1%}")
print(f"  fraccion en centro |x|<15   = {cen:.1%}")
verdict1 = "BORDE -> ESPURIO" if edge > cen else "CENTRO/HORIZONTE -> revisar"
print(f"  -> {verdict1}")

# ---------- PRUEBA 2: dependencia con el tamano de caja L ---------------------
print("\n" + "="*64)
print("PRUEBA 2 — max|Im(omega)| vs tamano de caja L (Nb fijo=800)")
print("="*64)
for Ltest in (40.0, 60.0, 80.0, 100.0):
    _, _, _, _, Mt = bdg_matrix(800, Ltest)
    wt = la.eigvals(Mt)
    print(f"  L={Ltest:5.0f} :  max|Im w| = {np.max(np.abs(wt.imag)):.3e}")
print("  (si cambia con L -> artefacto de frontera; si es constante -> fisico)")

# ------------------------------- Figura --------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].plot(x, np.abs(u)/amp.max(), label=r'$|u(x)|$')
ax[0].plot(x, np.abs(v)/amp.max(), label=r'$|v(x)|$')
ax[0].axvline(xH, color='k', ls='--', lw=1, label='horizonte')
ax[0].axvline(0.8*L, color='r', ls=':', lw=1, label=r'$0.8L$ (borde)')
ax[0].axvline(-0.8*L, color='r', ls=':', lw=1)
ax[0].set_xlabel('x'); ax[0].set_ylabel('amplitud (norm.)')
ax[0].set_title(f'Modo mas inestable  w={w_unst.real:.2f}+{w_unst.imag:.3f}i')
ax[0].legend(fontsize=9)

ax[1].scatter(wv.real, wv.imag, s=6, alpha=0.5, color='navy')
ax[1].scatter([w_unst.real], [w_unst.imag], s=80, color='red', marker='x',
              label='mas inestable')
ax[1].axhline(0, color='r', lw=0.8); ax[1].set_xlim(-5, 5)
ax[1].set_xlabel(r'Re($\omega$)'); ax[1].set_ylabel(r'Im($\omega$)')
ax[1].set_title('Espectro BdG (zoom cerca de 0)'); ax[1].legend(fontsize=9)
ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.savefig("out/bdg_mode_localization.png", dpi=120); plt.close()
print("\nFigura: out/bdg_mode_localization.png")
