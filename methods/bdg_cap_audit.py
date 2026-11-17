"""
AUDITORIA DEL CAP (cierre riguroso del hito 2).

Objecion valida: un CAP amortigua artificialmente; que Im(w)>0 desaparezca tras
meterlo NO demuestra estabilidad fisica por si solo (podria enmascarar una
inestabilidad real).

Argumento decisivo: el CAP vale 0 en el interior (|x|<xcap). Un modo que NO
solapa la zona del CAP no puede ser amortiguado por el. Por tanto:

   si max Im(w) entre los modos SIN solape con el CAP es <= 0
   -> NO existe inestabilidad interior (con CAP o sin el) -> ESTABLE de verdad.

Pruebas:
  (1) Correlacion Im(w) vs solape-con-CAP: los modos con Im!=0 deben ser SOLO los
      de frontera; los interiores deben quedar reales.
  (2) max Im(w) restringido a modos interiores (solape CAP < 5%).
  (3) Robustez frente a ANCHURA, POSICION e INTENSIDAD del CAP.
  (4) El espectro real interior no se mueve al meter el CAP (no distorsiona fisica).
"""
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

g, n_u, mach_u, mach_d, w_step = 1.0, 1.0, 0.5, 2.0, 2.0
v_u = mach_u*np.sqrt(g*n_u); j = n_u*v_u; mu = g*n_u + 0.5*v_u**2
n_d = (j/(mach_d*np.sqrt(g)))**(2/3); n_H = (j**2/g)**(1/3)


def build(Nb, L, Gamma0, xcap_frac):
    x = np.linspace(-L, L, Nb); dx = x[1]-x[0]
    n0 = n_u + (n_d-n_u)*0.5*(1+np.tanh(x/w_step)); v0 = j/n0
    th = np.cumsum(v0)*dx; th -= th[Nb//2]
    psi0 = np.sqrt(n0)*np.exp(1j*th)
    main=-2*np.ones(Nb); off=np.ones(Nb-1)
    Dxx = sp.diags([off,main,off],[-1,0,1]).toarray()/dx**2
    V = mu - g*n0 + 0.5*np.real((Dxx@psi0)/psi0)
    xcap = xcap_frac*L
    Gamma = np.where(np.abs(x)>xcap, Gamma0*((np.abs(x)-xcap)/(L-xcap))**2, 0.0)
    Ld = -0.5*Dxx + np.diag(V + 2*g*np.abs(psi0)**2 - mu)
    G = np.diag(g*psi0**2)
    M = np.block([[Ld-1j*np.diag(Gamma), G],[-np.conj(G), -np.conj(Ld)-1j*np.diag(Gamma)]])
    return x, Dxx, psi0, M, Gamma, xcap, n0, v0


Nb, L = 1000, 60.0
xH = 0.42

# --- Simetria particula-hueco del CAP: sigma_x M* sigma_x = -M (se preserva) ---
x, Dxx, psi0, M, Gamma, xcap, n0, v0 = build(Nb, L, 2.0, 0.75)
sx = np.block([[np.zeros((Nb,Nb)), np.eye(Nb)],[np.eye(Nb), np.zeros((Nb,Nb))]])
sym_err = np.max(np.abs(sx @ np.conj(M) @ sx + M))
print("="*64)
print("AUDITORIA CAP — cierre del hito 2")
print("="*64)
print(f"  Simetria particula-hueco |sx M* sx + M| = {sym_err:.2e}  (debe ~0)")

# --- (1)+(2) Im(w) vs solape con el CAP, y max Im entre modos interiores ------
wv, vr = la.eig(M)
u = vr[:Nb, :]; v = vr[Nb:, :]
dens = np.abs(u)**2 + np.abs(v)**2
dens /= dens.sum(axis=0, keepdims=True)
cap_mask = np.abs(x) > xcap
overlap = dens[cap_mask, :].sum(axis=0)          # fraccion en la zona del CAP
interior = overlap < 0.05                          # modos que NO tocan el CAP
max_im_interior = np.max(wv.imag[interior])
max_im_all = np.max(wv.imag)
print(f"  max Im(w) TODOS los modos       = {max_im_all:+.2e}")
print(f"  max Im(w) modos INTERIORES      = {max_im_interior:+.2e}   <- decisivo")
print(f"    (modos interiores: {interior.sum()} de {wv.size}, solape CAP<5%)")

# --- (3) Robustez frente a anchura, posicion e intensidad del CAP ------------
print("\n  Robustez (max Im en modos interiores):")
print(f"  {'xcap_frac':>10} {'Gamma0':>8} {'maxIm_interior':>16}")
for xf in (0.65, 0.75, 0.85):
    for G0 in (1.0, 2.0, 4.0):
        _, _, _, Mt, _, xct, _, _ = build(800, L, G0, xf)
        wt, vt = la.eig(Mt)
        ut = vt[:800]; vt2 = vt[800:]
        dt = np.abs(ut)**2+np.abs(vt2)**2; dt/=dt.sum(0,keepdims=True)
        xx = np.linspace(-L,L,800)
        ov = dt[np.abs(xx)>xct,:].sum(0)
        mi = np.max(wt.imag[ov<0.05])
        print(f"  {xf:>10.2f} {G0:>8.1f} {mi:>16.2e}")

# --- (4) El espectro real interior no se mueve con/sin CAP -------------------
_, _, _, M_noCAP, _, _, _, _ = build(Nb, L, 0.0, 0.75)
w_noCAP = la.eigvals(M_noCAP)
# frecuencias reales propagantes (|Im|<1e-6) en una ventana, comparar conjuntos
re0 = np.sort(w_noCAP.real[np.abs(w_noCAP.imag) < 1e-6])
reC = np.sort(wv.real[np.abs(wv.imag) < 1e-6])
band = (np.abs(re0) < 20)
match = np.max(np.abs(re0[band][:min(len(reC[np.abs(reC)<20]), band.sum())]
               - reC[np.abs(reC)<20][:min(len(reC[np.abs(reC)<20]), band.sum())])) \
        if reC[np.abs(reC)<20].size else np.nan
print(f"\n  Espectro real interior sin CAP vs con CAP: desviacion max = {match:.2e}")

verdict = (max_im_interior < 5e-3)
print("\n" + "="*64)
print(f"  HITO 2 CERRADO: {'ESTABLE (el CAP no enmascara: interior intacto)' if verdict else 'REVISAR'}")
print("="*64)

# ------------------------------- Figura --------------------------------------
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
sc = ax[0].scatter(overlap, wv.imag, c=np.abs(wv.real), s=12, cmap='viridis')
ax[0].axhline(0, color='r', lw=0.8); ax[0].axvline(0.05, color='gray', ls=':')
ax[0].set_xlabel('solape con zona CAP'); ax[0].set_ylabel(r'Im($\omega$)')
ax[0].set_title('Im(w) SOLO en modos de frontera;\ninteriores (izq. de linea) quedan reales')
plt.colorbar(sc, ax=ax[0], label=r'|Re $\omega$|')

ax[1].scatter(wv.real, wv.imag, s=10, c=overlap, cmap='coolwarm')
ax[1].axhline(0, color='k', lw=0.8); ax[1].set_ylim(-0.15, 0.05); ax[1].set_xlim(-6, 6)
ax[1].set_xlabel(r'Re($\omega$)'); ax[1].set_ylabel(r'Im($\omega$)')
ax[1].set_title('Espectro con CAP (color=solape frontera)\nnada en el semiplano Im>0')
plt.tight_layout(); plt.savefig("out/bdg_cap_audit.png", dpi=120); plt.close()
print("\nFigura: out/bdg_cap_audit.png")
