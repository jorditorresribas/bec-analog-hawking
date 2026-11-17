"""
Espectro de Bogoliubov-de Gennes (BdG) alrededor del fondo waterfall, para
determinar la ESTABILIDAD ESPECTRAL del horizonte sonico de forma rigurosa
(sin ambiguedades de fronteras de la evolucion temporal).

Convencion temporal  psi = e^{-i mu t}[ psi0 + u e^{-i w t} + v* e^{i w* t} ].
Problema de autovalores (un unico operador Dxx, 2o orden -> precaucion A):
    [[ L ,   g psi0^2   ],  [u]        [u]
     [ -g psi0*^2 , -L* ]]  [v]  =  w  [v]
con  L = -1/2 Dxx + V + 2 g |psi0|^2 - mu.

Lectura (convencion e^{-i w t}):
    Im(w) = 0  -> oscilatorio (estable)
    Im(w) > 0  -> inestable (crece)        -> max|Im w| es el indicador clave
    w ~ 0      -> modo Goldstone de fase

Precaucion B: se mide la parte IMAGINARIA de Dxx psi0 / psi0 (descartada al
definir V real). Si no es despreciable, psi0 no es estacionario exacto con V real.
Precaucion (convergencia): se corre a dos resoluciones y se compara max|Im w|.
"""
import numpy as np
import scipy.sparse as sp
import scipy.linalg as la
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------------- Parametros (= waterfall) ----------------------
g, n_u, mach_u, mach_d, w_step, L = 1.0, 1.0, 0.5, 2.0, 2.0, 60.0
c_u = np.sqrt(g * n_u); v_u = mach_u * c_u
j = n_u * v_u
mu = g * n_u + 0.5 * v_u**2
n_d = (j / (mach_d * np.sqrt(g)))**(2.0/3.0)
n_H = (j**2 / g)**(1.0/3.0)


def build_background(Nb):
    """Fondo analitico waterfall en malla de Nb puntos + V con UN solo Dxx (2o)."""
    x = np.linspace(-L, L, Nb)
    dx = x[1] - x[0]
    t = np.tanh(x / w_step)
    n0 = n_u + (n_d - n_u) * 0.5 * (1.0 + t)
    v0 = j / n0
    theta = np.cumsum(v0) * dx; theta -= theta[Nb // 2]
    psi0 = np.sqrt(n0) * np.exp(1j * theta)
    main = -2.0*np.ones(Nb); off = np.ones(Nb-1)
    Dxx = sp.diags([off, main, off], [-1, 0, 1]).toarray() / dx**2
    ratio = (Dxx @ psi0) / psi0
    V = mu - g*n0 + 0.5*np.real(ratio)
    im_part = np.imag(ratio)             # precaucion B
    return x, dx, n0, v0, psi0, Dxx, V, im_part


def bdg_spectrum(Nb):
    x, dx, n0, v0, psi0, Dxx, V, im_part = build_background(Nb)
    # Precaucion B: parte imaginaria del cociente (interior, lejos de bordes)
    core = np.abs(x) < 0.85 * L
    imB = np.max(np.abs(im_part[core]))
    # Residuo complejo completo de la GPE estacionaria (precaucion del registro)
    Hpsi = -0.5*(Dxx @ psi0) + (V + g*np.abs(psi0)**2)*psi0
    res_cplx = np.max(np.abs((Hpsi - mu*psi0)[core])) / np.max(np.abs(mu*psi0[core]))

    # Matriz BdG 2Nb x 2Nb
    Ld = -0.5*Dxx + np.diag(V + 2*g*np.abs(psi0)**2 - mu)
    G = np.diag(g * psi0**2)
    M = np.block([[Ld, G], [-np.conj(G), -np.conj(Ld)]])
    wvals = la.eigvals(M)                 # no hermitiano -> autovalores complejos

    # Indicadores
    im_sorted = np.sort(np.abs(wvals.imag))[::-1]
    max_imw = im_sorted[0]
    gold = np.min(np.abs(wvals))          # modo mas cercano a w=0 (Goldstone)
    return dict(x=x, n0=n0, v0=v0, wvals=wvals, imB=imB, res=res_cplx,
                max_imw=max_imw, gold=gold, Nb=Nb)


print("="*64)
print("ESPECTRO BdG — estabilidad espectral del horizonte waterfall")
print("="*64)
print(f"  n_horizonte = {n_H:.5f}   (contraste analitico)")

# --------- Precaucion B + residuo, y convergencia a dos resoluciones ---------
runs = {}
for Nb in (800, 1200):
    print(f"\n--- Nb = {Nb} ---")
    r = bdg_spectrum(Nb)
    runs[Nb] = r
    print(f"  [B] max|Im(Dxx psi0/psi0)| (interior) = {r['imB']:.2e}")
    print(f"      residuo complejo GPE estacionaria = {r['res']:.2e}")
    print(f"  max|Im(omega)|  = {r['max_imw']:.3e}   <- indicador de inestabilidad")
    print(f"  min|omega|      = {r['gold']:.3e}   <- modo Goldstone (debe ~0)")

# Convergencia de max|Im w| con la resolucion
c1, c2 = runs[800]['max_imw'], runs[1200]['max_imw']
print("\n" + "-"*64)
print("CONVERGENCIA de max|Im(omega)|")
print("-"*64)
print(f"  Nb=800 : {c1:.3e}")
print(f"  Nb=1200: {c2:.3e}")
trend = "DECRECE (artefacto de discretizacion -> estable)" if c2 < c1*0.7 else \
        ("ESTABLE/CRECE con N (posible modo fisico)" if c2 > c1*1.3 else
         "~constante (revisar localizacion de modos)")
print(f"  tendencia: {trend}")

# ------------------------------- Figuras -------------------------------------
r = runs[1200]; wv = r['wvals']
fig, ax = plt.subplots(1, 2, figsize=(13, 5))
ax[0].scatter(wv.real, wv.imag, s=6, alpha=0.5, color='navy')
ax[0].axhline(0, color='r', lw=0.8)
ax[0].set_xlabel(r'Re($\omega$)'); ax[0].set_ylabel(r'Im($\omega$)')
ax[0].set_title(f'Espectro complejo BdG (Nb={r["Nb"]})\nmax|Im w|={r["max_imw"]:.2e}')
ax[0].grid(alpha=0.3)

idx = np.argsort(wv.real)
ax[1].plot(np.sort(wv.real), '.', ms=3, color='navy', label=r'Re($\omega$)')
ax[1].plot(np.sort(np.abs(wv.imag)), '.', ms=3, color='crimson',
           label=r'|Im($\omega$)|')
ax[1].set_xlabel('indice modal (ordenado)'); ax[1].set_ylabel(r'$\omega$')
ax[1].legend(); ax[1].set_title('Frecuencias vs indice modal')
ax[1].grid(alpha=0.3)
plt.tight_layout(); plt.savefig("out/bdg_spectrum.png", dpi=120); plt.close()
print("\nFigura: out/bdg_spectrum.png")
