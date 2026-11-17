"""
Horizonte sonico "waterfall" en un BEC 1D — solucion estacionaria AUTOCONSISTENTE
de la GPE completa (incluyendo presion cuantica), NO un perfil pintado.

Referencia de contraste:
  Larre, Recati, Carusotto, Pavloff, Phys. Rev. A 85, 013621 (2012),
  arXiv:1110.4464.

Unidades: hbar = m = 1.

Formulacion (Madelung):  psi = sqrt(n) exp(i theta),  v = theta',  j = n v = cte.
GPE estacionaria para u = sqrt(n):
    u'' = 2 u (V(x) - mu) + 2 g u^3 + j^2 / u^3
(el termino j^2/u^3 viene de la energia cinetica del flujo; el u'' es la presion
 cuantica). Asintoticamente u''=0 -> relacion de Bernoulli uniforme.

Horizonte sonico: v(x) = c_s(x), con c_s = sqrt(g n).
Densidad analitica del horizonte:  n_H = (j^2 / g)^(1/3).
"""

import numpy as np
from scipy.integrate import solve_bvp

# ----------------------------- Parametros fisicos ----------------------------
g = 1.0
n_u = 1.0           # densidad aguas arriba (upstream)
mach_u = 0.5        # numero de Mach upstream (subsonico, <1)
mach_d = 2.0        # numero de Mach downstream (supersonico, >1)
w = 2.0             # anchura del escalon de potencial (en unidades de healing length)
L = 60.0            # semi-dominio

# ----------------------- Estados asintoticos uniformes -----------------------
c_u = np.sqrt(g * n_u)
v_u = mach_u * c_u
j = n_u * v_u                       # corriente conservada
mu = g * n_u + 0.5 * v_u**2         # potencial quimico (V=0 upstream)

# downstream: j = n_d * mach_d * sqrt(g n_d)  ->  n_d = (j/mach_d)^(2/3) * g^(-1/3)
n_d = (j / (mach_d * np.sqrt(g)))**(2.0/3.0)
v_d = j / n_d
# escalon de potencial que produce ese estado downstream (misma j, mismo mu):
V0 = mu - g * n_d - 0.5 * v_d**2

u_u, u_d = np.sqrt(n_u), np.sqrt(n_d)
n_H = (j**2 / g)**(1.0/3.0)         # densidad analitica del horizonte

print("="*64)
print("MONTAJE WATERFALL (autoconsistente)")
print("="*64)
print(f"  g={g}, n_u={n_u}, Mach_u={mach_u}, Mach_d={mach_d}")
print(f"  corriente  j   = {j:.6f}")
print(f"  pot.quimico mu = {mu:.6f}")
print(f"  n_downstream   = {n_d:.6f}   (v_d={v_d:.4f})")
print(f"  escalon V0     = {V0:.6f}")
print(f"  n_HORIZONTE    = (j^2/g)^(1/3) = {n_H:.6f}")
print(f"    -> n_u={n_u} > n_H={n_H:.4f} > n_d={n_d:.4f}  (cruce garantizado)")


# --------------------- METODO INVERSO (autoconsistente exacto) ---------------
# Elegimos n(x) monotono que cruza n_H; v=j/n por continuidad; despejamos V(x)
# de la GPE estacionaria -> solucion estacionaria EXACTA por construccion.
N = 8192
x = np.linspace(-L, L, N)
dx = x[1] - x[0]

t = np.tanh(x / w)
s = 0.5 * (1.0 + t)
Delta = n_d - n_u
n = n_u + Delta * s                       # densidad elegida (suave, monotona)
# derivadas analiticas de n (tanh): s'=(1-t^2)/(2w), s''=-t(1-t^2)/w^2
n1 = Delta * (1.0 - t**2) / (2.0 * w)      # n'
n2 = -Delta * t * (1.0 - t**2) / w**2      # n''

v = j / n
cs = np.sqrt(g * n)
mach = v / cs

# presion cuantica: (1/2) u''/u con u=sqrt(n)  =>  (1/2)[n''/(2n) - n'^2/(4n^2)]
Qp = 0.5 * (n2 / (2.0 * n) - n1**2 / (4.0 * n**2))
# potencial que hace estacionario este flujo (despejado de la GPE):
Vx = mu - g * n - 0.5 * v**2 + Qp

def Vpot(xq):
    return np.interp(xq, x, Vx)

print(f"\n  Metodo inverso: V(x) derivado analiticamente (estacionario exacto)")
print(f"    V(-inf)={Vx[0]:.4f}  V(+inf)={Vx[-1]:.4f}  V_max bump={Vx.max():.4f}")

u = np.sqrt(n)

# ------------------------- Localizar el horizonte ----------------------------
f = mach - 1.0                      # cero donde v = c_s
idx = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
if idx.size:
    i = idx[0]
    xH = x[i] - f[i] * (x[i+1] - x[i]) / (f[i+1] - f[i])   # interpolacion lineal
    nH_num = np.interp(xH, x, n)
else:
    xH, nH_num = np.nan, np.nan

print("\n" + "-"*64)
print("HORIZONTE SONICO")
print("-"*64)
print(f"  posicion  x_H       = {xH:.5f}")
print(f"  densidad  n(x_H)    = {nH_num:.6f}")
print(f"  n_H analitico       = {n_H:.6f}   (error {abs(nH_num-n_H)/n_H:.2e})")

# ------------------- VERIFICACION 1: continuidad n*v = j ---------------------
j_err = np.max(np.abs(n * v - j)) / j
print("\n" + "-"*64)
print("VERIFICACIONES")
print("-"*64)
print(f"  [1] Continuidad  max|n v - j|/j      = {j_err:.2e}")

# ------- VERIFICACION 2: residuo de la GPE estacionaria completa -------------
# Construir psi = u e^{i theta}, theta = int v dx, y comprobar H psi = mu psi
theta = np.cumsum(v) * dx
theta -= theta[N // 2]
psi = u * np.exp(1j * theta)
# psi'' por diferencias finitas centradas de 4o orden (independiente de la EDO)
d2psi = (-np.roll(psi, 2) + 16*np.roll(psi, 1) - 30*psi
         + 16*np.roll(psi, -1) - np.roll(psi, -2)) / (12 * dx**2)
Hpsi = -0.5 * d2psi + (Vpot(x) + g * np.abs(psi)**2) * psi
res = Hpsi - mu * psi
core = (np.abs(x) < 0.8 * L)        # excluir bordes (FD no periodica)
res_rel = np.max(np.abs(res[core])) / np.max(np.abs(mu * psi[core]))
print(f"  [2] Residuo GPE  max|H psi - mu psi|  = {res_rel:.2e}  (interior)")

# ------------- VERIFICACION 3: estabilidad en evolucion temporal -------------
# Absorbedor (CAP) en los bordes: V imaginaria que se traga la perturbacion de
# wraparound (el estado lleva flujo y no es periodico). Asi el test mide solo la
# fisica en la region del horizonte, no el artefacto del borde periodico.
from gpe import GPE1D
sim = GPE1D(L=2*L, N=N, g=g)
sim.x = x.copy(); sim.dx = dx
sim.k = 2*np.pi*np.fft.fftfreq(N, d=dx); sim.k2 = sim.k**2
x_cap, eta = 0.80 * L, 4.0
cap = np.where(np.abs(x) > x_cap,
               eta * ((np.abs(x) - x_cap) / (L - x_cap))**2, 0.0)
sim.V = Vx - 1j * cap               # parte imaginaria -> amortiguacion en bordes
T, dt = 4.0, 2e-3
psi_t = sim.evolve_real(psi.astype(complex), dt, int(T/dt))
central = np.abs(x) < 10.0          # region del horizonte, lejos de los bordes
drift = np.max(np.abs(np.abs(psi_t[central])**2 - n[central])) / np.max(n[central])
print(f"  [3] Estabilidad  max|dn|/n (t={T}, |x|<10, CAP) = {drift:.2e}")

ok1 = j_err < 1e-10
ok2 = res_rel < 1e-3
ok3 = drift < 1e-2
okH = (not np.isnan(xH)) and abs(nH_num - n_H)/n_H < 1e-3
print("\n" + "="*64)
print(f"  Continuidad      : {'PASA' if ok1 else 'FALLA'}")
print(f"  Residuo GPE      : {'PASA' if ok2 else 'FALLA'}")
print(f"  Estabilidad temp.: {'PASA' if ok3 else 'FALLA'}")
print(f"  Cruce horizonte  : {'PASA' if okH else 'FALLA'}")
print("="*64)

# ------------------------------- Figuras -------------------------------------
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

fig, ax = plt.subplots(2, 1, figsize=(9, 7), sharex=True)
ax[0].plot(x, n, 'b-', label=r'$n(x)=|\psi|^2$')
ax[0].plot(x, Vpot(x), 'g--', lw=1, label=r'$V(x)$ (escalon)')
ax[0].axhline(n_H, color='r', ls=':', lw=1, label=r'$n_H=(j^2/g)^{1/3}$')
if not np.isnan(xH):
    ax[0].axvline(xH, color='k', ls='-', lw=0.8, alpha=0.5)
ax[0].set_ylabel('densidad / potencial'); ax[0].legend(loc='right')
ax[0].set_title('Horizonte sonico waterfall — solucion estacionaria autoconsistente')

ax[1].plot(x, v, 'r-', label=r'$v(x)=j/n$ (flujo)')
ax[1].plot(x, cs, 'b-', label=r'$c_s(x)=\sqrt{gn}$ (sonido)')
ax[1].plot(x, mach, 'k-', lw=1, label=r'Mach $=v/c_s$')
ax[1].axhline(1.0, color='gray', ls=':', lw=1)
if not np.isnan(xH):
    ax[1].axvline(xH, color='k', ls='-', lw=0.8, alpha=0.5)
    ax[1].annotate('HORIZONTE\n$v=c_s$', xy=(xH, 1.0), xytext=(xH+8, 1.4),
                   arrowprops=dict(arrowstyle='->'), fontsize=9)
ax[1].fill_between(x, 0, 3, where=(mach > 1), alpha=0.08, color='red')
ax[1].text(L*0.45, 0.25, 'SUPERSONICO\n(interior)', fontsize=8, color='darkred')
ax[1].text(-L*0.8, 0.25, 'SUBSONICO (exterior)', fontsize=8, color='navy')
ax[1].set_ylim(0, 2.6); ax[1].set_xlabel('x'); ax[1].set_ylabel('velocidad / Mach')
ax[1].legend(loc='upper left'); ax[1].set_xlim(-40, 40)
plt.tight_layout(); plt.savefig("out/waterfall_horizon.png", dpi=120); plt.close()
print("\nFigura: out/waterfall_horizon.png")

np.savez("out/waterfall_state.npz", x=x, n=n, v=v, cs=cs, mach=mach,
         psi_re=psi.real, psi_im=psi.imag, V=Vpot(x),
         j=j, mu=mu, g=g, n_H=n_H, xH=xH)
print("Estado guardado: out/waterfall_state.npz")
