"""
FASE 3B — Espectro de radiacion de Hawking analoga desde la matriz S.

Produccion espontanea en un canal saliente de norma POSITIVA = suma de |S|^2
sobre los canales entrantes de norma NEGATIVA:
    N_i^out = sum_{j: eta_in(j)=-1} |S_ij|^2    (coeficiente de Bogoliubov |beta|^2)

Se barre 0<omega<omega_max=0.468, se verifica S^dag eta S = eta en cada omega,
y se comprueba (SIN imponerlo) si N_omega es termico a baja frecuencia:
    termico <=> N_omega = 1/(exp(omega/T_H) - 1)  =>  log(1+1/N) = omega/T_H  (lineal)
Se compara la T_H ajustada con la prediccion de gravedad superficial
    k_sg = (1/2) d/dx (c_s - v)|horizonte,   T_H = k_sg/(2 pi).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from bdg_smatrix import smatrix, xg, n0g, v0g, g

X0 = 12.0
omega_max = 0.468
omegas = np.linspace(0.02, 0.44, 40)

N_Lout, N_Rout, N_tot, resid, oms = [], [], [], [], []
for om in omegas:
    try:
        S, eta_in, eta_out, _ = smatrix(om, X0=X0)
    except StopIteration:
        continue
    Eout = np.diag(eta_out.astype(float)); Ein = np.diag(eta_in.astype(float))
    pu = np.max(np.abs(S.conj().T @ Eout @ S - Ein))
    jneg = np.where(eta_in < 0)[0]            # canal(es) IN de norma negativa
    # produccion en canales OUT de norma positiva
    ipos = np.where(eta_out > 0)[0]
    Ni = np.array([sum(abs(S[i, j])**2 for j in jneg) for i in ipos])
    oms.append(om); resid.append(pu)
    N_Lout.append(Ni[0]); N_Rout.append(Ni[1] if len(Ni) > 1 else np.nan)
    N_tot.append(Ni.sum())

oms = np.array(oms); resid = np.array(resid)
N_Lout = np.array(N_Lout); N_Rout = np.array(N_Rout); N_tot = np.array(N_tot)
good = resid < 5e-3                            # solo frecuencias con S fiable

print("="*68)
print("ESPECTRO DE HAWKING ANALOGO")
print("="*68)
print(f"  frecuencias con S fiable (resid<5e-3): {good.sum()}/{len(oms)}")
print(f"  residuo pseudo-unit mediano: {np.median(resid[good]):.2e}")

# --- gravedad superficial (prediccion teorica de T_H) ------------------------
cs = np.sqrt(g*n0g)
f = cs - v0g
iH = np.argmin(np.abs(f))                      # horizonte: c_s = v
ksg = 0.5*abs(np.gradient(f, xg)[iH])          # k_sg = (1/2)|d(c-v)/dx|
TH_sg = ksg/(2*np.pi)
print(f"\n  gravedad superficial k_sg = {ksg:.4f}")
print(f"  T_H (prediccion grav. sup.) = k_sg/2pi = {TH_sg:.4f}")

# --- ajuste termico (SIN imponer): log(1+1/N) = omega/T_H -------------------
om_f = oms[good]; N_f = N_tot[good]
y = np.log(1 + 1/N_f)                           # = omega/T_H si termico
lo = om_f < 0.22                                # regimen baja frecuencia
A = np.polyfit(om_f[lo], y[lo], 1)
TH_fit = 1/A[0]
print(f"\n  AJUSTE TERMICO (baja omega, sin imponer):")
print(f"  log(1+1/N) vs omega: pendiente {A[0]:.3f} -> T_H_fit = {TH_fit:.4f}")
print(f"  intercepto = {A[1]:.4f} (termico puro -> ~0)")
print(f"\n  COMPARACION:  T_H_fit={TH_fit:.4f}  vs  T_H_grav.sup.={TH_sg:.4f}"
      f"   (ratio {TH_fit/TH_sg:.2f})")
resterm = np.max(np.abs(y[lo] - np.polyval(A, om_f[lo])))
print(f"  residuo del ajuste lineal (baja omega): {resterm:.2e}")

# ------------------------------- Figuras -------------------------------------
fig, ax = plt.subplots(1, 3, figsize=(16, 5))
ax[0].semilogy(oms[good], N_tot[good], 'ko-', ms=4, label=r'$N_\omega$ total')
ax[0].semilogy(oms[good], N_Lout[good], 'b.-', ms=3, label='upstream (+)')
ax[0].semilogy(oms[good], N_Rout[good], 'r.-', ms=3, label='downstream (+)')
om_th = np.linspace(0.02, 0.44, 100)
ax[0].semilogy(om_th, 1/(np.exp(om_th/TH_sg)-1), 'g--', lw=1.5,
               label=f'termico T_H={TH_sg:.3f} (grav.sup.)')
ax[0].axvline(omega_max, color='gray', ls=':'); ax[0].set_xlabel(r'$\omega$')
ax[0].set_ylabel(r'$N_\omega=|\beta_\omega|^2$'); ax[0].legend(fontsize=8)
ax[0].set_title('Espectro de Hawking analogo'); ax[0].grid(alpha=0.3)

ax[1].plot(om_f, y, 'ko-', ms=4)
ax[1].plot(om_f[lo], np.polyval(A, om_f[lo]), 'r-',
           label=f'ajuste: T_H={TH_fit:.3f}')
ax[1].set_xlabel(r'$\omega$'); ax[1].set_ylabel(r'$\log(1+1/N_\omega)$')
ax[1].set_title('Test termico (lineal = Planck)\nsin imponer'); ax[1].legend(fontsize=9)
ax[1].grid(alpha=0.3)

ax[2].semilogy(oms, resid, 'mo-', ms=3); ax[2].axhline(5e-3, color='gray', ls=':')
ax[2].set_xlabel(r'$\omega$'); ax[2].set_ylabel(r'residuo $S^\dagger\eta S=\eta$')
ax[2].set_title('Pseudo-unitaridad vs frecuencia'); ax[2].grid(alpha=0.3)
plt.tight_layout(); plt.savefig("out/hawking_spectrum.png", dpi=120); plt.close()
print("\nFigura: out/hawking_spectrum.png")
