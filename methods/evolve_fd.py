"""
Evolucionador independiente de diferencias finitas para probar la estabilidad
dinamica del horizonte sonico waterfall de forma RIGUROSA.

Integrador: Crank-Nicolson semi-implicito (lineal implicito + no-lineal explicito)
  -> incondicionalmente estable para la parte cinetica (RK4 explicito explotaba).

Claves (siguiendo la revision):
  - El potencial se define con EL MISMO laplaciano discreto D_xx del evolucionador:
        V_i = mu - g n_i + (1/2) Re[(D_xx psi0)_i / psi0_i]
    asi psi0 es estacionario a precision de maquina CONTRA ESTE operador.
  - Bordes: se IMPONE el fondo estacionario psi0 (Dirichlet), que alimenta el
    flujo correcto sin absorberlo (nada de potencial absorbente).
  - Diagnosticos: densidad, corriente, x_H(t), y mapa espacio-tiempo Delta n(x,t).

Unidades hbar=m=1.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt

# ----------------------- Cargar el estado estacionario -----------------------
d = np.load("out/waterfall_state.npz")
x, n0, v0 = d["x"], d["n"], d["v"]
j, mu, g, n_H = float(d["j"]), float(d["mu"]), float(d["g"]), float(d["n_H"])
N = x.size
dx = x[1] - x[0]

theta = np.cumsum(v0) * dx
theta -= theta[N // 2]
psi0 = np.sqrt(n0) * np.exp(1j * theta)

# -------------------- Laplaciano discreto (2o orden, matriz) -----------------
main = -2.0 * np.ones(N)
off = np.ones(N - 1)
Dxx = sp.diags([off, main, off], [-1, 0, 1], format="csr") / dx**2

# ---------- Potencial consistente con ESTE D_xx (la correccion clave) --------
V = mu - g * n0 + 0.5 * np.real((Dxx @ psi0) / psi0)

# Residuo estacionario contra el operador DISCRETO del evolucionador
H_psi0 = -0.5 * (Dxx @ psi0) + (V + g * np.abs(psi0)**2) * psi0
core = np.abs(x) < 0.85 * x.max()
res0 = np.max(np.abs((H_psi0 - mu * psi0)[core])) / np.max(np.abs(mu * psi0[core]))
print("="*64)
print("EVOLUCIONADOR FD (Crank-Nicolson) — estabilidad rigurosa")
print("="*64)
print(f"  Residuo estacionario (operador discreto del evolucionador) = {res0:.2e}")

# ----------------- Operadores Crank-Nicolson (prefactorizados) ---------------
# Lineal L = -1/2 D_xx + diag(V).  (I + i dt/2 L) psi^{n+1} = (I - i dt/2 L) psi^n
#                                   - i dt g|psi^n|^2 psi^n   (no lineal explicito)
T, dt = 20.0, 2e-3
nsteps = int(T / dt)
I = sp.identity(N, format="csr")
L = (-0.5 * Dxx + sp.diags(V, 0, format="csr")).tocsr()
A = (I + 0.5j * dt * L).tolil()
B = (I - 0.5j * dt * L).tocsr()

nb = 6
clamp = np.zeros(N, bool); clamp[:nb] = True; clamp[-nb:] = True
cidx = np.where(clamp)[0]
for i in cidx:                           # filas Dirichlet: psi^{n+1}_i = psi0_i
    A.rows[i] = [i]; A.data[i] = [1.0]
A = A.tocsc()
lu = spla.splu(A)                        # prefactorizar una sola vez

# Capa esponja: amortigua la PERTURBACION (psi-psi0) hacia el fondo, NO hacia 0,
# en los bordes -> deja salir ondas sin reflejar y sin absorber el flujo de fondo.
x_sp, Gamma0 = 30.0, 3.0
Gamma = np.where(np.abs(x) > x_sp,
                 Gamma0 * ((np.abs(x) - x_sp) / (x.max() - x_sp))**2, 0.0)

# ------------------------------ Evolucionar ----------------------------------
cs0 = np.sqrt(g * n0)
xH0 = x[np.argmin(np.abs(v0 - cs0))]
central = np.abs(x) < 10.0
nsnap = 400
snap_every = max(1, nsteps // nsnap)

psi = psi0.copy()
times, dn_map, xH_t, jerr_t, drift_t = [], [], [], [], []
for step in range(nsteps + 1):
    if step % snap_every == 0:
        n = np.abs(psi)**2
        dpsi = (np.roll(psi, -1) - np.roll(psi, 1)) / (2*dx)
        jx = np.imag(np.conj(psi) * dpsi)
        v = jx / n
        cs = np.sqrt(g * n)
        f = (v - cs)[np.abs(x) < 15]; xr = x[np.abs(x) < 15]
        ci = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0]
        xH = xr[ci[0]] if ci.size else np.nan
        times.append(step * dt)
        dn_map.append((n - n0)[np.abs(x) < 20])
        xH_t.append(xH)
        jerr_t.append(np.max(np.abs((n*v - j))[central]) / j)
        drift_t.append(np.max(np.abs(n - n0)[central]) / np.max(n0[central]))
    if step < nsteps:
        rhs = (B @ psi - 1j * dt * g * np.abs(psi)**2 * psi
               - dt * Gamma * (psi - psi0))      # esponja sobre la perturbacion
        rhs[clamp] = psi0[clamp]
        psi = lu.solve(rhs)

times = np.array(times)
drift_t = np.array(drift_t); xH_t = np.array(xH_t); jerr_t = np.array(jerr_t)

# ------------------------------ Resultados -----------------------------------
norm0 = np.sum(np.abs(psi0[core])**2)*dx
normT = np.sum(np.abs(psi[core])**2)*dx
print("\n" + "-"*64)
print(f"RESULTADOS (central |x|<10, t hasta {T:.0f}, dt={dt})")
print("-"*64)
print(f"  [1] Residuo estacionario     = {res0:.2e}")
print(f"  [2] Continuidad  max|nv-j|/j  = {np.nanmax(jerr_t):.2e}")
print(f"  [3] Deriva densidad max|dn|/n = {np.nanmax(drift_t):.2e}")
print(f"  [4] Horizonte  x_H: {np.nanmin(xH_t):.4f} .. {np.nanmax(xH_t):.4f}  "
      f"(movimiento {np.nanmax(xH_t)-np.nanmin(xH_t):.2e})")
print(f"  [5] Norma central  |dN|/N     = {abs(normT-norm0)/norm0:.2e}")

ok = (np.nanmax(drift_t) < 1e-2 and (np.nanmax(xH_t)-np.nanmin(xH_t)) < 1e-1
      and np.nanmax(jerr_t) < 1e-3)
print("\n" + "="*64)
print(f"  HORIZONTE DINAMICAMENTE ESTABLE: {'SI' if ok else 'NO (ver mapa)'}")
print("="*64)

# -------------------------- Figura: mapa + curvas ----------------------------
dn_map = np.array(dn_map).T
xmap = x[np.abs(x) < 20]
fig, ax = plt.subplots(1, 2, figsize=(13, 5), gridspec_kw={'width_ratios':[1.4,1]})
vmax = max(np.max(np.abs(dn_map)), 1e-12)
im = ax[0].pcolormesh(times, xmap, dn_map, shading='auto', cmap='RdBu_r',
                      vmin=-vmax, vmax=vmax)
ax[0].axhline(xH0, color='k', ls='--', lw=1)
ax[0].set_xlabel('t'); ax[0].set_ylabel('x')
ax[0].set_title(r'$\Delta n(x,t)=n(x,t)-n(x,0)$   (horizonte: linea negra)')
plt.colorbar(im, ax=ax[0], label=r'$\Delta n$')

ax[1].plot(times, drift_t, 'r-', label=r'deriva densidad max$|\Delta n|/n$')
ax[1].plot(times, jerr_t, 'b-', label=r'error continuidad')
ax[1].axhline(1e-2, color='gray', ls=':', lw=1, label='umbral 1%')
ax[1].set_xlabel('t'); ax[1].set_yscale('log'); ax[1].legend(fontsize=9)
ax[1].set_title('Diagnosticos vs tiempo (central)')
plt.tight_layout(); plt.savefig("out/evolve_fd_stability.png", dpi=120); plt.close()
print("\nFigura: out/evolve_fd_stability.png")
