"""
BARRIDO DE LA BANDA 0<omega<0.468 — beta_omega por DOS metodos independientes.
Objetivo: determinar EN QUE PARTE de la banda beta es numericamente fiable.
NO interpretar temperatura. Tabla por frecuencia con los controles de ChatGPT.

Para cada omega: beta^2_MM (mode-matching, FD), beta^2_RK4 (transfer), discrepancia
relativa, pseudo-unitaridad de cada metodo, nº de canales propagantes. Una
frecuencia se marca FIABLE solo si ambos pseudo-unit < 1e-2 y discrepancia < 5%.
Convergencia (N) en una seleccion de frecuencias.

Cautela registrada: ambos metodos comparten operador BdG, raices asintoticas y
convenciones de norma -> corroboracion cruzada, NO excluye error sistematico comun.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import mode_matching as MM
import bdg_smatrix as B


def beta_RK4(om):
    S, ei, eo = B.smatrix(om)[:3]
    pu = np.max(np.abs(S.conj().T@np.diag(eo.astype(float))@S - np.diag(ei.astype(float))))
    b2 = abs(S[np.where(eo < 0)[0][0], np.where(ei > 0)[0][0]])**2
    nprop = len(ei) + len(eo)      # canales in+out
    return b2, pu, nprop


def beta_MM(om, L=24, N=4000):
    MM.omega = om
    r = MM.solve(L=L, N=N, mode='horizon')
    pu = abs(r['flux_in'] - r['flux_out'])
    nprop = sum(m['prop'] for m in r['Laway']) + sum(m['prop'] for m in r['Raway']) + 1
    return r['beta2'], pu, nprop


print("="*94)
print("BARRIDO DE BANDA — beta por mode-matching (MM) y transfer (RK4). Fiable si ambos ok.")
print("="*94)
print(f"{'omega':>6} {'b2_MM':>10} {'b2_RK4':>10} {'d_rel':>9} {'eU_MM':>9} {'eU_RK4':>9} {'nch':>4}  estado")
print("-"*94)

omegas = np.concatenate([[0.1], np.linspace(0.03, 0.44, 19)])
omegas = np.unique(np.round(omegas, 4))
rows = []
for om in omegas:
    try:
        b_mm, pu_mm, nc_mm = beta_MM(om)
    except Exception as e:
        print(f"{om:6.3f}  MM fallo: {type(e).__name__}"); continue
    try:
        b_rk, pu_rk, nc_rk = beta_RK4(om)
    except Exception as e:
        print(f"{om:6.3f}  RK4 fallo: {type(e).__name__}"); b_rk, pu_rk, nc_rk = np.nan, np.nan, 0
    drel = abs(b_mm - b_rk)/max(abs(b_mm), abs(b_rk), 1e-12)
    ok = (pu_mm < 1e-2) and (pu_rk < 1e-2) and (drel < 0.05)
    estado = "FIABLE" if ok else "flag"
    rows.append((om, b_mm, b_rk, drel, pu_mm, pu_rk, nc_mm, ok))
    print(f"{om:6.3f} {b_mm:10.5f} {b_rk:10.5f} {drel:9.1e} {pu_mm:9.1e} {pu_rk:9.1e} {nc_mm:4d}  {estado}")

nfiable = sum(r[7] for r in rows)
print("-"*94)
print(f"  frecuencias FIABLES (ambos metodos + discrepancia<5%): {nfiable}/{len(rows)}")
if nfiable:
    fi = [r for r in rows if r[7]]
    print(f"  rango fiable: omega in [{min(r[0] for r in fi):.3f}, {max(r[0] for r in fi):.3f}]")
    print(f"  discrepancia MM-RK4 en zona fiable: max {max(r[3] for r in fi):.1e}, "
          f"mediana {np.median([r[3] for r in fi]):.1e}")

# ---- convergencia en N en frecuencias seleccionadas (baja, media, cerca de umbral)
print("\n  CONVERGENCIA EN N (beta^2_MM):")
for om in (0.05, 0.20, 0.40):
    vals = []
    for N in (2000, 4000, 6000):
        MM.omega = om
        vals.append(MM.solve(L=24, N=N, mode='horizon')['beta2'])
    print(f"    omega={om}: " + " -> ".join(f"{val:.5f}" for val in vals)
          + f"   (variacion {abs(vals[-1]-vals[0])/vals[-1]:.1e})")
print("\n  (sin interpretar temperatura; solo fiabilidad numerica de beta en la banda)")
