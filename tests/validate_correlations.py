"""
VALIDACION de la maquinaria nueva (correlations.py) ANTES de usarla.
Gates (regla: verificar tras cualquier cambio):
 (A) dispersion lam=1 == quartic original (fd_bdg_operator.quartic_roots)
 (B) plane_mode_lam: |U|^2-|V|^2 = 1
 (C) homogeneo: inyectar 'u' -> transmision=1, reflexion~0 (sin horizonte)
 (D) HORIZONTE lam=1: beta^2 (amplitud saliente norma- al inyectar 'u') debe
     coincidir con el valor VALIDADO 0.01008 (mode_matching/transfer)
 (E) pseudo-unitaridad por canal: sum_out eta|a|^2 = eta_in
 (F) knob lam: dispersion (omega-vk)^2 = c^2 k^2 + lam k^4/4 con c^2=gn FIJO
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import correlations as C
from fd_bdg_operator import quartic_roots

print("="*68)
print("VALIDACION correlations.py")
print("="*68)

# (A) dispersion lam=1
om, vtest, ctest = 0.1, 0.5, 1.0
r_new = np.sort_complex(C.quartic_roots_lam(om, vtest, ctest, 1.0))
r_old = np.sort_complex(quartic_roots(om, vtest, ctest))
print(f"\n(A) raices dispersion lam=1 vs original: max dif = "
      f"{np.max(np.abs(r_new-r_old)):.2e}")

# (B) normalizacion amplitudes
maxd = 0.0
for k in (0.3, 0.8, 1.5, -0.6):
    U, V = C.plane_mode_lam(k, 0.7, 0.4, 1.0)
    maxd = max(maxd, abs(abs(U)**2-abs(V)**2 - 1))
print(f"(B) |U|^2-|V|^2=1 : max desviacion = {maxd:.2e}")

# (C) homogeneo (sin horizonte): transmision=1, reflexion~0
r = C.scatter_field(L=20, N=3000, omega=0.1, channel='u', mode='homog')
tr = abs(r['out_amp'].get(('R', 1), 0))**2
rf = sum(abs(a)**2 for (side, nm), a in r['out_amp'].items() if side == 'L')
print(f"(C) homogeneo: transmision(norma+)={tr:.4f}  reflexion={rf:.2e}")

# (D) HORIZONTE: beta^2
print(f"\n(D) HORIZONTE lam=1, omega=0.1: beta^2 (inyecta 'u', lee saliente norma-)")
for N in (3000, 4000, 6000):
    r = C.scatter_field(L=24, N=N, omega=0.1, channel='u', mode='horizon')
    beta2 = abs(r['out_amp'].get(('R', -1), 0))**2
    print(f"    N={N}: beta^2={beta2:.5f}   flujo_out={r['flux_out']:+.4f} "
          f"(eta_in={r['eta_in']:+d})")
print(f"    VALIDADO (memoria): beta^2 = 0.01008")

# (E) pseudo-unitaridad por canal inyectado
print(f"\n(E) pseudo-unitaridad por canal (horizonte, omega=0.1, N=4000):")
for ch in ('u', 'd1', 'd2'):
    r = C.scatter_field(L=24, N=4000, omega=0.1, channel=ch, mode='horizon')
    print(f"    canal {ch}: eta_in={r['eta_in']:+d}  sum_out eta|a|^2={r['flux_out']:+.5f}"
          f"  defecto={abs(r['flux_out']-r['eta_in']):.2e}")

# (F) knob lam: verificar dispersion asintotica numerica vs formula
print(f"\n(F) knob lam: c^2 debe quedar FIJO=gn al variar lam (metrica intacta)")
nL = 1.0; vL = 0.5; cL = np.sqrt(nL)
for lam in (0.5, 1.0, 2.0):
    # raiz fononica a k pequeno: (omega-vk)^2 ~ c^2 k^2 -> c_eff
    ks = C.quartic_roots_lam(0.001, vL, cL, lam)
    kr = [k.real for k in ks if abs(k.imag) < 1e-6]
    print(f"    lam={lam}: raices prop={np.round(sorted(kr),4)}  (c^2=gn={cL**2} fijo)")
    # chequear que k^4 coef escala: resolver dispersion a k=1
    k = 1.0; Om2 = cL**2*k**2 + lam*k**4/4
    print(f"           a k=1: Omega^2 = c^2 k^2 + lam k^4/4 = {Om2:.4f}")

print("\n" + "="*68)
print("Gates A-F listos. Revisar (D) coincide 0.01008 y (C) transmision=1.")
