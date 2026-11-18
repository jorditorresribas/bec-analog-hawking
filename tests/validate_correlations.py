"""
VALIDACION de la maquinaria nueva (correlations.py) ANTES de usarla.
Gates (regla: verificar tras cualquier cambio):
 (A) dispersion lam=1 == quartic original (fd_bdg_operator.quartic_roots)
 (B) plane_mode_lam: |U|^2-|V|^2 = sign(Omega)  (norma de Bogoliubov, +1 o -1)
 (C) homogeneo: inyectar 'u' -> transmision=1, reflexion~0 (sin horizonte)
 (D) HORIZONTE lam=1: beta^2 (amplitud saliente norma- al inyectar 'u') debe
     coincidir con el valor VALIDADO 0.01008 (mode_matching/transfer)
 (E) pseudo-unitaridad por canal: sum_out eta|a|^2 = eta_in
 (F) knob lam: las raices propagantes cumplen (omega-vk)^2 = c^2 k^2 + lam k^4/4

Auditoria fase 1 (2026-11-17): cada gate tiene ahora una tolerancia y el script
sale con codigo 1 si alguna falla (antes solo imprimia y salia siempre con 0).
Correccion de (B): la version anterior comparaba |U|^2-|V|^2 con +1 para todos los
modos y por eso imprimia "max desviacion = 2.00": 3 de los 4 k de prueba eligen la
rama de norma NEGATIVA, cuyo valor correcto es -1 (= sign(Omega), Omega = omega - v k).
No era un fallo del codigo sino del test. El "2e-16" de docs/CORRELATIONS.md §2 no
correspondia a esta comprobacion.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import correlations as C
from fd_bdg_operator import quartic_roots

FAILS = []


def check(name, ok, detail):
    print(f"    [{'PASA' if ok else 'FALLA'}] {name}: {detail}")
    if not ok:
        FAILS.append(name)


print("="*68)
print("VALIDACION correlations.py")
print("="*68)

# (A) dispersion lam=1
om, vtest, ctest = 0.1, 0.5, 1.0
r_new = np.sort_complex(C.quartic_roots_lam(om, vtest, ctest, 1.0))
r_old = np.sort_complex(quartic_roots(om, vtest, ctest))
dA = np.max(np.abs(r_new-r_old))
print(f"\n(A) raices dispersion lam=1 vs original")
check("A", dA < 1e-12, f"max dif = {dA:.2e} (tol 1e-12)")

# (B) normalizacion de Bogoliubov: |U|^2-|V|^2 = sign(Omega)
print(f"\n(B) |U|^2-|V|^2 = sign(Omega) (n=0.7, v=0.4, lam=1)")
maxd = 0.0
for k in (0.3, 0.8, 1.5, -0.6):
    U, V = C.plane_mode_lam(k, 0.7, 0.4, 1.0)
    ev = np.linalg.eigvals(C._symbol_lam(k, 0.7, 0.4, 1.0, 0.0))
    e = ev[np.argmin(np.abs(ev))]               # misma eleccion que plane_mode_lam
    Om = (e - 0.4*k).real
    nb = abs(U)**2 - abs(V)**2
    maxd = max(maxd, abs(nb - np.sign(Om)))
    print(f"    k={k:+.1f}: Omega={Om:+.4f}  |U|^2-|V|^2={nb:+.15f}")
check("B", maxd < 1e-12, f"max |(|U|^2-|V|^2) - sign(Omega)| = {maxd:.2e} (tol 1e-12)")

# (C) homogeneo (sin horizonte): transmision=1, reflexion~0
print(f"\n(C) homogeneo (sin horizonte)")
r = C.scatter_field(L=20, N=3000, omega=0.1, channel='u', mode='homog')
tr = abs(r['out_amp'].get(('R', 1), 0))**2
rf = sum(abs(a)**2 for (side, nm), a in r['out_amp'].items() if side == 'L')
check("C", abs(tr-1) < 1e-6 and rf < 1e-10,
      f"transmision={tr:.10f} (tol |T-1|<1e-6)  reflexion={rf:.2e} (tol 1e-10)")

# (D) HORIZONTE: beta^2
print(f"\n(D) HORIZONTE lam=1, omega=0.1: beta^2 (inyecta 'u', lee saliente norma-)")
REF_B2 = 0.01008     # valor validado documentado (5 decimales -> +-0.5e-5)
for N in (3000, 4000, 6000):
    r = C.scatter_field(L=24, N=N, omega=0.1, channel='u', mode='horizon')
    beta2 = abs(r['out_amp'].get(('R', -1), 0))**2
    check(f"D N={N}", abs(beta2-REF_B2) < 5e-6 and abs(r['flux_out']-r['eta_in']) < 1e-3,
          f"beta^2={beta2:.7f} (ref {REF_B2} +-5e-6)  flujo_out={r['flux_out']:+.6f}")

# (E) pseudo-unitaridad por canal inyectado
print(f"\n(E) pseudo-unitaridad por canal (horizonte, omega=0.1, N=4000)")
TOL_E = {'u': 1e-6, 'd1': 5e-4, 'd2': 5e-4}   # documentado: 2.9e-7 / 1.2e-4 / 1.2e-4
for ch in ('u', 'd1', 'd2'):
    r = C.scatter_field(L=24, N=4000, omega=0.1, channel=ch, mode='horizon')
    dE = abs(r['flux_out']-r['eta_in'])
    check(f"E {ch}", dE < TOL_E[ch], f"eta_in={r['eta_in']:+d}  sum_out eta|a|^2={r['flux_out']:+.6f}"
          f"  defecto={dE:.2e} (tol {TOL_E[ch]:.0e})")

# (F) knob lam: las raices propagantes cumplen la dispersion modificada con c^2=gn
print(f"\n(F) knob lam: (omega-vk)^2 = c^2 k^2 + lam k^4/4 con c^2=gn fijo")
nL = 1.0; vL = 0.5; cL = np.sqrt(nL); omF = 0.1
maxF = 0.0
for lam in (0.5, 1.0, 2.0):
    for k in C.quartic_roots_lam(omF, vL, cL, lam):
        if abs(k.imag) < 1e-9:
            kr = k.real
            maxF = max(maxF, abs((omF - vL*kr)**2 - (cL**2*kr**2 + lam*kr**4/4)))
check("F", maxF < 1e-12, f"max residuo de la dispersion = {maxF:.2e} (tol 1e-12)")
print("    (nota: esto valida la DISPERSION homogenea; no dice nada del fondo. El operador")
print("     build_op_lambda con lam != 1 tiene el defecto de Qp documentado en CORRELATIONS.md.)")

print("\n" + "="*68)
if FAILS:
    print(f"FALLAN {len(FAILS)} gates: {', '.join(FAILS)}")
    _sys.exit(1)
print("TODOS LOS GATES PASAN")
_sys.exit(0)
