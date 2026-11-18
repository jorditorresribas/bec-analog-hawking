"""
VALIDACION del operador modificado build_op_lambda (knob k^4) ANTES de usarlo.

!! AVISO (auditoria fase 1, 2026-11-17) !!
build_op_lambda mantiene la presion cuantica del fondo Qp INDEPENDIENTE de lam.
Para lam != 1 ese fondo NO es estacionario para el operador modificado (rompe el
modo de Goldstone; residuo 1/2|1-sqrt(lam)| ||(sqrt n)''|| medido con cociente 1.0000,
residual_anatomy/diag_goldstone.py). Los gates de abajo solo comprueban propiedades
MECANICAS (fondo homogeneo, lam=1, pseudo-unitariedad), que el defecto NO detecta:
por eso pasaban. El apartado (D) es la FIRMA DEL ARTEFACTO, no una validacion:
beta^2(lam) no monotona con minimo en lam=1. Ver docs/CORRELATIONS.md (aviso inicial
y §9) y residual_anatomy/diag_scaling_equivalence.py.

Gates con tolerancia (el script sale con codigo 1 si alguno falla):
 (A) homogeneo op='lambda': transmision=1, reflexion~0 (cualquier lam)
 (B) lam=1 op='lambda' vs op='validated': beta^2 identico (mismo operador en lam=1)
 (C) pseudo-unitaridad por canal (op='lambda', varios lam)
Informativos (sin tolerancia; NO son validacion fisica):
 (D) beta^2(lam) con el operador defectuoso (firma del artefacto)
 (E) ventana de Hawking omega_max(lam): depende solo de la dispersion homogenea,
     no del fondo -> no afectada por el defecto
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import correlations as C

FAILS = []


def check(name, ok, detail):
    print(f"    [{'PASA' if ok else 'FALLA'}] {name}: {detail}")
    if not ok:
        FAILS.append(name)


print("="*68); print("VALIDACION knob k^4 (build_op_lambda) — ver AVISO en la cabecera"); print("="*68)

def beta2(lam, N=4000, op='lambda', om=0.1):
    r = C.scatter_field(L=24, N=N, omega=om, channel='u', w=2.0, lam=lam,
                        mode='horizon', op=op)
    return abs(r['out_amp'].get(('R', -1), 0))**2, r['flux_out'], r['eta_in']

# (A) homogeneo
print("\n(A) homogeneo op='lambda': transmision=1, reflexion~0")
for lam in (0.5, 1.0, 2.0):
    r = C.scatter_field(L=20, N=3000, omega=0.1, channel='u', lam=lam,
                        mode='homog', op='lambda')
    tr = abs(r['out_amp'].get(('R', 1), 0))**2
    rf = sum(abs(a)**2 for (s, nm), a in r['out_amp'].items() if s == 'L')
    check(f"A lam={lam}", abs(tr-1) < 1e-6 and rf < 1e-10,
          f"transmision={tr:.10f} (tol 1e-6)  reflexion={rf:.2e} (tol 1e-10)")

# (B) lam=1 modificado vs validado (en lam=1 ambos operadores son el mismo)
print("\n(B) lam=1: beta^2 op='lambda' vs op='validated' (ambos con Qp; en lam=1 coinciden)")
for N in (4000, 6000):
    b_l, _, _ = beta2(1.0, N, 'lambda')
    b_v, _, _ = beta2(1.0, N, 'validated')
    d = abs(b_l-b_v)/b_v
    # Historial de la tolerancia (declarado): 1e-8 -> FALLO (dif real 3e-8: mismo operador, distinto orden
    # de operaciones + mal condicionamiento); 1e-5 (suelo de redondeo ~7e-6) -> demasiado generosa (~300x la dif
    # real); 1e-6 (2026-11-17): ~30x la dif real, sigue detectando fallos. Los dos ultimos pasos los propuso el
    # asesor del chat (Claude en claude.ai, misma familia que el implementador): NO es una revision independiente.
    # (Hasta el 2026-11-18 este comentario lo atribuia por error "al revisor".)
    check(f"B N={N}", d < 1e-6, f"beta^2(lambda)={b_l:.8f}  beta^2(validado)={b_v:.8f}  dif rel={d:.1e} (tol 1e-6)")

# (C) pseudo-unitaridad
print("\n(C) pseudo-unitaridad por canal (op='lambda', N=4000, omega=0.1)")
print("    (comprobacion MECANICA: el operador defectuoso tambien la cumple)")
for lam in (0.5, 1.0, 2.0):
    defs = []
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L=24, N=4000, omega=0.1, channel=ch, lam=lam,
                            mode='horizon', op='lambda')
        if r['Rin_neg'] is not None:
            defs.append(abs(r['flux_out']-r['eta_in']))
    check(f"C lam={lam}", len(defs) == 3 and max(defs) < 1e-3,
          f"max defecto pseudo-unit = {max(defs):.2e} (tol 1e-3)  ({len(defs)}/3 canales con d2)")

# (D) INFORMATIVO: firma del artefacto
print("\n(D) INFORMATIVO — beta^2(lam) a omega=0.1 con el operador DEFECTUOSO (Qp fijo).")
print("    No monotona con minimo en lam=1: es la firma del artefacto, NO fisica.")
for lam in (0.25, 0.5, 1.0, 2.0, 4.0):
    b, fo, ei = beta2(lam, 4000, 'lambda')
    print(f"    lam={lam:4}: beta^2={b:.5f}  flujo_out={fo:+.4f} (eta_in={ei:+d})")

# (E) INFORMATIVO: ventana de Hawking omega_max(lam) (donde el canal d2 desaparece)
print("\n(E) INFORMATIVO — omega_max(lam): borde superior de la ventana de Hawking")
def omega_max(lam):
    lo, hi = 0.05, 1.2
    for _ in range(30):
        mid = 0.5*(lo+hi)
        r = C.scatter_field(L=24, N=2500, omega=mid, channel='d2', lam=lam,
                            mode='horizon', op='lambda')
        if r['Rin_neg'] is not None and r['nchan'] == 2:
            lo = mid
        else:
            hi = mid
    return lo
for lam in (0.5, 1.0, 2.0):
    print(f"    lam={lam}: omega_max = {omega_max(lam):.3f}  (lam=1 validado: 0.468)")

print("\n" + "="*68)
if FAILS:
    print(f"FALLAN {len(FAILS)} gates: {', '.join(FAILS)}")
    _sys.exit(1)
print("GATES MECANICOS (A-C) PASAN — recuerda: no validan la fisica para lam != 1")
_sys.exit(0)
