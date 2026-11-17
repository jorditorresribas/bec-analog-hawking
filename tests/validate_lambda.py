"""
VALIDACION del operador modificado build_op_lambda (knob k^4) ANTES de usarlo.
Gates:
 (A) homogeneo op='lambda': transmision=1, reflexion~0 (cualquier lam)
 (B) lam=1 op='lambda' vs op='validated': beta^2 debe coincidir salvo O(Qp)
 (C) pseudo-unitaridad por canal (op='lambda', varios lam)
 (D) escalado: beta^2(lam) y ventana de Hawking omega_max(lam) con c^2=gn FIJO
 (E) metrica fija: velocidad del sonido asintotica identica para todo lam
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import correlations as C

print("="*68); print("VALIDACION knob k^4 (build_op_lambda)"); print("="*68)

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
    print(f"    lam={lam}: transmision={tr:.4f}  reflexion={rf:.2e}")

# (B) lam=1 modificado vs validado
print("\n(B) lam=1: beta^2 op='lambda' (sin Qp) vs op='validated' (con Qp)")
for N in (4000, 6000):
    b_l, _, _ = beta2(1.0, N, 'lambda')
    b_v, _, _ = beta2(1.0, N, 'validated')
    print(f"    N={N}: beta^2(lambda)={b_l:.5f}  beta^2(validado)={b_v:.5f}  "
          f"dif rel={abs(b_l-b_v)/b_v:.1%}  (= efecto Qp)")

# (C) pseudo-unitaridad
print("\n(C) pseudo-unitaridad por canal (op='lambda', N=4000, omega=0.1):")
for lam in (0.5, 1.0, 2.0):
    defs = []
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L=24, N=4000, omega=0.1, channel=ch, lam=lam,
                            mode='horizon', op='lambda')
        if r['Rin_neg'] is not None:
            defs.append(abs(r['flux_out']-r['eta_in']))
    print(f"    lam={lam}: max defecto pseudo-unit = {max(defs):.2e}  "
          f"({len(defs)}/3 canales con d2)")

# (D) escalado beta^2(lam) y ventana
print("\n(D) beta^2(lam) a omega=0.1 (c^2=gn FIJO, solo cambia dispersion UV):")
for lam in (0.25, 0.5, 1.0, 2.0, 4.0):
    b, fo, ei = beta2(lam, 4000, 'lambda')
    print(f"    lam={lam:4}: beta^2={b:.5f}  flujo_out={fo:+.4f} (eta_in={ei:+d})")

# (E) ventana de Hawking: omega_max(lam) (donde el canal d2 desaparece)
print("\n(E) omega_max(lam): borde superior de la ventana de Hawking")
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
print("Revisar: (A) transmision=1, (B) dif Qp pequena y estable, (C) pseudo-unit<1e-2")
