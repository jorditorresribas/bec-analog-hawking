"""
Ejecuta todas las validaciones con codigo de salida y resume PASA/FALLA.
Uso (desde la raiz del repositorio):  python3 tests/run_all.py
Sale con codigo 1 si alguna falla. Duracion aproximada: ~2 min.
(Auditoria fase 1, 2026-11-17: antes los tests no tenian aserciones y salian siempre con 0.)
"""
import os, subprocess, sys, time

ROOT = os.path.dirname(os.path.dirname(os.path.abspath(__file__)))
CHECKS = [
    'tests/validate.py',
    'tests/validate_correlations.py',
    'tests/validate_correlation_norm.py',
    'tests/validate_lambda.py',
    'benchmarks/fp_v1/FP_check_reproduction.py',
    'benchmarks/fp_v1/FP_d_sensitivity.py',
    'residual_anatomy/diag_scaling_equivalence.py',
]
results = []
for s in CHECKS:
    t0 = time.time()
    p = subprocess.run([sys.executable, s], cwd=ROOT, stdout=subprocess.PIPE, stderr=subprocess.STDOUT, text=True)
    results.append((s, p.returncode, time.time() - t0))
    print(f"[{'PASA' if p.returncode == 0 else 'FALLA'}] {s}  (exit {p.returncode}, {time.time()-t0:.0f} s)")
    if p.returncode != 0:
        print('\n'.join('    ' + l for l in p.stdout.strip().splitlines()[-15:]))
bad = [s for s, rc, _ in results if rc != 0]
print(f"\n{len(results)-len(bad)}/{len(results)} pasan" + (f" — FALLAN: {', '.join(bad)}" if bad else ""))
sys.exit(1 if bad else 0)
