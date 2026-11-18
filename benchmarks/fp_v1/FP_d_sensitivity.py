"""
SENSIBILIDAD del cociente T_upstream/(kbar/2pi) a la anchura de ventana d (auditoria fase 1, 2026-11-17).
Reutiliza T_upstream congelado (reference/fp_benchmark_reference.csv) y recalcula kbar(d) con la MISMA
geometria y ventana simetrica que FP_benchmark_table.geom (alli d = 1.58 fijo). No resuelve BdG (~1 s).
Comprueba que d = 1.58 reproduce la columna Tup_over_kbar2pi de la referencia (tol 1e-6); sale con 1 si no.
Tabla documentada en docs/FP_BENCHMARK.md §5b.
"""
import csv, os, sys
import numpy as np

ROOT = os.path.dirname(os.path.dirname(os.path.dirname(os.path.abspath(__file__))))
ref = {float(r['w']): r for r in csv.DictReader(open(os.path.join(ROOT, 'reference', 'fp_benchmark_reference.csv')))}
g, n_u, mach_u, md = 1.0, 1.0, 0.5, 2.0
jc = n_u*mach_u*np.sqrt(g*n_u); n_d = (jc/(md*np.sqrt(g)))**(2/3)
WS = (2.0, 3.0, 4.0, 6.0, 8.0)
DS = (1e-4, 1.0, 1.58, 2.0, 3.0)


def kbar(w, d):   # = FP_benchmark_table.geom con d variable
    L = max(40., 8*w+20.); x = np.linspace(-L, L, 80000); n = n_u+(n_d-n_u)*0.5*(1+np.tanh(x/w))
    c = np.sqrt(n); v = jc/n; f = c-v; i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    xH = x[i]+f[i]/(f[i]-f[i+1])*(x[i+1]-x[i]); cH = np.interp(xH, x, c)
    kK = abs(np.interp(xH, x, np.gradient(f, x)))
    dxi = d*cH/(kK*(2*cH**2)**2)**(1/3)
    return abs(np.interp(xH+0.5*dxi, x, f)-np.interp(xH-0.5*dxi, x, f))/dxi


print("T_upstream/(kbar(d)/2pi), ventana simetrica (d->0 equivale a kappa_local):")
print("  d     " + "  ".join(f"w={int(w)}  " for w in WS))
worst = 0.0
for d in DS:
    vals = [float(ref[w]['T_upstream'])/(kbar(w, d)/(2*np.pi)) for w in WS]
    print(f"  {d:<5g} " + "  ".join(f"{v:.4f}" for v in vals))
    if d == 1.58:
        worst = max(abs(v-float(ref[w]['Tup_over_kbar2pi'])) for v, w in zip(vals, WS))
ok = worst < 1e-6
print(f"control d=1.58 vs referencia: max dif = {worst:.1e} (tol 1e-6) -> {'PASA' if ok else 'FALLA'}")
sys.exit(0 if ok else 1)
