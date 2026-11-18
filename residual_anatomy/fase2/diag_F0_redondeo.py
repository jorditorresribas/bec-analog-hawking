"""
DIAGNOSTICO A POSTERIORI de F0 (no cambia su veredicto). F0 fallo en (a) y (b): con gamma = 10 y N = 8000 el residuo
de Goldstone (5.2e-7) supera al de v1.0 (9.4e-8). Hipotesis: redondeo de la resta P sqrt(n) - P sqrt(n) en el producto
matriz-vector (entradas de P ~ gamma g / dx^4 ~ 1e9 -> error ~ 1e9 * eps). Prueba: el mismo producto M Phi_G con
acumulacion en long double (80 bits), y la norma del bloque P aplicado por separado a cada columna.
"""
import os, sys, numpy as np
sys.path.insert(0, os.path.dirname(os.path.abspath(__file__)))
import sextic as S
def matvec_ld(M, z):
    M = M.tocsr(); zl = z.astype(np.clongdouble); out = np.zeros(M.shape[0], dtype=np.clongdouble)
    rows = np.repeat(np.arange(M.shape[0]), np.diff(M.indptr)); np.add.at(out, rows, M.data.astype(np.clongdouble)*zl[M.indices])
    return out
print(" N    gamma | residuo double | residuo long double | v1.0 (long double) | ||gamma g P sqrt n|| (escala de la resta)")
for N in (2000, 4000, 8000, 16000):
    x = np.linspace(-40, 40, N); n, v, V0, n4 = S.background(x, 2.0); s = np.sqrt(n)
    Phi = np.concatenate([s, -s]).astype(complex); mk = np.abs(x) <= 30.0; mk2 = np.concatenate([mk, mk])
    nrm = np.linalg.norm(Phi[mk2])
    r0 = matvec_ld(S.operator(x, n, v, V0, n4, 0.0), Phi)
    for gam in (10.0,):
        M = S.operator(x, n, v, V0, n4, gam)
        rd = M @ Phi; rl = matvec_ld(M, Phi)
        Ponly = (M - S.operator(x, n, v, V0, n4, 0.0)) @ np.concatenate([s, 0*s]).astype(complex)
        print(f"{N:6d} {gam:5.0f} | {np.linalg.norm(rd[mk2])/nrm:.3e}      | {float(np.linalg.norm(rl[mk2].astype(complex)))/nrm:.3e}"
              f"           | {float(np.linalg.norm(r0[mk2].astype(complex)))/nrm:.3e}          | {np.linalg.norm(Ponly[mk2])/nrm:.2e}")
