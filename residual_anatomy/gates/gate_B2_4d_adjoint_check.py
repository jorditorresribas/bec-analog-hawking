"""
BLOQUE 2 — B2.4d: verificacion INDEPENDIENTE del adjunto (sin tocar B2.2/B2.4 ni registros).

B2.4b (ruta antigua) dio ||y||=4.79e8 a X0=28; B2.4c (lu_solve trans=2) dio ||y||=108, r_adj=5e-14.
Aqui se reproducen VARIAS rutas sobre EXACTAMENTE la misma Mrep, RHS y selector c, en X0=12 y 28,
para localizar la causa (NO se asume que sea la no-contiguidad). Se comparan:
  - r_adj relativo y ABSOLUTO,
  - identidad  s = y^dag b  vs  s0 = c^dag x,
  - ||y|| de cada ruta y ||y_ruta - y_lu||.

Rutas:
  A  : la.solve(Mrep.conj().T, crep)                 (exactamente como g2.kappa_out / B2.4b)
  A2 : la.solve(ascontiguousarray(Mrep.conj().T), c) (test de la hipotesis de contiguidad)
  B  : lu_solve(lu_factor(Mrep), crep, trans=2)       (ruta validada en B2.4c)
NO modifica B2.1/2/3/v1.0. Proceso fresco y aislado.
"""
import numpy as np
import scipy.linalg as la
import gate_B2_1_algebraic as g1
import gate_B2_2_conditioning as g2


def check(w, X0, dx=0.02, out_sel=0, in_prob=0):
    N = int(round(2*X0/dx))
    print("="*92)
    print(f"[B2.4d] w={w} X0={X0} N={N} ndof={4*(N+1)}")
    print("="*92)
    M, info = g1.assemble(w, N, X0)
    dr, dc = g2.physical_scaling(info)
    Mrep = dr[:, None] * M * dc[None, :]; del M
    sin = np.zeros(3); sin[in_prob] = 1.0
    brep = dr * g1.rhs_for(info, sin)
    outs = g2.out_indices(info)
    c = np.zeros(info["ndof"], dtype=complex); c[outs[out_sel]] = 1.0
    crep = np.conj(dc) * c
    # primal (una LU reutilizable)
    lu = la.lu_factor(Mrep)
    z = la.lu_solve(lu, brep)
    s0 = np.vdot(crep, z)                 # s = c^dag x (x=dc z ; aqui en vars escaladas c^dag z)
    MH = Mrep.conj().T                    # vista no contigua
    def report(name, y):
        res = MH @ y - crep
        r_abs = la.norm(res); r_rel = r_abs/la.norm(crep)
        spred = np.vdot(y, brep)          # y^dag brep (debe == s0 si y correcto)
        print(f"  [{name:3s}] ||y||={la.norm(y):.4e}  r_adj_abs={r_abs:.2e}  r_adj_rel={r_rel:.2e}  "
              f"s=y^dag b={spred:+.4e}  (s0={s0:+.4e})  |s-s0|={abs(spred-s0):.2e}")
        return y
    yB  = report("B",  la.lu_solve(lu, crep, trans=2))
    yA  = report("A",  la.solve(MH, crep))
    yA2 = report("A2", la.solve(np.ascontiguousarray(MH), crep))
    print(f"  ||yA - yB|| = {la.norm(yA-yB):.3e}   ||yA2 - yB|| = {la.norm(yA2-yB):.3e}")
    print(f"  (B es la ruta validada por SVD en B2.4c. Si A difiere mucho de B, A es la defectuosa.)")


if __name__ == "__main__":
    check(0.10, X0=12, dx=0.02)
    check(0.10, X0=28, dx=0.02)
    print("\n>>> B2.4d COMPLETO")
