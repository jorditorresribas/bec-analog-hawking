"""
BLOQUE 2 — Gate B2.3: BARRIDO FINO de CONVERGENCIA de S (caja de 2o orden, X0 FIJO).

Caracteriza la convergencia en dx de la matriz S del BVP (representacion fisica), manteniendo:
  - X0 FIJO (el eje X0 = gate aparte; aqui NO se toca),
  - el MISMO esquema de caja (no se cambia esquema y resolucion a la vez),
  - ventana 3x3 (guard de B2.1).

Mide: orden observado por Richardson (ratio dx=2, por elemento y Frobenius), construye una referencia
S_ref por extrapolacion de Richardson orden 2, y reporta Delta_S(N) vs S_ref (por elemento + Frobenius)
y el dx/N necesario para alcanzar precisiones objetivo.

OJO (alcance): converge S_box(X0=12, dx->0), que NO es la S fisica (X0->inf). Eso es el gate de X0.
NO interpreta fisica, NO pseudo-unitaridad. Reutiliza por import el ensamblado auditado de B2.1 y la
representacion fisica de B2.2 (no los modifica).
"""
import numpy as np
import scipy.linalg as la
import gate_B2_1_algebraic as g1
import gate_B2_2_conditioning as g2


def S_box(w, N, X0=12.0):
    """S (3x3, convencion fisica) del BVP de caja en malla N. Reusa B2.1 (assemble) + B2.2 (fisica)."""
    M, info = g1.assemble(w, N, X0)
    labL, labR = info["labL"], info["labR"]
    window = (labL.count("in")==1 and labL.count("out")==1 and labL.count("dec")==1 and
              labL.count("grow")==1 and labR.count("in")==2 and labR.count("out")==2 and
              labR.count("dec")==0 and labR.count("grow")==0)
    if not window:
        raise RuntimeError(f"w={w} fuera de ventana 3x3")
    dr, dc = g2.physical_scaling(info)
    S, _, _, _ = g2.smatrix_and_solves(info, M, dr, dc)
    return S


def richardson_order(S_h, S_h2, S_h4):
    """Orden observado por elemento y Frobenius, con triple (h, h/2, h/4), ratio dx=2.
       p = log2( |S_h - S_h2| / |S_h2 - S_h4| )."""
    d1 = np.abs(S_h - S_h2); d2 = np.abs(S_h2 - S_h4)
    # orden por elemento: indefinido si d1==0 o d2==0 (NO es orden infinito) -> nan
    good = (d1 > 0) & (d2 > 0)
    p_elem = np.full(d1.shape, np.nan, dtype=float)
    p_elem[good] = np.log2(d1[good]/d2[good])
    n1 = la.norm(S_h - S_h2); n2 = la.norm(S_h2 - S_h4)
    pF = np.log2(n1/n2) if (n1 > 0 and n2 > 0) else np.nan   # indefinido si diferencia nula
    return p_elem, pF


def richardson_extrap(S_h, S_h2, p=2.0, r=2.0):
    """Extrapolacion de Richardson a dx->0 asumiendo orden p y ratio r: S_ref=(r^p S_h2 - S_h)/(r^p-1)."""
    f = r**p
    return (f*S_h2 - S_h)/(f - 1.0)


def convergence(w, X0=12.0, base_triples=((100,200,400),(150,300,600),(200,400,800)),
                seq=(100,150,200,300,400,600,800,1200)):
    print("="*86)
    print(f"GATE B2.3 — CONVERGENCIA de S   w={w}  X0={X0}  (caja 2o orden, X0 fijo)")
    print("="*86)
    Sc = {N: S_box(w, N, X0) for N in sorted(set(seq) | {n for t in base_triples for n in t})}

    # 1) orden observado por Richardson (ratio dx=2)
    print("\n  [1] Orden observado (Richardson, ratio dx=2):")
    pF_list = []
    for (a,b,c) in base_triples:
        p_elem, pF = richardson_order(Sc[a], Sc[b], Sc[c])
        pF_list.append(pF)
        finite = p_elem[np.isfinite(p_elem)]
        print(f"      N=({a:4d},{b:4d},{c:4d}) dx=({2*X0/a:.4f},{2*X0/b:.4f},{2*X0/c:.4f})  "
              f"p_Frob={pF:.3f}  p_elem in [{finite.min():.2f},{finite.max():.2f}]")
    print(f"      -> compatible con O(dx^2); p_Frob -> {pF_list[-1]:.3f} en el triple mas fino")

    # 2) referencia por extrapolacion de Richardson (orden 2) con el par ratio-2 mas fino disponible
    pairs = sorted((h, 2*h) for h in seq if 2*h in Sc)
    if len(pairs) < 2:
        raise ValueError("Se necesitan al menos dos pares de mallas con ratio dx=2.")
    h, h2 = pairs[-1]
    S_ref = richardson_extrap(Sc[h], Sc[h2], p=2.0, r=2.0)
    # cross-check: referencia desde el par ratio-2 inmediatamente mas grueso
    hC, h2C = pairs[-2]
    S_ref2 = richardson_extrap(Sc[hC], Sc[h2C], p=2.0, r=2.0)
    selfcons = float(np.max(np.abs(S_ref - S_ref2)))   # autoconsistencia entre extrapolaciones
    print(f"\n  [2] Referencia S_ref = Richardson(orden 2) de (N={h},{h2}); "
          f"cross-check vs (N={hC},{h2C}): max|dif| = {selfcons:.2e}  "
          f"(autoconsistencia, NO cota al continuo)")

    # 3) Delta_S(N) vs S_ref (Frobenius relativo + max abs por elemento)
    den = max(1.0, la.norm(S_ref))
    print("\n  [3] Delta_S(N) vs S_ref (extrapolada):")
    print(f"      {'N':>5} {'dx':>8} {'DeltaS_Frob':>12} {'max|elem|':>11}")
    rows = []
    for N in seq:
        dF = la.norm(Sc[N]-S_ref)/den
        dM = float(np.max(np.abs(Sc[N]-S_ref)))
        rows.append((N, 2.0*X0/N, dF, dM))
        print(f"      {N:5d} {2.0*X0/N:8.4f} {dF:12.3e} {dM:11.3e}")

    # 4) N/dx para alcanzar objetivos (medido si se alcanza; extrapolado O(dx^2) si no)
    print("\n  [4] N/dx para precision objetivo (Frobenius):")
    Ns = np.array([r[0] for r in rows]); dxs = 2.0*X0/Ns
    dF = np.array([r[2] for r in rows])
    # constante del ajuste dF ~ A dx^2 usando EXPLICITAMENTE las 3 mallas mas finas (menor dx)
    fine = np.argsort(dxs)[:3]
    if len(fine) < 3 or not np.all(np.isfinite(dF[fine]) & (dF[fine] > 0)):
        raise ValueError("No hay tres valores Delta_S positivos y finitos para ajustar A.")
    A = np.exp(np.mean(np.log(dF[fine]) - 2*np.log(dxs[fine])))
    Nfit = sorted(int(n) for n in Ns[fine])
    print(f"      (A estimada con las 3 mallas mas finas: N={Nfit}, "
          f"dx={[round(2.0*X0/n,4) for n in Nfit]})")
    for tol in (1e-2, 1e-3, 1e-4):
        reached = Ns[dF <= tol]
        if reached.size:
            print(f"      tol={tol:.0e}: ALCANZADO en N={int(reached.min())} (dx={2.0*X0/reached.min():.4f})")
        else:
            dx_need = np.sqrt(tol/A); N_need = 2.0*X0/dx_need
            flag = ("  [<< por debajo de la autoconsistencia de S_ref "
                    f"({selfcons:.1e}): NO respaldado]") if tol < selfcons else ""
            print(f"      tol={tol:.0e}: EXTRAPOLADO (O(dx^2), A={A:.2e}) -> "
                  f"dx~{dx_need:.4f}, N~{N_need:.0f}{flag}")
    print("\n  NOTA: S_box(X0=12). Converge la S del dominio TRUNCADO, NO la S fisica (X0->inf = gate aparte).")
    return dict(S=Sc, S_ref=S_ref, rows=rows, A=A)


if __name__ == "__main__":
    convergence(w=0.10)
    print("\n")
    convergence(w=0.05)
