"""
FASE 3B-0 — Infraestructura del scattering BdG (SIN matriz S todavia).

CORRECCION (ChatGPT, acertada): las ecuaciones BdG de la GPE son de 2o orden en
u y v -> sistema de PRIMER ORDEN de dimension 4, NO de orden 8. El k^4 de la
dispersion sale de elevar al cuadrado el k^2 en el DETERMINANTE, no de derivadas
de 4o orden.

BdG en fondo con flujo (hbar=m=1), modos con la fase del fondo factorizada:
  w u = [-1/2 (d_x + i v)^2 + 2 g n + V - mu] u + g n v
 -w v = [-1/2 (d_x - i v)^2 + 2 g n + V - mu] v + g n u

Expandidas (n,v,V funciones de x):
  u'' = -2 i v u' + (-i v' + v^2 + 4 g n + 2V - 2 mu - 2 w) u + 2 g n v
  v'' = +2 i v v' + (+i v' + v^2 + 4 g n + 2V - 2 mu + 2 w) v + 2 g n u

Sistema Y=(u,u',v,v')^T,  dY/dx = A(x,w) Y,  A es 4x4.

Validacion central: en region homogenea, los autovalores de A son {i k} con k las
4 raices del cuartico (w - v k)^2 = c_s^2 k^2 + k^4/4 de la fase 3A.
"""
import numpy as np
import scipy.linalg as la

g = 1.0
d = np.load("out/waterfall_state.npz")
x, n0, v0, Varr = d["x"], d["n"], d["v"], d["V"]
mu = float(d["mu"])
nL, vL, VL = n0[:50].mean(), v0[:50].mean(), Varr[:50].mean()
nR, vR, VR = n0[-50:].mean(), v0[-50:].mean(), Varr[-50:].mean()
cL, cR = np.sqrt(g*nL), np.sqrt(g*nR)

print("="*68)
print("FASE 3B-0 — sistema de primer orden BdG (dimension 4)")
print("="*68)
print(f"  Upstream : n={nL:.4f} v={vL:.4f} V={VL:.4f} c={cL:.4f}")
print(f"  Downstream: n={nR:.4f} v={vR:.4f} V={VR:.4f} c={cR:.4f}")
print(f"  mu={mu:.4f}")


def A_matrix(n, v, V, dv, omega):
    """Matriz 4x4 A(x,w) para Y=(u,u',v,v'). dv = v'(x)."""
    a_uu = -1j*dv + v**2 + 4*g*n + 2*V - 2*mu - 2*omega   # coef de u en u''
    a_uv = 2*g*n                                           # coef de v en u''
    a_uup = -2j*v                                          # coef de u' en u''
    b_vv = 1j*dv + v**2 + 4*g*n + 2*V - 2*mu + 2*omega     # coef de v en v''
    b_vu = 2*g*n                                           # coef de u en v''
    b_vvp = 2j*v                                           # coef de v' en v''
    return np.array([
        [0,      1,     0,     0    ],
        [a_uu,   a_uup, a_uv,  0    ],
        [0,      0,     0,     1    ],
        [b_vu,   0,     b_vv,  b_vvp]], dtype=complex)


def quartic_roots(omega, v, c):
    return np.roots([0.25, 0.0, c**2 - v**2, 2*omega*v, -omega**2])


# ---- VALIDACION 1: autovalores de A (homogeneo) vs raices del cuartico -------
print("\n" + "-"*68)
print("VALIDACION 1 — char.poly de A (4x4) == cuartico de dispersion")
print("-"*68)
maxerr = 0.0
for name, n, v, V, c in (("upstream", nL, vL, VL, cL), ("downstream", nR, vR, VR, cR)):
    for omega in (0.1, 0.2, 0.4):
        A = A_matrix(n, v, V, 0.0, omega)         # dv=0 en region homogenea
        k_fromA = -1j*la.eigvals(A)               # lambda=i k -> k=-i lambda
        k_quart = quartic_roots(omega, v, c)
        # emparejar cada k de A con su mas cercano del cuartico (robusto al orden)
        for kk in k_fromA:
            maxerr = max(maxerr, np.min(np.abs(k_quart - kk)))
print(f"  dim(A) = 4x4  (NO 8x8)  -> correccion confirmada")
print(f"  max |k(autovalores de A) - k(cuartico)| = {maxerr:.2e}")
print(f"  -> {'OK: sistema de orden 4, reproduce la dispersion' if maxerr<1e-10 else 'REVISAR'}")

# ---- Bases asintoticas: clasificar los 4 modos por lado ----------------------
def classify_modes(n, v, V, c, omega, side):
    A = A_matrix(n, v, V, 0.0, omega)
    lam, vec = la.eig(A)
    k = -1j*lam
    modes = []
    for i in range(4):
        kr = k[i]
        prop = abs(kr.imag) < 1e-9
        if prop:
            Om = omega - v*kr.real
            vg = v + (c**2*kr.real + kr.real**3/2)/Om
            # in/out: upstream(x<0): vg>0 -> hacia horizonte = IN; downstream(x>0): vg>0 = OUT
            if side == 'L':
                ch = 'IN' if vg > 0 else 'OUT'
            else:
                ch = 'OUT' if vg > 0 else 'IN'
            modes.append(dict(k=kr.real, prop=True, norm=int(np.sign(Om)),
                              vg=vg, ch=ch, Y=vec[:, i]))
        else:
            # evanescente: decae hacia el interior segun signo de Im(k)
            modes.append(dict(k=kr, prop=False, norm=0, vg=np.nan,
                              ch='evan', Y=vec[:, i]))
    return modes


print("\n" + "-"*68)
print("BASES ASINTOTICAS y ORDEN DE CANALES (omega=0.1)")
print("-"*68)
omega = 0.1
chan = []
for side, (n, v, V, c), lab in (('L', (nL, vL, VL, cL), 'UPSTREAM(sub)'),
                                 ('R', (nR, vR, VR, cR), 'DOWNSTREAM(sup)')):
    ms = classify_modes(n, v, V, c, omega, side)
    nprop = sum(m['prop'] for m in ms)
    print(f"  {lab}: {nprop} propagantes + {4-nprop} evanescentes")
    for m in ms:
        if m['prop']:
            sgn = '+' if m['norm'] > 0 else '-'
            print(f"      k={m['k']:+.4f} vg={m['vg']:+.4f} norma[{sgn}] {m['ch']}")
            chan.append((lab, m['ch'], sgn, m['k'], m['vg']))
        else:
            print(f"      k={m['k'].real:+.3f}{m['k'].imag:+.3f}i  evanescente")

ins = [c for c in chan if c[1] == 'IN']
outs = [c for c in chan if c[1] == 'OUT']
print(f"\n  Canales IN ({len(ins)}): " + ", ".join(f"{c[0][:4]}[{c[2]}]" for c in ins))
print(f"  Canales OUT({len(outs)}): " + ", ".join(f"{c[0][:4]}[{c[2]}]" for c in outs))
print(f"  -> matriz S sera {len(outs)}x{len(ins)};  eta = diag de signos de norma OUT")
etaout = [c[2] for c in outs]
print(f"  eta(OUT) = {etaout}  (hay norma '-': S pseudo-unitaria, |beta|^2>0 = Hawking)")

# ---- VALIDACION 2: la propagacion reproduce la solucion homogenea -----------
print("\n" + "-"*68)
print("VALIDACION 2 — propagar un modo plano homogeneo (RK4) conserva la forma")
print("-"*68)
# tomar un modo propagante upstream, propagarlo con A homogenea y comparar e^{ikx}
m0 = next(m for m in classify_modes(nL, vL, VL, cL, omega, 'L') if m['prop'])
Y = m0['Y'].astype(complex); k0 = m0['k']
dx = 0.01; N = 500
A0 = A_matrix(nL, vL, VL, 0.0, omega)
Yx = Y.copy()
for _ in range(N):
    k1 = A0@Yx; k2 = A0@(Yx+0.5*dx*k1); k3 = A0@(Yx+0.5*dx*k2); k4 = A0@(Yx+dx*k3)
    Yx = Yx + (dx/6)*(k1+2*k2+2*k3+k4)
expected = Y * np.exp(1j*k0*N*dx)        # solucion exacta e^{ikx}
err = np.max(np.abs(Yx - expected)) / np.max(np.abs(expected))
print(f"  propagado {N*dx:.1f} en x;  error vs e^(ikx) exacta = {err:.2e}")
print(f"  -> {'OK: la base y el propagador son consistentes' if err<1e-6 else 'REVISAR'}")

print("\n" + "="*68)
print("FASE 3B-0 COMPLETADA. Siguiente (sesion dedicada): matriz S + S^dag eta S = eta")
print("="*68)
