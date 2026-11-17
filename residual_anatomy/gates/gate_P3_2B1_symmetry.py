"""
GATE P3.2-B1 (Control 0) — SIMETRIA por operador y por discretizacion.
ALCANCE (ChatGPT, diseño v2 aprobado): mide la autoadjuncion euclidea y PONDERADA (peso trapezoidal W)
de los bloques constituyentes D, A, K_P, C POR SEPARADO y por discretizacion, para LOCALIZAR defectos.
NO demuestra pseudo-herm de M_eff (eso es B3). NO valida scattering. NO toca v1.0.

NORMAS (DEFECTO RELATIVO EN NORMA MAXIMO-DE-ENTRADA, no una norma matricial inducida):
  defecto euclideo  abs = max_ij|X - X^dagger|_ij ;  rel = abs / max_ij|X_ij|
  defecto ponderado abs = max_ij|W X - X^dagger W|_ij ; rel = abs / max_ij|(W X)_ij| ; W=diag(1/2,1,...,1,1/2)
DISCRETIZACIONES:
  'plana'  : 2G + D_plano  (SIN filas de BC; truncamiento finito de la recta, NO = Neumann demostrado)
  'robin'  : 2G + D con Robin de 1er orden en bordes (P3.1)   -> W-simetrica, NO euclidea
  'dir-int': Dirichlet (filas identidad) RESTRINGIDO al interior (nodos 1..N-2)
HALLAZGO ESPERADO (verificado en diseño): plana -> K_P euclideo-s.a.; robin -> A es W-s.a. pero K_P NO.
"""
import numpy as np
import sys, os, contextlib
from pathlib import Path
root = Path(__file__).resolve().parents[2]
with open(os.devnull, 'w') as _dn, contextlib.redirect_stdout(_dn):   # descriptor cerrado explicitamente
    sys.path.insert(0, str(root / 'solver')); from correlations import bg, g

lam = 2.0; c1 = (lam-1)/4; c2 = -(lam-1)**2/16

def build_D_A(Gx, N, dx, disc):
    """Devuelve D (plano, -d^2, SIN BC), D2, y A=(2G+D_variant)^-1 segun la discretizacion."""
    main = 2/dx**2*np.ones(N); off = -1/dx**2*np.ones(N-1)
    D = np.diag(main) + np.diag(off, 1) + np.diag(off, -1)     # D PLANO (comun a todas: factor externo)
    D2 = D @ D
    twoGpD = 2*np.diag(Gx) + D.copy()                          # base
    if disc == 'robin':
        kR = np.sqrt(2*Gx[-1]); kL = np.sqrt(2*Gx[0])
        twoGpD[0, 0] = 2*Gx[0] + 2/dx**2 + 2*kL/dx; twoGpD[0, 1] = -2/dx**2
        twoGpD[-1, -1] = 2*Gx[-1] + 2/dx**2 + 2*kR/dx; twoGpD[-1, -2] = -2/dx**2
    elif disc in ('dirichlet', 'dir-int'):
        twoGpD[0, :] = 0; twoGpD[0, 0] = 1; twoGpD[-1, :] = 0; twoGpD[-1, -1] = 1
    A = np.linalg.inv(twoGpD)
    return D, D2, A

def C_operator(v, dx, form='mbdg'):
    """C = -i (v d_x + v'/2). form='mbdg': diag(v)@D1 + 1/2 diag(v') (COMO build_MBdG v1.0;
    hermitica solo a O(dx) en norma relativa de operador). form='exact': -i/2 (vD1+D1v)
    (antisimetrica EXACTA -> C hermitica a maquina), pero NO coincide con v1.0."""
    N = len(v); e = np.ones(N-1)
    D1 = (np.diag(e, 1) - np.diag(e, -1))/(2*dx)
    if form == 'exact':
        return -0.5j*(np.diag(v)@D1 + D1@np.diag(v))
    vp = np.gradient(v, dx)
    return -1j*(np.diag(v) @ D1 + 0.5*np.diag(vp))

def defects(X, W):
    eA = np.abs(X - X.conj().T).max(); eR = eA/np.abs(X).max()
    Wm = np.diag(W)
    wA = np.abs(Wm@X - X.conj().T@Wm).max(); wR = wA/np.abs(Wm@X).max()
    return eA, eR, wA, wR

def run():
    print("GATE P3.2-B1 (Control 0) — simetria por operador y discretizacion\n" + "="*72)
    N = 160; L = 16.; x = np.linspace(-L, L, N); dx = x[1]-x[0]
    n, v, V = bg(x, w=2.0, mode='horizon'); Gx = g*n
    W = np.ones(N); W[0] = 0.5; W[-1] = 0.5
    C = C_operator(v, dx)
    print(f"fondo tanh w=2, N={N}, L={L}, dx={dx:.4f}; norma: max|.|; rel = abs/max|.|")
    print(f"{'disc':>8} {'op':>4} {'eucl_abs':>10} {'eucl_rel':>9} {'W_abs':>10} {'W_rel':>9}")
    for disc in ['plana', 'robin', 'dir-int']:
        D, D2, A = build_D_A(Gx, N, dx, disc)
        KP = 0.5*D + c1*(D@A@D) + c2*(D2@(A@A@A)@D2)
        ops = {'D': D, 'A': A, 'K_P': KP, 'C': C}
        for nm, X in ops.items():
            if disc == 'dir-int' and nm in ('A', 'K_P'):      # Dirichlet: medir en INTERIOR
                sl = slice(1, N-1); Xi = X[sl, sl]; Wi = W[sl]
                eA, eR, wA, wR = defects(Xi, Wi); tag = nm+'(int)'
            else:
                eA, eR, wA, wR = defects(X, W); tag = nm
            print(f"{disc:>8} {tag:>4} {eA:10.2e} {eR:9.2e} {wA:10.2e} {wR:9.2e}")
        print()
    # EXPERIMENTO SEPARADO (homog G=1, NO continuacion de la tabla tanh de arriba):
    # crecimiento del defecto numerico de K_P plana con N. Causa PROBABLE (no demostrada en
    # exclusiva por este barrido): condicionamiento de formar D2A3D2 densa (~1/dx^4).
    print("Experimento SEPARADO (homog G=1) — crecimiento numerico de K_P plana eucl_rel con N:")
    for Nt in [120, 240, 480]:
        xt = np.linspace(-L, L, Nt); dxt = xt[1]-xt[0]; Gt = np.ones(Nt)
        D, D2, A = build_D_A(Gt, Nt, dxt, 'plana')
        KP = 0.5*D + c1*(D@A@D) + c2*(D2@(A@A@A)@D2)
        print(f"  N={Nt:4d} dx={dxt:.4f} K_P eucl_rel={np.abs(KP-KP.conj().T).max()/np.abs(KP).max():.2e}")
    # convergencia de C: hermiticidad O(dx) (forma mbdg) vs exacta
    print("\nConvergencia de C (hermiticidad): forma 'mbdg' (v1.0) es O(dx); 'exact' es ~maquina:")
    for Nt in [160, 320, 640]:
        xt = np.linspace(-L, L, Nt); dxt = xt[1]-xt[0]; _, vt, _ = bg(xt, w=2.0, mode='horizon')
        Cm = C_operator(vt, dxt, 'mbdg'); Ce = C_operator(vt, dxt, 'exact')
        rm = np.abs(Cm-Cm.conj().T).max()/np.abs(Cm).max()
        re = np.abs(Ce-Ce.conj().T).max()/np.abs(Ce).max()
        print(f"  N={Nt:4d} dx={dxt:.4f}  C_mbdg eucl_rel={rm:.2e} (O(dx))  C_exact eucl_rel={re:.1e}")
    print("\nLectura (CORREGIDA):")
    print(" - plana: K_P euclideo-s.a.; el 'defecto' crece con N = PROBABLE condicionamiento de D2A3D2")
    print("   DENSA (1/dx^4), no demostrado en exclusiva aqui; en banda disperso se evita.")
    print(" - robin: A es W-s.a. (W_rel~1e-16) pero K_P NO (ni eucl ni W).")
    print(" - dir-int: interior(K_P) NO simetrico = acoplo interior-borde (interior(prod)!=prod(interior)).")
    print(" - C (forma build_MBdG v1.0): hermitico solo a O(dx) en norma rel de operador, NO de borde,")
    print("   NO O(dx^2). Converge linealmente. Forma exacta 1/2(vD1+D1v) seria maquina pero != v1.0.")
    print("ALCANCE: localiza simetrias; NO demuestra pseudo-herm de M_eff (B3).")

if __name__ == '__main__':
    run()
