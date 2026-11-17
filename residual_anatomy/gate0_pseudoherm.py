"""
GATE 0 (RUMBO #2) — operador no local discreto: pseudo-hermiticidad + consistencia/convergencia.

Ruta: cuantizacion de WEYL (pseudo-espectral exacta). Simbolo cinetico
    eps_lam(x,k) = -G(x) + sqrt(G(x)^2 + c(x)^2 k^2 + (lam/4) k^4),   G(x)=c(x)^2=g n(x).
Operador cinetico por nucleo de Weyl (ordenamiento simetrico, punto medio):
    K(x,y) = (1/2pi) int dk e^{ik(x-y)} eps_lam((x+y)/2, k)
           = (1/pi) int_0^{kmax} dk cos(k(x-y)) eps_lam((x+y)/2, k)   [eps par en k]
-> real y simetrico en (x,y) => cinetico HERMITICO por construccion. BdG:
    L_± = K ∓ i(v d + v'/2) + (G + Qp),   M = [[L+, G],[-G, -L-]].

NO modifica la entrega congelada ni build_op_lambda. Implementa SOLO el gate 0 de NONLOCAL_SPEC.md.

Dos partes (la 2a es la no trivial, ChatGPT):
  A. pseudo-hermiticidad M†η=ηM  (debe salir ~maquina por construccion simetrica de K).
  B. consistencia/convergencia: K reproduce el simbolo eps_lam(k) en fondo homogeneo;
     estudiar escalado del residuo y del error de dispersion con N, L (no un solo tamano).
"""
import numpy as np

g = 1.0; n_u = 1.0; mach_u = 0.5
vLb = mach_u*np.sqrt(g*n_u); jcur = n_u*vLb; mu = g*n_u + 0.5*vLb**2

def nval(x, w=2.0, mode='horizon'):
    if mode == 'homog':
        return np.full_like(np.asarray(x, float), n_u)
    md = 2.0; n_d = (jcur/(md*np.sqrt(g)))**(2/3)
    return n_u + (n_d - n_u)*0.5*(1+np.tanh(x/w))

def bg(x, w=2.0, mode='horizon'):
    if mode == 'homog':   # upstream CONSISTENTE: n=1,v=vLb,V=0 -> Qp=V-mu+gn+v^2/2=0 (ChatGPT)
        return np.full_like(x, n_u), np.full_like(x, vLb), np.zeros_like(x)
    md = 2.0; n_d = (jcur/(md*np.sqrt(g)))**(2/3)
    t = np.tanh(x/w); s = 0.5*(1+t); Dn = n_d-n_u
    n = n_u+Dn*s; n1 = Dn*(1-t**2)/(2*w); n2 = -Dn*t*(1-t**2)/w**2
    v = jcur/n; Qp = 0.5*(n2/(2*n)-n1**2/(4*n**2))
    return n, v, mu-g*n-0.5*v**2+Qp

def weyl_kinetic(x, lam, w=2.0, mode='horizon', Nk=None):
    """Nucleo cinetico de Weyl K(x_i,x_j) (N x N, real simetrico).
    Nk (cuadratura en k) debe escalar con N: el integrando cos(k*D)*eps oscila con
    k*D (D hasta ~2L), y resolverlo exige dk < pi/(2L) => Nk >~ N. Si no, el error de
    cuadratura CRECE con N (el kernel 'converge' y luego se degrada)."""
    N = len(x); dx = x[1]-x[0]; kmax = np.pi/dx
    if Nk is None:
        Nk = max(800, 3*N)
    kk = np.linspace(0.0, kmax, Nk); dk = kk[1]-kk[0]
    wq = np.ones(Nk); wq[0] = wq[-1] = 0.5                 # trapezoidal
    XM = 0.5*(x[:, None] + x[None, :])                     # punto medio (N,N)
    D = x[:, None] - x[None, :]                            # x_i - x_j
    nXM = nval(XM, w, mode); G = g*nXM; c2 = g*nXM         # G=c^2=g n(xm)
    K = np.zeros((N, N))
    for ik, k in enumerate(kk):
        eps = -G + np.sqrt(G**2 + c2*k*k + lam*k**4/4.0)
        K += (wq[ik]*dk) * np.cos(k*D) * eps
    # K(x_i,x_j) continuo = (1/pi) int_0^kmax ...; matvec discreto necesita la medida dx
    # de la suma sobre j:  (K u)_i = sum_j [K(x_i,x_j) dx] u_j
    return (K/np.pi) * dx

def build_nonlocal_M(x, lam, w=2.0, mode='horizon', Nk=None):
    N = len(x); dx = x[1]-x[0]
    n, v, V = bg(x, w, mode)
    K = weyl_kinetic(x, lam, w, mode, Nk)
    e = np.ones(N-1)
    D1 = (np.diag(e, 1) - np.diag(e, -1))/(2*dx)           # exactamente antisimetrico (skew)
    # conveccion simetrizada en forma EXACTAMENTE skew: 0.5(v D1 + D1 v) -> v d + v'/2 en el
    # continuo, pero skew-hermitica a maquina en discreto => M pseudo-hermitico EXACTO (residuo 0).
    conv = 0.5*(np.diag(v)@D1 + D1@np.diag(v))
    G = g*n; Qp = V - mu + g*n + 0.5*v**2                  # G=gn FIJO; Qp del fondo
    diagGQ = np.diag(G + Qp)
    Lp = K - 1j*conv + diagGQ
    Lm = K + 1j*conv + diagGQ
    Gd = np.diag(G)
    M = np.block([[Lp, Gd], [-Gd, -Lm]])
    return M

def pseudoherm_residual(M):
    N2 = M.shape[0]; N = N2//2
    eta = np.concatenate([np.ones(N), -np.ones(N)])
    R = M.conj().T @ np.diag(eta) - np.diag(eta) @ M       # M†η - ηM  (=0 si pseudo-herm)
    return np.abs(R).max()/np.abs(M).max()

# ============================================================
print("="*70); print("GATE 0 — operador no local (Weyl): pseudo-herm + consistencia"); print("="*70)

# --- A. Pseudo-hermiticidad (fondo con horizonte) ---
print("\n[A] Pseudo-hermiticidad M†η=ηM  (fondo horizonte, L=40)")
print(f"    {'lam':>5} {'N':>6} {'residuo rel':>13}")
for lam in (1.0, 2.0, 4.0):
    for N in (200, 400):
        x = np.linspace(-40, 40, N)
        M = build_nonlocal_M(x, lam, w=2.0)
        print(f"    {lam:5.1f} {N:6d} {pseudoherm_residual(M):13.2e}")

# --- B. Consistencia de dispersion en fondo HOMOGENEO ---
# K aplicado a onda plana e^{ikx} debe dar eps_lam(k) e^{ikx} en el interior.
print("\n[B] Consistencia: autovalor de K sobre onda plana vs simbolo eps_lam(k)  (homogeneo n=1)")
def eps_symbol(k, lam):
    G = g*n_u; c2 = g*n_u
    return -G + np.sqrt(G**2 + c2*k*k + lam*k**4/4.0)
print(f"    {'lam':>5} {'k':>7} {'K·u/u (interior)':>18} {'eps_lam(k)':>12} {'err rel':>10}")
for lam in (1.0, 2.0):
    x = np.linspace(-40, 40, 800);
    K = weyl_kinetic(x, lam, mode='homog')
    for kphys in (0.3, 1.0, 2.0):
        u = np.exp(1j*kphys*x)
        Ku = K @ u
        c = N0 = len(x)//2
        sl = slice(len(x)//2-60, len(x)//2+60)            # interior, lejos de bordes
        ratio = np.real(np.mean(Ku[sl]/u[sl]))
        es = eps_symbol(kphys, lam)
        print(f"    {lam:5.1f} {kphys:7.3f} {ratio:18.6f} {es:12.6f} {abs(ratio-es)/es:10.2e}")

# --- C. Convergencia del residuo de pseudo-herm con N y L (no un solo tamano) ---
print("\n[C] Convergencia pseudo-herm con N (L=40) y con L (N=400), lam=2  [horizonte]")
print(f"    {'N':>6} {'L':>5} {'residuo rel':>13}")
for N in (150, 300, 600):
    x = np.linspace(-40, 40, N); M = build_nonlocal_M(x, 2.0)
    print(f"    {N:6d} {40:5d} {pseudoherm_residual(M):13.2e}")
for L in (20, 40, 80):
    x = np.linspace(-L, L, 400); M = build_nonlocal_M(x, 2.0)
    print(f"    {400:6d} {L:5d} {pseudoherm_residual(M):13.2e}")

# --- B2. Convergencia del error de dispersion con N (consistencia al continuo) ---
print("\n[B2] Convergencia del error de dispersion |K·u/u - eps_lam(k)| con N (homog, lam=2, k=1.0)")
print(f"    {'N':>6} {'err rel':>10}")
for N in (200, 400, 800):
    x = np.linspace(-40, 40, N); K = weyl_kinetic(x, 2.0, mode='homog')
    u = np.exp(1j*1.0*x); Ku = K@u
    sl = slice(N//2-N//12, N//2+N//12)
    ratio = np.real(np.mean(Ku[sl]/u[sl])); es = eps_symbol(1.0, 2.0)
    print(f"    {N:6d} {abs(ratio-es)/es:10.2e}")
print("\n(nota: pseudo-herm debe ser ~maquina por construccion simetrica de K;")
print(" el contenido NO trivial del gate es la CONSISTENCIA [B/B2]: que K reproduzca")
print(" eps_lam(k) y converja al continuo. Un residuo pequeno en una malla NO basta.)")
