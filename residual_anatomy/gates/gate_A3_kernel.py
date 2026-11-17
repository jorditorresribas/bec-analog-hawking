"""
GATE A3 (Paso 3, primer gate) — la realizacion dispersa de A^3 reproduce el inverso
de Green de recta completa, y su error depende de BC y de L como debe.

HIPOTESIS: A=(2G+D)^-1 band-solve con BC autoadjunta reproduce, en el interior, el
resolvente de recta completa (simbolo 1/(k^2+kappa^2), kappa=sqrt(2G)). A^3 = simbolo
1/(k^2+kappa^2)^3. NO basta coincidir puntualmente: hay que ver dependencia de BC y L.

REFERENCIAS EXACTAS (homogeneo G const):
  - onda plana:  A^3 (e^{ikx}) = e^{ikx}/(k^2+kappa^2)^3   (analitico, sin malla)
  - gaussiana:   A^3 via FFT de simbolo en caja periodica GRANDE (semilongitud 200), que APROXIMA la
                 recta completa (periodicidad despreciable para kappa=sqrt2 en la region examinada)
TOLERANCIA: en el interior (|x|<L/2), error -> piso de discretizacion FD O(dx^2) al crecer L. El
truncamiento de caja para A^3 es O(L^2 e^{-kappa L}) (las colas de A^3 son polinomio*e^{-kappa|x|}).
Dirichlet y Robin convergen; Robin (de 1er orden) MEJORA el borde frente a Dirichlet pero NO es
exactamente transparente para A^3 (la condicion (d+kappa)q3=0 ignora los factores polinomicos de las
colas; discrepancia relativa O(1/L), magnitud absoluta exp. suprimida).
NOTA (para P3.2-B): la matriz Robin NO es simetrica en el euclideo ordinario, pero SI respecto al peso
trapezoidal W=diag(1/2,1,...,1,1/2): W M = M^dagger W.
NO toca v1.0.
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def build_2GplusD(Gc, N, dx, bc):
    """(2G + D) con D=-d^2, FD centrada; bc in {'dirichlet','robin'}. kappa=sqrt(2G)."""
    kap = np.sqrt(2*Gc)
    main = (2*Gc + 2.0/dx**2)*np.ones(N)
    off  = (-1.0/dx**2)*np.ones(N-1)
    M = sp.diags([off, main, off], [-1,0,1]).tolil()
    if bc == 'dirichlet':
        # p=0 en bordes: filas identidad (resolvemos con rhs=0 en bordes)
        M[0,:]=0; M[0,0]=1.0
        M[-1,:]=0; M[-1,-1]=1.0
    elif bc == 'robin':
        # Robin de 1er orden (aprox. para colas de A^3, NO transparente exacta).
        # ghost: (d+kap)p=0 dcha -> p_{N}=p_{N-2}-2dx kap p_{N-1}; (d-kap)p=0 izda
        # nodo 0 (izda, (d-kap)p=0 -> p'=(p1-p_{-1})/2dx=kap p0 -> p_{-1}=p1-2dx kap p0)
        # stencil: (2G+2/dx^2)p0 -1/dx^2 p_{-1} -1/dx^2 p1 = f0
        #  -> (2G+2/dx^2)p0 -1/dx^2(p1-2dx kap p0) -1/dx^2 p1
        M[0,0]  = 2*Gc + 2.0/dx**2 + (2*kap/dx)
        M[0,1]  = -2.0/dx**2
        M[-1,-1]= 2*Gc + 2.0/dx**2 + (2*kap/dx)
        M[-1,-2]= -2.0/dx**2
    return M.tocsc(), kap

def make_solver(M, bc):
    """Factoriza M UNA vez (LU) y devuelve un solve que reutiliza la factorizacion
    (eficiencia: las 3 aplicaciones de A comparten el LU). Proyeccion Dirichlet + complejos."""
    lu = spla.splu(M)
    def solve(f):
        rhs = f.copy()
        if bc == 'dirichlet':
            rhs = rhs.astype(f.dtype); rhs[0]=0; rhs[-1]=0
        if np.iscomplexobj(rhs):
            return lu.solve(rhs.real) + 1j*lu.solve(rhs.imag)
        return lu.solve(rhs)
    return solve

def A3_bandsolve(f, M, bc):
    solve = make_solver(M, bc)           # LU factorizado UNA vez, reutilizado 3 veces
    q1 = solve(f); q2 = solve(q1); q3 = solve(q2)
    return q3

def grid(L, dx):
    """Malla centrada SIMETRICA: N impar, ambos extremos +-L incluidos, espaciado EXACTO dx, 0 incluido."""
    M = int(round(L/dx))
    return (np.arange(2*M+1) - M)*dx     # [-L, ..., 0, ..., +L]

def gauss_reference(x, dx, kap, sigma=1.0, Lr=200.0):
    """A^3 (aprox. recta completa) para gaussiana, via FFT en caja periodica GRANDE
    ALINEADA con x (mismo dx, misma malla simetrica) -> referencia por SLICE, sin interp."""
    xr = grid(Lr, dx); Nr = xr.size
    fr = np.exp(-(xr**2)/(2*sigma**2))
    kk = 2*np.pi*np.fft.fftfreq(Nr, d=dx)
    ref_full = np.fft.ifft(np.fft.fft(fr)/(kk**2+kap**2)**3).real
    i0 = Nr//2 + np.round(x/dx).astype(int)      # indices de x dentro de xr (alineados)
    return ref_full[i0]

def run():
    Gc = 1.0; kap = np.sqrt(2*Gc)
    dx = 0.02
    print(f"GATE A3  (G={Gc}, kappa={kap:.4f}, dx={dx})")
    print(f"{'src':>8} {'bc':>10} {'L':>5} {'N':>6} {'err_interior':>14} {'err_borde':>12}")
    # fuente gaussiana (localizada) y onda plana k=0.8
    for src in ['gauss','plane']:
        for L in [10., 20., 40.]:
            x = grid(L, dx); N = x.size
            if src=='gauss':
                f = np.exp(-(x**2)/2.0)
                ref = gauss_reference(x, dx, kap)       # recta completa, sin interp
            else:
                kw=0.8; f=np.exp(1j*kw*x)
                ref=np.exp(1j*kw*x)/(kw**2+kap**2)**3   # analitico exacto
            for bc in ['dirichlet','robin']:
                M,_=build_2GplusD(Gc,N,dx,bc)
                q3=A3_bandsolve(f.astype(complex) if src=='plane' else f, M, bc)
                err=np.abs(q3-ref)
                msk=np.abs(x)<L*0.5
                ei=err[msk].max(); eb=err[:5].max()+err[-5:].max()
                print(f"{src:>8} {bc:>10} {L:5.0f} {N:6d} {ei:14.3e} {eb:12.3e}")
        print()
    print("Lectura: err_interior -> piso FD (O(dx^2)) al crecer L; borde de Dirichlet mayor")
    print("que Robin a L fijo; ambos convergen. Confirma que la BC implementa el inverso correcto.")

def convergence_dx():
    """Confirma que el error INTERIOR escala O(dx^2): L=40 fijo (borde despreciable),
    gaussiana, Dirichlet, barrido en dx. pendiente log-log ~ 2."""
    Gc=1.0; kap=np.sqrt(2*Gc); L=40.
    print("\nConvergencia en dx (L=40 fijo, gauss, Dirichlet):")
    print(f"{'dx':>8} {'N':>7} {'err_interior':>14} {'ratio':>8}")
    errs=[]; dxs=[0.08,0.04,0.02,0.01]
    for dx in dxs:
        x=grid(L,dx); N=x.size
        f=np.exp(-(x**2)/2.0)
        ref=gauss_reference(x, dx, kap)
        M,_=build_2GplusD(Gc,N,dx,'dirichlet')
        q3=A3_bandsolve(f,M,'dirichlet')
        ei=np.abs(q3-ref)[np.abs(x)<L*0.5].max(); errs.append(ei)
        r = errs[-2]/ei if len(errs)>1 else float('nan')
        print(f"{dx:8.3f} {N:7d} {ei:14.3e} {r:8.2f}")
    slope=np.polyfit(np.log(dxs),np.log(errs),1)[0]
    print(f"pendiente log-log = {slope:.2f}  (esperado ~2 => O(dx^2))")

if __name__=='__main__':
    run()
    convergence_dx()
