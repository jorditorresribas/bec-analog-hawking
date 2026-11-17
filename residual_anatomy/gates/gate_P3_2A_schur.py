"""
GATE P3.2-A (Paso 3) — EQUIVALENCIA ALGEBRAICA del sistema aumentado con K_P.
ALCANCE (ChatGPT): comprueba EXCLUSIVAMENTE que el sistema aumentado disperso, al eliminar los
auxiliares (Schur), reproduce K_P. NO es validacion de scattering ni de pseudo-hermiticidad (eso es
P3.2-B / P3.4). NO toca v1.0.

K_P a = 1/2 D a + c1 D p + c2 D^2 q3, con auxiliares resueltas por el bloque aumentado:
  (2G+D)p = Da ; (2G+D)q1 = D^2 a ; (2G+D)q2 = q1 ; (2G+D)q3 = q2.
Matriz aumentada auxiliar M_aux (4N x 4N), dispersa de banda ( (2G+D) tridiag; acoplos -I ).

Controles (ChatGPT):
 1. Residuos de las 4 restricciones auxiliares.
 2. Equivalencia de Schur: aug (solve del bloque) vs secuencial A (definicion de K_P).
 3. Referencia INDEPENDIENTE: K_P denso (A explicita, otra ruta de algebra) en N pequeno.
 4. Convergencia: variar N, L, BC.  Separar error algebraico / caja / discretizacion.
 5. Condicionamiento: cond(2G+D) y residuos de los solves de la cadena A^3.
CAVEAT (ChatGPT): comparar dos rutas con los MISMOS solves detecta errores de ensamblaje, no demuestra
que la BC discreta represente el resolvente fisico (eso = equivalencia de dominios, pendiente).
"""
import numpy as np
import scipy.sparse as sp
import scipy.sparse.linalg as spla

def build_ops(Gx, N, dx, bc='dirichlet'):
    """D=-d^2 (tridiag), D2=D@D (pentadiag), (2G+D) con BC. kappa en bordes para Robin."""
    main=2.0/dx**2*np.ones(N); off=-1.0/dx**2*np.ones(N-1)
    D=sp.diags([off,main,off],[-1,0,1]).tolil()
    twoGpD=(sp.diags(2*Gx)+D).tolil()
    if bc=='dirichlet':
        for M in (twoGpD,):
            M[0,:]=0; M[0,0]=1.0; M[-1,:]=0; M[-1,-1]=1.0
    elif bc=='robin':
        kR=np.sqrt(2*Gx[-1]); kL=np.sqrt(2*Gx[0])
        twoGpD[0,0]=2*Gx[0]+2.0/dx**2+2*kL/dx; twoGpD[0,1]=-2.0/dx**2
        twoGpD[-1,-1]=2*Gx[-1]+2.0/dx**2+2*kR/dx; twoGpD[-1,-2]=-2.0/dx**2
    D=D.tocsc(); twoGpD=twoGpD.tocsc(); D2=(D@D).tocsc()
    return D, D2, twoGpD

def rhs_bc(v, bc):
    if bc=='dirichlet':
        v=v.copy(); v[0]=0; v[-1]=0
    return v

def kp_sequential(a, D, D2, lu, bc, c1, c2):
    """Definicion de K_P via solves secuenciales de A=(2G+D)^-1."""
    def A(f): return lu.solve(rhs_bc(f,bc))
    p=A(D@a); q1=A(D2@a); q2=A(q1); q3=A(q2)
    return 0.5*(D@a)+c1*(D@p)+c2*(D2@q3), (p,q1,q2,q3)

def kp_augmented(a, D, D2, twoGpD, bc, c1, c2):
    """K_P via la MATRIZ AUMENTADA dispersa M_aux (4N) y su solve (Schur del bloque auxiliar)."""
    N=len(a); I=sp.identity(N,format='csc'); Z=sp.csc_matrix((N,N))
    # incognitas (p,q1,q2,q3); bloque lower-triangular por cadena
    M_aux=sp.bmat([[twoGpD, Z,      Z,      Z     ],
                   [Z,      twoGpD, Z,      Z     ],
                   [Z,     -I,      twoGpD, Z     ],
                   [Z,      Z,     -I,      twoGpD]],format='csc')
    rhs=np.concatenate([rhs_bc(D@a,bc), rhs_bc(D2@a,bc), rhs_bc(np.zeros(N),bc), rhs_bc(np.zeros(N),bc)])
    sol=spla.spsolve(M_aux, rhs)
    p,q1,q2,q3=sol[:N],sol[N:2*N],sol[2*N:3*N],sol[3*N:4*N]
    kp=0.5*(D@a)+c1*(D@p)+c2*(D2@q3)
    resid=M_aux@sol-rhs                     # control 1: residuos POR BLOQUE (p,q1,q2,q3)
    res4=[np.abs(resid[i*N:(i+1)*N]).max() for i in range(4)]
    return kp, (p,q1,q2,q3), res4, M_aux

def kp_dense(a, D, D2, twoGpD, bc, c1, c2):
    """Referencia INDEPENDIENTE: A explicita (inversa densa) y K_P = 1/2 D + c1 DAD + c2 D^2 A^3 D^2.
    A_eff = Adense @ Pbc aplica la MISMA BC que las rutas dispersas (Pbc anula el borde en Dirichlet;
    identidad en Robin) -> comparacion consistente."""
    N=len(a); Adense=np.linalg.inv(twoGpD.toarray())
    Pbc=np.eye(N)
    if bc=='dirichlet': Pbc[0,0]=0.0; Pbc[-1,-1]=0.0
    Aeff=Adense@Pbc
    Dd=D.toarray(); D2d=D2.toarray()
    KP=0.5*Dd + c1*(Dd@Aeff@Dd) + c2*(D2d@(Aeff@Aeff@Aeff)@D2d)
    return KP@a

def run():
    lam=2.0; c1=(lam-1)/4; c2=-(lam-1)**2/16
    print("GATE P3.2-A — equivalencia algebraica aumentado <-> K_P\n"+"="*64)
    # raiz del proyecto resuelta de forma PORTABLE (no ruta absoluta)
    import sys, os, contextlib
    from pathlib import Path
    root=Path(__file__).resolve().parents[2]       # gates -> residual_anatomy -> raíz del repo
    with contextlib.redirect_stdout(open(os.devnull,'w')):
        sys.path.insert(0,str(root / 'solver')); from correlations import bg, g
    # --- controles 1,2,3 en homogeneo y tanh ---
    for prof in ['homog','tanh']:
        for bc in ['dirichlet','robin']:
            N=300; L=20.; x=np.linspace(-L,L,N); dx=x[1]-x[0]
            Gx=np.ones(N) if prof=='homog' else g*bg(x,w=2.0,mode='horizon')[0]
            D,D2,twoGpD=build_ops(Gx,N,dx,bc)
            lu=spla.splu(twoGpD)
            a=np.exp(-(x**2)/8)+0.1*np.sin(1.1*x)      # campo de prueba suave
            kp_seq,_=kp_sequential(a,D,D2,lu,bc,c1,c2)
            kp_aug,_,res4,_=kp_augmented(a,D,D2,twoGpD,bc,c1,c2)
            kp_den=kp_dense(a,D,D2,twoGpD,bc,c1,c2)
            e_aug_seq=np.abs(kp_aug-kp_seq).max()/np.abs(kp_seq).max()
            e_aug_den=np.abs(kp_aug-kp_den).max()/np.abs(kp_den).max()
            print(f"[{prof:5s}/{bc:9s}] res(p,q1,q2,q3)=[{res4[0]:.0e},{res4[1]:.0e},{res4[2]:.0e},"
                  f"{res4[3]:.0e}]  aug~seq={e_aug_seq:.1e}  aug~denso={e_aug_den:.1e}")
    # --- control 4: convergencia en N (aug vs denso), homog Dirichlet ---
    print("\nControl 4 — convergencia en N (homog, Dirichlet), error aug~denso:")
    for N in [100,200,400,800]:
        L=20.; x=np.linspace(-L,L,N); dx=x[1]-x[0]; Gx=np.ones(N)
        D,D2,twoGpD=build_ops(Gx,N,dx,'dirichlet')
        a=np.exp(-(x**2)/8)
        kp_aug,_,_,_=kp_augmented(a,D,D2,twoGpD,'dirichlet',c1,c2)
        kp_den=kp_dense(a,D,D2,twoGpD,'dirichlet',c1,c2)
        print(f"  N={N:4d}  aug~denso(rel)={np.abs(kp_aug-kp_den).max()/np.abs(kp_den).max():.2e}")
    # --- control 4b: barrido en L (N/L ~ const => dx ~ const), fuente localizada, Dir+Robin ---
    print("\nControl 4b — barrido en L (dx~const), fuente localizada:")
    print(f"  {'bc':>9} {'L':>4} {'N':>5} {'aug~seq':>10} {'aug~denso':>10} {'res_q3':>9}")
    for bc in ['dirichlet','robin']:
        for L in [10.,20.,40.]:
            N=int(L*25); x=np.linspace(-L,L,N); dx=x[1]-x[0]; Gx=np.ones(N)
            D,D2,twoGpD=build_ops(Gx,N,dx,bc); lu=spla.splu(twoGpD)
            a=np.exp(-(x**2)/8)                      # localizada (centro), lejos de bordes al crecer L
            kp_seq,_=kp_sequential(a,D,D2,lu,bc,c1,c2)
            kp_aug,_,res4,_=kp_augmented(a,D,D2,twoGpD,bc,c1,c2)
            kp_den=kp_dense(a,D,D2,twoGpD,bc,c1,c2)
            es=np.abs(kp_aug-kp_seq).max()/np.abs(kp_seq).max()
            ed=np.abs(kp_aug-kp_den).max()/np.abs(kp_den).max()
            print(f"  {bc:>9} {L:4.0f} {N:5d} {es:10.1e} {ed:10.1e} {res4[3]:9.0e}")
    # --- control 5: condicionamiento de (2G+D) y residuos de la cadena A^3 ---
    print("\nControl 5 — cond(2G+D) y residuos de la cadena A^3 (homog, Dirichlet):")
    for N in [200,400,800]:
        L=20.; x=np.linspace(-L,L,N); dx=x[1]-x[0]; Gx=np.ones(N)
        D,D2,twoGpD=build_ops(Gx,N,dx,'dirichlet'); lu=spla.splu(twoGpD)
        a=np.exp(-(x**2)/8)
        def A(f): return lu.solve(rhs_bc(f,'dirichlet'))
        s=rhs_bc(D2@a,'dirichlet'); q1=A(s); q2=A(q1); q3=A(q2)
        r1=np.abs(twoGpD@q1-s).max(); r2=np.abs(twoGpD@q2-q1).max(); r3=np.abs(twoGpD@q3-q2).max()
        cond=np.linalg.cond(twoGpD.toarray())
        print(f"  N={N:4d}  cond(2G+D)={cond:.2e}  res_q1={r1:.1e} res_q2={r2:.1e} res_q3={r3:.1e}")
    print("\nAlcance: solo equivalencia algebraica. NO valida scattering/pseudo-herm (P3.2-B,P3.4).")
    print("Caveat: aug~seq usa los MISMOS solves (detecta ensamblaje); aug~denso es ruta independiente.")

if __name__=='__main__':
    run()
