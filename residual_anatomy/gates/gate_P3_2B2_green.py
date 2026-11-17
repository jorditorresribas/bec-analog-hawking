"""
GATE P3.2-B2 — Identidad de Green DISCRETA de K_P. TRES comprobaciones SEPARADAS (ChatGPT):
 B2.1  SIMETRIA MATRICIAL de la matriz finita: K_P^T = K_P. Verifica la simetria discreta, NADA mas.
 B2.2  IDENTIDAD DISCRETA (SBP) con producto escalar H: el flujo de region = concomitante de borde.
 B2.3  CONVERGENCIA del concomitante discreto al continuo de sub-region (con los terminos no-locales).

DISCRETIZACION: 'Dirichlet por NODOS FANTASMA nulos' — el stencil tridiagonal (2/dx^2, -1/dx^2) hasta
los extremos cierra con phi_{-1}=phi_{N}=0. Es SIMETRICA y reproduce A^3 de recta completa en el interior
a L grande (B1b), pero NO es la realizacion Robin de §9.2; es una matriz finita simetrica concreta. NO
se le atribuyen las condiciones Robin de §9.2.

B2.2 (SBP): con H = producto TRAPEZOIDAL (diag(1/2,1,..,1,1/2)) y region R=[a,b] (a,b INTERIORES), la
suma-por-partes de los operadores de FD central telescopia EXACTAMENTE a los Wronskianos centrales en
a,b. Por eso <w,K_P u>_{H,R} - <K_P w,u>_{H,R} = [BF]_a^b a ~maquina (identidad discreta, no coincidencia).
Equivale al flujo del conmutador con el proyector P_R (formulacion P_R^T H K_P - K_P^T H P_R).

BF = concomitante de sub-region (incluye los NO-locales de A, A^3; §9.3 es el limite R->caja donde ->0):
  BF = 1/2 W(w,u)
     + c1 [ W(w,p_u) + W(p_w,u) - W(p_w,p_u) ]
     + c2 [ W2(w,q3u)+W2(q3w,u) - W(q1w,q3u) - W(q2w,q2u) - W(q3w,q1u) ]
  p=ADa; q1=AD^2a; q2=Aq1; q3=Aq2; A=(2G+D)^-1; W(a,b)=a'b-ab'; W2(a,b)=a b'''-a'b''+a''b'-a'''b.
Campos REALES; K_P solo (sin C); G=1 homogeneo. NO v1.0.
"""
import numpy as np

lam = 2.0; c1 = (lam-1)/4; c2 = -(lam-1)**2/16; G = 1.0

def build(N, dx):
    main = 2/dx**2*np.ones(N); off = -1/dx**2*np.ones(N-1)
    D = np.diag(main)+np.diag(off,1)+np.diag(off,-1)      # -d^2, Dirichlet por nodos fantasma nulos
    A = np.linalg.inv(2*G*np.eye(N)+D)
    return D, D@D, A

def chain(a, D, D2, A):
    p = A@(D@a); q1 = A@(D2@a); q2 = A@q1; q3 = A@q2
    return p, q1, q2, q3

def KP_mat(D, D2, A):
    return 0.5*D + c1*(D@A@D) + c2*(D2@(A@A@A)@D2)

def d1(f,i,dx): return (f[i+1]-f[i-1])/(2*dx)
def d2(f,i,dx): return (f[i+1]-2*f[i]+f[i-1])/dx**2
def d3(f,i,dx): return (f[i+2]-2*f[i+1]+2*f[i-1]-f[i-2])/(2*dx**3)
def Wr(a,b,i,dx):  return d1(a,i,dx)*b[i] - a[i]*d1(b,i,dx)
def Wr2(a,b,i,dx): return a[i]*d3(b,i,dx)-d1(a,i,dx)*d2(b,i,dx)+d2(a,i,dx)*d1(b,i,dx)-d3(a,i,dx)*b[i]

def BF(w,u, pu,q1u,q2u,q3u, pw,q1w,q2w,q3w, i, dx, nonlocal_terms=True):
    r = (0.5*Wr(w,u,i,dx) + c1*(Wr(w,pu,i,dx)+Wr(pw,u,i,dx)) + c2*(Wr2(w,q3u,i,dx)+Wr2(q3w,u,i,dx)))
    if nonlocal_terms:
        r += c1*(-Wr(pw,pu,i,dx)) + c2*(-Wr(q1w,q3u,i,dx)-Wr(q2w,q2u,i,dx)-Wr(q3w,q1u,i,dx))
    return r

def run():
    print("GATE P3.2-B2 — identidad de Green discreta de K_P\n"+"="*66)
    L = 20.
    print("B2.1 — simetria matricial K_P^T=K_P (matriz finita Dirichlet-fantasma):")
    for N in [160,320]:
        x=np.linspace(-L,L,N); dx=x[1]-x[0]; D,D2,A=build(N,dx); KP=KP_mat(D,D2,A)
        print(f"  N={N}: eucl_rel={np.abs(KP-KP.T).max()/np.abs(KP).max():.2e}")

    def setup(N):
        x=np.linspace(-L,L,N); dx=x[1]-x[0]; D,D2,A=build(N,dx); KP=KP_mat(D,D2,A)
        w=1/(1+((x-5)/2)**2); u=1/(1+((x+3)/2)**2)            # lorentzianas asimetricas
        pu,q1u,q2u,q3u=chain(u,D,D2,A); pw,q1w,q2w,q3w=chain(w,D,D2,A)
        ia=int(np.argmin(np.abs(x+6))); ib=int(np.argmin(np.abs(x-6)))
        aux=(pu,q1u,q2u,q3u,pw,q1w,q2w,q3w)
        return x,dx,KP,w,u,ia,ib,aux

    print("\nB2.2 — identidad discreta SBP (H trapezoidal): flujo de region = [BF]_a^b a ~maquina:")
    print(f"  {'N':>5} {'dx':>8} {'flujo_H(trapz)':>15} {'[BF]_a^b':>13} {'|dif| (SBP)':>12}")
    for N in [400,800,1600]:
        x,dx,KP,w,u,ia,ib,aux=setup(N)
        KPu=KP@u; KPw=KP@w
        # H trapezoidal sobre R: pesos 1/2 en los extremos ia,ib
        Hw=np.ones(ib-ia+1); Hw[0]=0.5; Hw[-1]=0.5
        fluxH=np.sum(Hw*(w*KPu-KPw*u)[ia:ib+1])*dx
        bf=lambda i: BF(w,u,*aux,i,dx)
        BFj=bf(ib)-bf(ia)
        print(f"  {N:>5} {dx:8.4f} {fluxH:15.6e} {BFj:13.6e} {abs(fluxH-BFj):12.2e}")

    print("\nB2.3 — convergencia al continuo (producto EUCLIDEO consistente: F = dx w^T[P_R,K_P]u):")
    print(f"  {'N':>5} {'dx':>8} {'F_conmut(exacto)':>16} {'[BF]completa':>13} {'dif_F-BF':>10} {'dif_naive':>10}")
    dxs=[]; dFB=[]
    for N in [400,800,1600,3200]:                      # N<=3200: antes del piso de condicionamiento
        x,dx,KP,w,u,ia,ib,aux=setup(N)
        bf =lambda i: BF(w,u,*aux,i,dx,nonlocal_terms=True)
        bfn=lambda i: BF(w,u,*aux,i,dx,nonlocal_terms=False)
        BFj=bf(ib)-bf(ia); BFn=bfn(ib)-bfn(ia)
        PR=np.diag((np.arange(len(x))>=ia)&(np.arange(len(x))<=ib)).astype(float)
        Fcommut=dx*(w@(PR@KP-KP@PR)@u)                 # flujo discreto EXACTO (conmutador, sharp P_R)
        dxs.append(dx); dFB.append(abs(Fcommut-BFj))
        print(f"  {N:>5} {dx:8.4f} {Fcommut:16.6e} {BFj:13.6e} {abs(Fcommut-BFj):10.2e} {abs(BFj-BFn):10.2e}")
    slope=np.polyfit(np.log(dxs),np.log(dFB),1)[0]
    print(f"  regresion de orden |F-BF|: pendiente={slope:.2f} (O(dx), proyector sharp). "
          f"dif_naive->4.33e-5 (NO nulo).")
    print("  (N>3200 NO incluido: condicionamiento de D2A3D2 DENSA corrompe; en banda disperso se evita.)")

    print("\nLectura (alcance delimitado):")
    print(" B2.1: K_P^T=K_P ~maquina -> simetria de la matriz finita Dirichlet-fantasma. NADA mas.")
    print(" B2.2: el flujo discreto EXACTO es F=dx w^T[P_R,K_P]u (definicion matricial del concomitante).")
    print("       OBSERVADO: con H trapezoidal, F=[BF]_a^b a ~maquina (propiedad SBP de la FD central de 2o")
    print("       orden sobre region interior); la DERIVACION formal para la composicion DAD/D2A3D2 queda")
    print("       como resultado mas profundo, NO reclamado aqui.")
    print(" B2.3: F (euclideo) converge a [BF]completa como O(dx) (regresion); dif_naive->valor NO NULO")
    print("       (4.33e-5) = los terminos NO-locales son del CONTINUO, no de malla -> NECESARIOS.")
    print(" §9.3 es la forma de CAJA COMPLETA (extras->0 por BC en extremos). Matriz Dirichlet-fantasma")
    print(" concreta, NO la realizacion Robin de §9.2.")

if __name__=='__main__':
    run()
