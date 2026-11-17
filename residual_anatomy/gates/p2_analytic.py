"""
PASO 2 — verificaciones analiticas consolidadas (reproducibles).
Documentacion: PADE_SPEC.md §9. NO toca v1.0.
Ejecuta P2.1..P2.7 y reporta PASS/valores. Requiere numpy, scipy, sympy.
"""
import numpy as np, scipy.sparse as sp, scipy.sparse.linalg as spla, sys, os, contextlib
G=1.0; lam=2.0; c1=(lam-1)/4.0; c2=-(lam-1)**2/16.0

def P2_1_bc_count():
    k=1.3; K2,K3=k**2,k**3
    M=np.array([[0,1,0,0,0,0],[K2,0,0,0,0,0],[0,0,0,1,0,0],
                [-1,0,K2,0,0,0],[0,0,0,0,0,1],[0,0,-1,0,K2,0]],float)
    ev=np.linalg.eigvals(M)
    Gm=np.linalg.matrix_power(M+k*np.eye(6),3); ns=6-np.linalg.matrix_rank(Gm)
    _,s,vt=np.linalg.svd(Gm); stable=vt[6-ns:].conj().T
    B=np.array([[k,1,0,0,0,0],[-1,0,2*K2,2*k,0,0],[0,0,-3*k,-1,4*K3,4*K2]],float)
    ok = (ns==3) and (np.linalg.matrix_rank(B)==3) and (np.abs(B@stable).max()<1e-12)
    return ok, f"dim_estable={ns}, rank(B)={np.linalg.matrix_rank(B)}, |B*estable|={np.abs(B@stable).max():.1e}"

def P2_2_positivity():
    # homogeneo G=1
    L=30.;N=1500;x=np.linspace(-L,L,N);dx=x[1]-x[0]
    off=-1/dx**2*np.ones(N-1)
    homog=sp.diags([off,(2*G+2/dx**2)*np.ones(N),off],[-1,0,1]).tocsc()
    ev_h=spla.eigsh(homog,k=1,which='SA',return_eigenvectors=False)[0]
    # tanh real
    with contextlib.redirect_stdout(open(os.devnull,'w')):
        sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[2] / 'solver')); from correlations import bg, g
    L2=40.;N2=3000;x2=np.linspace(-L2,L2,N2);dx2=x2[1]-x2[0]
    n,v,V=bg(x2,w=2.0,mode='horizon'); Gx=g*n
    off2=-1/dx2**2*np.ones(N2-1)
    tanh=sp.diags([off2,2*Gx+2/dx2**2,off2],[-1,0,1]).tocsc()
    ev_t=spla.eigsh(tanh,k=1,which='SA',return_eigenvectors=False)[0]
    ok=(ev_h>0) and (ev_t>0)
    return ok, f"min_autoval homog(G=1)={ev_h:.3f}, tanh_real={ev_t:.3f} (2G_min={2*min(Gx):.3f})"

def P2_3_symmetry():
    L=40.;N=2000;x=np.linspace(-L,L,N);dx=x[1]-x[0]
    m=2*np.ones(N);o=-1*np.ones(N-1)
    D=(np.diag(m)+np.diag(o,1)+np.diag(o,-1))/dx**2
    A=np.linalg.inv(2*G*np.eye(N)+D);A3=A@A@A;D2=D@D
    KP=0.5*D+c1*(D@A@D)+c2*(D2@A3@D2)
    rel=np.abs(KP-KP.T).max()/np.abs(KP).max()
    return rel<1e-7, f"asim_rel(K_P)={rel:.1e} (cond. de D2 densa, ~maquina en disperso)"

def P2_4_sign():
    Lam=4.;N=8000;x=np.linspace(0,Lam,N+1);h=x[1]-x[0]
    u=np.sin(1.3*x)+0.5*x; w=np.cos(0.7*x)+0.2*x**2
    def d2(f):
        g=np.zeros_like(f);g[1:-1]=-(f[2:]-2*f[1:-1]+f[:-2])/h**2;return g
    main=(2*G+2/h**2)*np.ones(N-1);off=(-1/h**2)*np.ones(N-2)
    Lm=sp.diags([off,main,off],[-1,0,1]).tocsc()
    def A(r):
        p=np.zeros(N+1);p[1:-1]=spla.spsolve(Lm,r[1:-1]);return p
    pu=A(d2(u));pw=A(d2(w));DADu=d2(pu);DADw=d2(pw)
    sl=slice(2,N-1)
    LHS=np.trapezoid((w*DADu)[sl],x[sl])-np.trapezoid((DADw*u)[sl],x[sl])
    def dv(f,i,s):
        return (3*f[i]-4*f[i-1]+f[i-2])/(2*h) if s=='R' else (-3*f[i]+4*f[i+1]-f[i+2])/(2*h)
    def W(a,b):
        return (dv(a,N,'R')*b[N]-a[N]*dv(b,N,'R'))-(dv(a,0,'L')*b[0]-a[0]*dv(b,0,'L'))
    Wp=W(w,pu)+W(pw,u); Wm=W(w,pu)-W(pw,u)
    ok=abs(LHS-Wp)<abs(LHS-Wm) and abs(LHS-Wp)<1e-2
    return ok, f"LHS={LHS:.4f} RHS+={Wp:.4f}(dif {abs(LHS-Wp):.1e}) RHS-={Wm:.4f}"

def P2_5_w2_symbolic():
    import sympy as sp2
    k,Gg,Gp,q3,q2,u,up,upp,uppp=sp2.symbols("kappa G Gp q3 q2 u up upp uppp")
    q3p=-k*q3; q3pp=2*Gg*q3-q2; q3ppp=2*Gp*q3+2*Gg*q3p-(-k*q2)
    W2=u*q3ppp-up*q3pp+upp*q3p-uppp*q3
    W2=W2.subs({2*Gg:k**2}).subs({up:-k*u,uppp:-k*upp})
    res_Gp0=sp2.simplify(W2.subs(Gp,0)); res_full=sp2.simplify(W2)
    ok=(res_Gp0==0) and (sp2.simplify(res_full-2*Gp*q3*u)==0)
    return ok, f"W2|G'=0={res_Gp0}; W2|G'!=0={res_full}"

def P2_6_gprime_convergence():
    """Residuo de borde completo R_c2 = 2 c2 G'(x_b)[ubar q3u + q3wbar u] del perfil tanh
    real. Para AISLAR G' se usa una ONDA PLANA de amplitud FIJA (u(x_b)~O(1), no decae con L):
    entonces la unica dependencia de L viene de G'(x_b)~e^{-L/w}, y R/|G'(x_b)| ~ const.
    Se incluye tambien la gaussiana como control de fuente localizada (doble supresion).
    NOTA: A aqui es tipo DIRICHLET (tridiagonal sin filas Robin); este test prueba el termino
    de G' (externo), independiente de la BC interna de A -> no es un test de la Robin de P3.1."""
    with contextlib.redirect_stdout(open(os.devnull,'w')):
        sys.path.insert(0,str(__import__('pathlib').Path(__file__).resolve().parents[2] / 'solver')); from correlations import bg, g
    w=2.0; kw=0.8
    def Gprime(xb):
        t=np.tanh(xb/w); n_d=(0.5/(2*np.sqrt(g)))**(2/3); Dn=n_d-1.0
        return g*Dn*(1-t**2)/(2*w)
    def q3_of(u,A,D2): return A(A(A(D2(u))))
    res_g=[]; res_p=[]; ratio_p=[]
    for L in [10.,20.,30.]:
        N=int(L/0.04); x=np.linspace(-L,L,N); dx=x[1]-x[0]
        n,v,V=bg(x,w=w,mode='horizon'); Gx=g*n
        main=2*Gx+2/dx**2; off=-1/dx**2*np.ones(N-1)
        lu=spla.splu(sp.diags([off,main,off],[-1,0,1]).tocsc())
        def A(f): return lu.solve(np.asarray(f.real))+1j*lu.solve(np.asarray(f.imag)) if np.iscomplexobj(f) else lu.solve(f)
        def D2(f):
            f1=np.zeros_like(f); f1[1:-1]=-(f[2:]-2*f[1:-1]+f[:-2])/dx**2
            f2=np.zeros_like(f); f2[1:-1]=-(f1[2:]-2*f1[1:-1]+f1[:-2])/dx**2; return f2
        # fuente localizada (gaussiana): doble supresion (control)
        ug=np.exp(-(x**2)/(2*1.5**2)); q3g=q3_of(ug,A,D2)
        res_g.append(abs(2*Gprime(x[-1])*ug[-1]*q3g[-1])+abs(2*Gprime(x[0])*ug[0]*q3g[0]))
        # onda plana amplitud FIJA: residual completo de 2 terminos (u=w) => 2 Re(ubar q3u)
        up=np.exp(1j*kw*x).astype(complex); q3p=q3_of(up,A,D2)
        def Rfull(xb_i):
            return abs(2*(-(lam-1)**2/16)*Gprime(x[xb_i])*2*np.real(np.conj(up[xb_i])*q3p[xb_i]))
        Rp=Rfull(-1)+Rfull(0); res_p.append(Rp)
        ratio_p.append(Rp/(abs(Gprime(x[-1]))+abs(Gprime(x[0]))))   # R / |G'| : debe ser ~const
    ok = (res_p[0]>res_p[1]>res_p[2]) and (max(ratio_p)/min(ratio_p) < 3.0)  # R sigue a G'
    return ok, (f"ONDA PLANA R_c2 L=10,20,30: {res_p[0]:.1e},{res_p[1]:.1e},{res_p[2]:.1e}; "
                f"R/|G'|={ratio_p[0]:.2e},{ratio_p[1]:.2e},{ratio_p[2]:.2e} (~const=>aisla G'); "
                f"gauss(control): {res_g[0]:.1e},{res_g[1]:.1e},{res_g[2]:.1e}")

def P2_7_flux():
    def eps_l(k):return -G+np.sqrt(G**2+G*k**2+(lam/4)*k**4)
    def eps_lp(k):return (G*k+(lam/2)*k**3)/(eps_l(k)+G)
    def eps_P(s):return s/2+(lam-1)*s**2/(4*(2*G+s))-(lam-1)**2*s**4/(16*(2*G+s)**3)
    h=1e-6;ks=np.linspace(0.02,1.6,400)
    dJ=np.array([abs((eps_P((k+h)**2)-eps_P((k-h)**2))/(2*h)-eps_lp(k))/abs(eps_lp(k)) for k in ks])
    return dJ.max()<0.01, f"max_k delta_J={dJ.max()*100:.3f}% en k={ks[dJ.argmax()]:.2f} (<1%)"

if __name__=='__main__':
    print("PASO 2 — verificaciones analiticas (PADE_SPEC §9)\n"+"="*60)
    for name,fn in [("P2.1 recuento BC",P2_1_bc_count),("P2.2 positividad",P2_2_positivity),
                    ("P2.3 simetria K_P",P2_3_symmetry),("P2.4 signo Green",P2_4_sign),
                    ("P2.5 residuo W2",P2_5_w2_symbolic),("P2.6 converg G'",P2_6_gprime_convergence),
                    ("P2.7 flujo dJ",P2_7_flux)]:
        ok,msg=fn()
        print(f"[{'PASS' if ok else 'FALLA':>5}] {name:22s} {msg}")
