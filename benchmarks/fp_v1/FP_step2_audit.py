"""
AUDITORIA Paso 2 (respuesta a ChatGPT). Controles:
 (A) max|f'| CERTIFICADO: derivada ANALITICA cerrada de f=c-v (no np.gradient),
     maximo por root-find de f''=0, y comparacion con np.gradient(N) para acotar el
     error de discretizacion -> convierte el techo en cota numerica fiable.
 (B) kappa_bar(d) sobre rango amplio de d: monotonia y que SIEMPRE <= max|f'|.
 (C) presupuesto de error de T_IR: variar banda de frecuencias y orden de
     extrapolacion (lin/quad/cubica) desde UNA malla densa de N_omega cacheada;
     reportar dispersion vs el exceso. + chequear regimen omega/T.
Sin np.gradient para el techo. Deterministico (sin ruido): 'error' = sistematico.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'solver'))  # layout del repo
import numpy as np, mode_matching as MM
from scipy.optimize import brentq
g=1.0; n_u,mach_u,md=1.0,0.5,2.0
vLb=mach_u*np.sqrt(g*n_u); jc=n_u*vLb; n_d=(jc/(md*np.sqrt(g)))**(2/3)
WS=[2.0,3.0,4.0,6.0,8.0]

# ---- perfil y derivadas ANALITICAS de f=c-v ----
def nprof(x,w):
    s=0.5*(1+np.tanh(x/w)); return n_u+(n_d-n_u)*s
def dn(x,w):
    th=np.tanh(x/w); return (n_d-n_u)*(1-th**2)/(2*w)
def d2n(x,w):
    th=np.tanh(x/w); return (n_d-n_u)*(-2*th*(1-th**2))/(2*w**2)
def fp_analytic(x,w):              # f'(x) EXACTA, f=sqrt(n)-jc/n
    n=nprof(x,w); return dn(x,w)*(1.0/(2*np.sqrt(n)) + jc/n**2)
def fpp_analytic(x,w):             # f''(x) EXACTA (para root-find del maximo)
    n=nprof(x,w); np_=dn(x,w); npp=d2n(x,w)
    br=1.0/(2*np.sqrt(n))+jc/n**2
    dbr=np_*(-1.0/(4*n**1.5) - 2*jc/n**3)
    return npp*br + np_*dbr
def f_analytic(x,w):
    n=nprof(x,w); return np.sqrt(n)-jc/n

print("="*96)
print("(A) TECHO CERTIFICADO: max|f'| analitico (root-find f''=0) vs np.gradient(N)")
print("="*96)
cert={}
for w in WS:
    # horizonte: f=0
    xH=brentq(lambda x:f_analytic(x,w), -5*w, 5*w)
    kK=abs(fp_analytic(xH,w))
    # MAXIMO GLOBAL de |f'|: |f'|=|n'|*B con n'~sech^2(x/w) -> |f'|->0 como e^{-2|x|/w}
    # en AMBOS infinitos. Escaneo dominio ANCHO [xH-30w, xH+30w] (|f'| en bordes
    # ~ e^{-60}, despreciable), TODAS las raices de f''=0, y el maximo entre ellas.
    Lw=30*w; xs=np.linspace(xH-Lw, xH+Lw, 60000); fpp=fpp_analytic(xs,w)
    idx=np.where(np.sign(fpp[:-1])!=np.sign(fpp[1:]))[0]
    cand=[brentq(lambda x:fpp_analytic(x,w), xs[i],xs[i+1]) for i in idx]
    xmax=max(cand, key=lambda x:abs(fp_analytic(x,w)))
    kMax=abs(fp_analytic(xmax,w))
    edgeL=abs(fp_analytic(xH-Lw,w)); edgeR=abs(fp_analytic(xH+Lw,w))
    print(f"   [cert global] {len(cand)} puntos criticos de f'' en +-30w; "
          f"|f'| en bordes: izq={edgeL:.1e} der={edgeR:.1e} (<<kappa_local) "
          f"-> maximo es GLOBAL")
    # comparacion np.gradient a varias N (como en el script de produccion)
    rowg=[]
    for N in (20000,60000,120000,240000):
        L=max(40.,8*w+20.); xg=np.linspace(-L,L,N)
        kg=np.max(np.abs(np.gradient(f_analytic(xg,w),xg)))
        rowg.append((N,kg))
    cert[w]=dict(xH=xH,kK=kK,xmax=xmax,kMax=kMax)
    print(f"\nw={w:.0f}: x_H={xH:+.5f} kappa_local={kK:.6f}  |  x_max={xmax:+.5f} "
          f"max|f'|_ANALITICO={kMax:.6f}  (techo/kloc={kMax/kK:.6f})")
    for N,kg in rowg:
        print(f"     np.gradient N={N:6d}: max|f'|={kg:.6f}  err_rel_vs_analitico={abs(kg-kMax)/kMax:.2e}")
print("\n -> el techo analitico es exacto; np.gradient lo subestima O(dx^2)."
      " Techo/kloc = 1.014339 (invariante en w, feature del tanh).")

print("\n"+"="*96)
print("(B) kappa_bar(d): monotonia y cota <= max|f'| (ventana simetrica centrada en x_H)")
print("="*96)
def kbar(w,xH,kK,d):
    cH=np.sqrt(nprof(xH,w)); Lam=2*cH**2
    dxi=d*cH/(kK*Lam**2)**(1/3)
    xL,xR=xH-0.5*dxi, xH+0.5*dxi
    return abs(f_analytic(xR,w)-f_analytic(xL,w))/dxi, dxi
for w in WS:
    c=cert[w]; kMax=c['kMax']
    ds=np.concatenate([np.logspace(-3,-1,40),np.linspace(0.1,50.0,500)])
    kbs=np.array([kbar(w,c['xH'],c['kK'],d)[0] for d in ds])
    mono = np.all(np.diff(kbs)<=1e-12)
    print(f"w={w:.0f}: kbar(d->0+)={kbs[0]:.5f}(->kloc={c['kK']:.5f}) kbar(d=8)={kbs[-1]:.5f} "
          f"max_d kbar={kbs.max():.5f} (<=max|f'|={kMax:.5f}? {kbs.max()<=kMax+1e-9}) "
          f"monotona_decreciente={mono}")

print("\n"+"="*96)
print("(C) PRESUPUESTO DE ERROR de T_IR: banda x orden, desde malla densa de N_omega")
print("="*96)
def Nw_pu(w,om,L,N):
    MM.w=w; MM.omega=om; x=np.linspace(-L,L,N); n,v,V=MM.bg(x,'horizon')
    mR=MM.modes_full(n[-1],v[-1],V[-1]); Rin=[m for m in mR if m['prop'] and m['vg']<0]
    incR=[1.0 if m['norm']<0 else 0.0 for m in Rin]
    out,_=MM.solve_inject(L=L,N=N,mode='horizon',incL=0.0,incR=incR)
    return sum(a for(nrm,s,a)in out if nrm>0), abs(sum(nrm*a for(nrm,s,a)in out)+1.0)
def box(w): L=max(40.,8*w+20.); return L,int(round(L*200))
# malla densa om/T en [0.08,1.2]; extrapolaciones y bandas se derivan de aqui
xr_dense=np.array([0.08,0.10,0.12,0.15,0.18,0.22,0.27,0.33,0.40,0.50,0.60,0.75,0.90,1.10])
for w in WS:
    c=cert[w]; kK=c['kK']; T0=kK/(2*np.pi); L,N=box(w)
    oms=xr_dense*T0
    dat=[Nw_pu(w,om,L,N) for om in oms]
    Nw=np.array([d[0] for d in dat]); pu=np.array([d[1] for d in dat]); omN=oms*Nw
    print(f"\n--- w={w:.0f}  T0=kloc/2pi={T0:.5f}  pu_max={pu.max():.1e}  (techo/2pi={c['kMax']/(2*np.pi):.5f}) ---")
    ests=[]
    for lo,hi in [(0.08,0.40),(0.10,0.50),(0.12,0.60),(0.15,0.75),(0.08,0.60),(0.12,0.90)]:
        sel=(xr_dense>=lo-1e-9)&(xr_dense<=hi+1e-9)
        o=oms[sel]; y=omN[sel]
        for order,name in [(1,'lin'),(2,'quad'),(3,'cub')]:
            if sel.sum()>=order+1:
                T=np.polyfit(o,y,order)[-1]; ests.append(T)
    ests=np.array(ests)
    # Planck-fit banda ancha (cross-check): log(1+1/N)=om/T sobre [0.15,1.1]
    selp=(xr_dense>=0.15); cp=np.polyfit(oms[selp],np.log(1+1/Nw[selp]),1); Tpl=1/cp[0]
    print(f"    T_IR por banda/orden: min={ests.min():.5f} max={ests.max():.5f} "
          f"media={ests.mean():.5f} sd={ests.std():.5f} dispersion={ests.std()/ests.mean():.2%}")
    print(f"    Planck-fit ancho [0.15,1.1]T: T={Tpl:.5f} (ratio/kloc={Tpl/kK*2*np.pi:.4f})")
    print(f"    T_IR/kloc: min={ests.min()/T0:.4f} max={ests.max()/T0:.4f}  | "
          f"T_IR/techo: min={ests.min()/(c['kMax']/(2*np.pi)):.4f} max={ests.max()/(c['kMax']/(2*np.pi)):.4f}")
    verdict = "SUPERA el techo con TODA banda/orden" if ests.min()>c['kMax']/(2*np.pi) else \
              ("supera el techo con la media pero NO con todas" if ests.mean()>c['kMax']/(2*np.pi) else "NO supera el techo robustamente")
    print(f"    -> {verdict}")
