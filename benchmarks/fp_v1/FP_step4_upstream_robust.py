"""
PASO 4 - robustez del observable de FP (ChatGPT §5.4): T_upstream = T(n_L) = T_0 del
canal u (n_omega=|beta|^2), estable al variar BANDA y RESOLUCION. FP define T_0 via
n_omega=1/(e^{omega/T}-1) -> T_omega=omega/ln(1+1/n_L); T_0=lim_{omega->0}.
Comparacion principal con FP: T_upstream / (kappa_bar/2pi).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'solver'))  # layout del repo
import numpy as np, mode_matching as MM
g=MM.g; n_u,mach_u,md=1.0,0.5,2.0
vLb=mach_u*np.sqrt(g*n_u); jc=n_u*vLb; n_d=(jc/(md*np.sqrt(g)))**(2/3)
def box(w): L=max(40.,8*w+20.); return L,int(round(L*240))
def geom(w):
    L,_=box(w); x=np.linspace(-L,L,60000); s=0.5*(1+np.tanh(x/w)); n=n_u+(n_d-n_u)*s
    c=np.sqrt(n); v=jc/n; f=c-v; i=np.where(np.sign(f[:-1])!=np.sign(f[1:]))[0][0]
    t=f[i]/(f[i]-f[i+1]); xH=x[i]+t*(x[i+1]-x[i]); cH=np.interp(xH,x,c)
    kK=abs(np.interp(xH,x,np.gradient(f,x)))
    Lam=2*cH**2; dxi=1.58*cH/(kK*Lam**2)**(1/3)
    kbar=abs(np.interp(xH+0.5*dxi,x,f)-np.interp(xH-0.5*dxi,x,f))/dxi
    return kK, kbar
def nL(w,om,L,N):
    MM.w=w; MM.omega=om; x=np.linspace(-L,L,N); n,v,V=MM.bg(x,'horizon')
    Rin=[m for m in MM.modes_full(n[-1],v[-1],V[-1]) if m['prop'] and m['vg']<0]
    incR=[1.0 if m['norm']<0 else 0.0 for m in Rin]
    out,_=MM.solve_inject(L=L,N=N,mode='horizon',incL=0.0,incR=incR)
    return sum(a for(q,s,a)in out if q>0 and s=='L')
print("="*96)
print("PASO 4 - T_upstream=T(n_L) estable en banda y resolucion; vs kappa_bar/2pi (FP)")
print("="*96)
for w in (2.0,4.0,8.0):
    kK,kbar=geom(w); T0=kK/(2*np.pi); Tkbar=kbar/(2*np.pi); L,N=box(w)
    xr=np.array([0.05,0.065,0.08,0.10,0.12,0.15,0.18,0.22,0.27,0.33,0.40,0.50,0.60]); oms=xr*T0
    nv=np.array([nL(w,om,L,N) for om in oms])
    Tw=oms/np.log(1+1/nv)           # T_omega del canal u (def FP)
    print(f"\n=== w={w:.0f}  kappa_local/2pi={T0:.5f}  kappa_bar/2pi={Tkbar:.5f} (kbar/kloc={kbar/kK:.4f}) ===")
    print(f"  {'banda/T0':>12} {'npts':>5} {'T0_slope/kloc':>13} {'T0_extrap/kloc':>14} {'T_up/(kbar/2pi)':>15}")
    for lo,hi in [(0.05,0.25),(0.08,0.40),(0.12,0.60),(0.05,0.60)]:
        sel=(xr>=lo-1e-9)&(xr<=hi+1e-9); o=oms[sel]
        Tslope=1/np.polyfit(o,np.log(1+1/nv[sel]),1)[0]      # pendiente (def FP Ec.11)
        Textr=np.polyfit(o,Tw[sel],2)[-1]                    # T_0=lim omega->0 (Ec.12)
        print(f"  [{lo:.2f},{hi:.2f}] {sel.sum():5d} {Tslope/T0:13.4f} {Textr/T0:14.4f} {Textr/Tkbar:15.4f}")
    # convergencia en N (banda baja) en el omega mas bajo
    o0=oms[0]; a1=nL(w,o0,L,N); a2=nL(w,o0,L,int(N*1.5))
    print(f"  conv N en om/T0={xr[0]:.3f}: T_omega(n_L)  N={N}:{o0/np.log(1+1/a1)/T0:.4f}  "
          f"N={int(N*1.5)}:{o0/np.log(1+1/a2)/T0:.4f}  (x kloc)")
