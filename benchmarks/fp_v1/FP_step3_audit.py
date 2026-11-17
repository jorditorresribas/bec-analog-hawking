"""
PASO 3 - AUDITORIA de la matriz S (pasos 1-2 de ChatGPT), familia tanh, SIN tocar el
solver. Corrige la sobre-correccion de la 1a version de FP_STEP3 ("exceso = artefacto"):
 (1) enumerar TODOS los canales IN (ambos bordes) con su norma -> ¿cuantos eta=-1?
     Si hay UNO, la inyeccion unica da la ocupacion espontanea COMPLETA n_i=|S_{i,p}|^2.
 (2) verificar |alpha|^2 - n_tot = 1 (normalizacion Bogoliubov) -> n_tot=|beta_total|^2,
     luego T_eff,tot (=estimador H4) es la T de BALANCE DETALLADO, no un artefacto.
 (3) termalidad: ln(1+1/n)=omega/T lineal por el origen (Planck) para n_tot y n_L.
Conclusion: hay DOS temperaturas Planckianas legitimas (n_tot>kappa_local, n_L<kappa_local)
que CABALGAN kappa_local; cual mapea al observable de FP queda ABIERTO (ver FP_STEP3.md).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'solver'))  # layout del repo
import numpy as np, mode_matching as MM
g=MM.g; n_u,mach_u,md=1.0,0.5,2.0
vLb=mach_u*np.sqrt(g*n_u); jc=n_u*vLb; n_d=(jc/(md*np.sqrt(g)))**(2/3)
def box(w): L=max(40.,8*w+20.); return L,int(round(L*240))
def kloc(w):
    L,_=box(w); x=np.linspace(-L,L,60000); s=0.5*(1+np.tanh(x/w)); n=n_u+(n_d-n_u)*s
    c=np.sqrt(n); v=jc/n; f=c-v; i=np.where(np.sign(f[:-1])!=np.sign(f[1:]))[0][0]
    t=f[i]/(f[i]-f[i+1]); xH=x[i]+t*(x[i+1]-x[i]); return abs(np.interp(xH,x,np.gradient(f,x)))
print("="*92); print("PASO 3 - AUDITORIA S: canales IN, Bogoliubov |a|^2-n_tot=1, termalidad"); print("="*92)
for w in (2.0,4.0,8.0):
    T0=kloc(w)/(2*np.pi); L,N=box(w)
    MM.w=w; MM.omega=0.3*T0; x=np.linspace(-L,L,N); n,v,V=MM.bg(x,'horizon')
    mL=MM.modes_full(n[0],v[0],V[0]); mR=MM.modes_full(n[-1],v[-1],V[-1])
    Lin=[m for m in mL if m['prop'] and m['vg']>0]; Rin=[m for m in mR if m['prop'] and m['vg']<0]
    n_neg=sum(1 for m in Lin+Rin if m['norm']<0)
    print(f"\n=== w={w:.0f}  kappa_local/2pi={T0:.5f} ===")
    print(f"  canales IN: izq(vg>0) {[m['norm'] for m in Lin]}, der(vg<0) {[m['norm'] for m in Rin]} "
          f"-> #IN norma -1 = {n_neg} ({'UNO: inyeccion COMPLETA' if n_neg==1 else 'MAS DE UNO!'})")
    xr=np.array([0.05,0.08,0.12,0.18,0.27,0.40,0.60]); oms=xr*T0
    print(f"  {'om/kloc':>8} {'n_L(+1,L)':>10} {'n_R+(+1,R)':>11} {'|a|^2(-1,R)':>12} {'n_tot':>8} {'|a|^2-ntot':>11}")
    nL_l=[]; ntot_l=[]
    for om in oms:
        MM.omega=om; Rink=[m for m in MM.modes_full(n[-1],v[-1],V[-1]) if m['prop'] and m['vg']<0]
        incR=[1.0 if m['norm']<0 else 0.0 for m in Rink]
        out,_=MM.solve_inject(L=L,N=N,mode='horizon',incL=0.0,incR=incR)
        nL=sum(a for(q,s,a)in out if q>0 and s=='L'); nRp=sum(a for(q,s,a)in out if q>0 and s=='R')
        alpha2=sum(a for(q,s,a)in out if q<0); ntot=nL+nRp
        print(f"  {om/(T0*2*np.pi):8.4f} {nL:10.4f} {nRp:11.4f} {alpha2:12.4f} {ntot:8.4f} {alpha2-ntot:11.5f}")
        nL_l.append(nL); ntot_l.append(ntot)
    oms=np.array(oms); nL_l=np.array(nL_l); ntot_l=np.array(ntot_l)
    for name,nv in [('n_tot',ntot_l),('n_L  ',nL_l)]:
        y=np.log(1+1/nv); c1=np.polyfit(oms,y,1); res=np.sqrt(np.mean((y-np.polyval(c1,oms))**2))
        print(f"    termalidad {name}: T_slope/kloc={1/c1[0]/T0:.4f}  intercept={c1[1]:.2e}(~0) rms={res:.1e}")
print("\n-> UN solo IN de norma -1 (inyeccion completa); |a|^2-n_tot=1 (Bogoliubov);")
print("   n_tot y n_L ambos Planckianos a T distintas que CABALGAN kappa_local.")
print("   Cual es el observable de FP: ABIERTO (FP_STEP3.md §4-5).")
