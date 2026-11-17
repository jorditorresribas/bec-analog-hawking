"""
TABLA AUTORITATIVA del benchmark FP (alta precision, 5 anchuras). Consolida:
 kappa_local, kappa_local/2pi, kappa_bar(sim,d=1.58)/2pi,
 T_upstream=T_0(n_L) y T_total=T_0(n_tot) (ambos = lim omega->0 de T_omega),
 y los cocientes clave T_upstream/(kappa_bar/2pi) [comparacion principal con FP].
Observable FP = n_L (canal u upstream, Ec.42). Guarda out/fp_benchmark.csv.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'solver'))  # layout del repo
import numpy as np, csv, mode_matching as MM
g=MM.g; n_u,mach_u,md=1.0,0.5,2.0
vLb=mach_u*np.sqrt(g*n_u); jc=n_u*vLb; n_d=(jc/(md*np.sqrt(g)))**(2/3)
def box(w): L=max(40.,8*w+20.); return L,int(round(L*240))
def geom(w):
    L,_=box(w); x=np.linspace(-L,L,80000); s=0.5*(1+np.tanh(x/w)); n=n_u+(n_d-n_u)*s
    c=np.sqrt(n); v=jc/n; f=c-v; i=np.where(np.sign(f[:-1])!=np.sign(f[1:]))[0][0]
    t=f[i]/(f[i]-f[i+1]); xH=x[i]+t*(x[i+1]-x[i]); cH=np.interp(xH,x,c)
    kK=abs(np.interp(xH,x,np.gradient(f,x)))
    Lam=2*cH**2; dxi=1.58*cH/(kK*Lam**2)**(1/3)
    kbar=abs(np.interp(xH+0.5*dxi,x,f)-np.interp(xH-0.5*dxi,x,f))/dxi
    return kK,kbar,cH
def occ(w,om,L,N):
    MM.w=w; MM.omega=om; x=np.linspace(-L,L,N); n,v,V=MM.bg(x,'horizon')
    Rin=[m for m in MM.modes_full(n[-1],v[-1],V[-1]) if m['prop'] and m['vg']<0]
    incR=[1.0 if m['norm']<0 else 0.0 for m in Rin]
    out,_=MM.solve_inject(L=L,N=N,mode='horizon',incL=0.0,incR=incR)
    nL=sum(a for(q,s,a)in out if q>0 and s=='L'); nRp=sum(a for(q,s,a)in out if q>0 and s=='R')
    pu=abs(sum(q*a for(q,s,a)in out)+1.0)
    return nL,nL+nRp,pu
rows=[]
print(f"{'w':>3} {'kappa_local':>11} {'kloc/2pi':>9} {'kbar/2pi':>9} {'T_up':>9} {'T_tot':>9} "
      f"{'T_up/kloc':>9} {'T_up/kbar':>10} {'T_tot/kloc':>10} {'pu_max':>8}")
for w in (2.0,3.0,4.0,6.0,8.0):
    kK,kbar,cH=geom(w); T0=kK/(2*np.pi); Tkbar=kbar/(2*np.pi); L,N=box(w)
    xr=np.array([0.05,0.065,0.08,0.10,0.12,0.15,0.18,0.22,0.27]); oms=xr*T0
    data=[occ(w,om,L,N) for om in oms]
    nL=np.array([d[0] for d in data]); nt=np.array([d[1] for d in data]); pu=np.array([d[2] for d in data])
    TwL=oms/np.log(1+1/nL); Twt=oms/np.log(1+1/nt)
    c2=np.polyfit(oms,TwL,2); Tup=c2[-1]; Ttot=np.polyfit(oms,Twt,2)[-1]   # lim omega->0 (FP Ec.12)
    # error de aproximar el limite: residuo RMS del ajuste cuadratico y variacion por orden
    rms=np.sqrt(np.mean((TwL-np.polyval(c2,oms))**2))
    Tup_lin=np.polyfit(oms,TwL,1)[-1]; Tup_cub=np.polyfit(oms,TwL,3)[-1]
    dord=max(abs(Tup-Tup_lin),abs(Tup-Tup_cub))/Tup      # dispersion relativa por orden
    print(f"{int(w):3d} {kK:11.6f} {T0:9.6f} {Tkbar:9.6f} {Tup:9.6f} {Ttot:9.6f} "
          f"{Tup/T0:9.4f} {Tup/Tkbar:10.4f} {Ttot/T0:10.4f} {pu.max():8.1e} {rms:8.1e} {dord:7.2%}")
    rows.append([w,kK,cH,T0,Tkbar,Tup,Ttot,Tup/T0,Tup/Tkbar,Ttot/T0,pu.max(),rms,dord,Tup_lin,Tup_cub])
with open('out/fp_benchmark.csv','w',newline='') as f:
    wr=csv.writer(f); wr.writerow(['w','kappa_local','cH','kloc_over_2pi','kbar_over_2pi',
        'T_upstream','T_total','Tup_over_kloc2pi','Tup_over_kbar2pi','Ttot_over_kloc2pi',
        'pseudo_u_max','fit_rms_resid','order_dispersion','T_up_linear','T_up_cubic']); wr.writerows(rows)
print("\nCSV: out/fp_benchmark.csv  (observable FP = T_upstream; T_total secundario)")
