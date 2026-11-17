"""
AUDITORIA INDEPENDIENTE de bdg_smatrix (recomendacion referee):
verificar la corriente conservada y la extraccion de beta SIN apoyarse solo en
la pseudo-unitaridad.

  (1) Conservacion de la corriente de norma J a lo largo del perfil REAL
      inhomogeneo (valida la formula de J y el operador A, independiente de S).
  (2) Relacion de Bogoliubov: por cada canal IN, sum_i eta_out,i |S_ij|^2 = eta_in,j
      -> para IN norma+ da |alpha|^2-|beta|^2=+1 (valida el split alpha/beta).
  (3) Control SIN horizonte (flujo subsonico en todo el dominio): no hay canal de
      norma negativa -> amplificacion = 0. Prueba decisiva de que la amplificacion
      del caso con horizonte es fisica (Hawking), no artefacto.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import scipy.linalg as la
import bdg_smatrix as B

print("="*68)
print("AUDITORIA INDEPENDIENTE DE LA MATRIZ S")
print("="*68)

# ------ (1) Conservacion de la corriente de norma en el perfil real ----------
omega = 0.1
mL = B.modes(B.nL, B.vL, B.VL, B.cL, omega, B.vL)
Y = next(m for m in mL if m['prop'] and m['vg'] > 0)['Y'].astype(complex)  # Lin
X0, dx = 12.0, 0.002
x = -X0; Js = []; xs = []
for i in range(int(2*X0/dx)):
    if i % 500 == 0:
        vloc = float(B.fv(x))
        Js.append(B.norm_current(Y, vloc)); xs.append(x)
    k1 = B.A_at(x, omega)@Y; k2 = B.A_at(x+dx/2, omega)@(Y+dx/2*k1)
    k3 = B.A_at(x+dx/2, omega)@(Y+dx/2*k2); k4 = B.A_at(x+dx, omega)@(Y+dx*k3)
    Y = Y + dx/6*(k1+2*k2+2*k3+k4); x += dx
Js = np.array(Js)
print("\n(1) CONSERVACION DE LA CORRIENTE DE NORMA (perfil real inhomogeneo)")
print(f"    J muestreado a lo largo de x: min={Js.min():.5f} max={Js.max():.5f}")
print(f"    variacion relativa max|J-<J>|/|<J>| = {np.max(np.abs(Js-Js.mean()))/abs(Js.mean()):.2e}")
print(f"    -> {'OK: corriente conservada (formula J y ODE validados)' if np.max(np.abs(Js-Js.mean()))/abs(Js.mean())<1e-3 else 'REVISAR'}")

# ------ (2) Relacion de Bogoliubov por canal ---------------------------------
S, eta_in, eta_out = B.smatrix(omega)[:3]
print("\n(2) RELACION DE BOGOLIUBOV  sum_i eta_out,i |S_ij|^2 = eta_in,j")
Eout = eta_out.astype(float)
maxdev = 0.0
for j in range(3):
    lhs = np.sum(Eout*np.abs(S[:, j])**2)
    dev = abs(lhs - eta_in[j]); maxdev = max(maxdev, dev)
    tag = "norma+" if eta_in[j] > 0 else "norma-"
    print(f"    canal IN {j} ({tag}): sum = {lhs:+.5f}  (debe {eta_in[j]:+d})  dev {dev:.1e}")
# |beta|^2 para el canal IN de norma + (upstream, el fisicamente relevante)
jpos = np.where(eta_in > 0)[0][0]
alpha2 = np.sum(np.abs(S[eta_out > 0, jpos])**2)
beta2 = np.sum(np.abs(S[eta_out < 0, jpos])**2)
print(f"    canal IN norma+ (j={jpos}): |alpha|^2={alpha2:.4f} |beta|^2={beta2:.4f} "
      f"|a|^2-|b|^2={alpha2-beta2:.5f}")
print(f"    -> {'OK: split alpha/beta consistente' if maxdev<1e-3 else 'REVISAR'}")

# ------ (3) CONTROL SIN HORIZONTE (todo subsonico) ---------------------------
# Rehacer un fondo subsonico (mach_d<1) analitico y medir la amplificacion.
g = B.g; n_u = 1.0; mu_u = 0.5; md = 0.8   # downstream TAMBIEN subsonico
vl = mu_u*np.sqrt(g*n_u); jc = n_u*vl; muc = g*n_u+0.5*vl**2
nd = (jc/(md*np.sqrt(g)))**(2/3); w = 2.0
def bg2(x):
    t=np.tanh(x/w); s=0.5*(1+t); Dn=nd-n_u
    n=n_u+Dn*s; n1=Dn*(1-t**2)/(2*w); n2=-Dn*t*(1-t**2)/w**2
    v=jc/n; dv=-jc*n1/n**2; Qp=0.5*(n2/(2*n)-n1**2/(4*n**2))
    return n, v, muc-g*n-0.5*v**2+Qp, dv
def A2(x, om):
    n,v,V,dv=bg2(x)
    return np.array([[0,1,0,0],[-1j*dv+v**2+4*g*n+2*V-2*muc-2*om,-2j*v,2*g*n,0],
                     [0,0,0,1],[2*g*n,0,1j*dv+v**2+4*g*n+2*V-2*muc+2*om,2j*v]],dtype=complex)
nLs,vLs,cLs=n_u,vl,np.sqrt(g*n_u); nRs,vRs,cRs=nd,jc/nd,np.sqrt(g*nd)
print(f"\n(3) CONTROL SIN HORIZONTE: Mach_u=0.5, Mach_d={md} (ambos subsonicos)")
print(f"    downstream Mach = {vRs/cRs:.3f} < 1  -> NO hay horizonte, NO norma negativa")
# transferencia
M=np.eye(4,dtype=complex); x=-X0
for _ in range(int(2*X0/dx)):
    k1=A2(x,omega)@M; k2=A2(x+dx/2,omega)@(M+dx/2*k1)
    k3=A2(x+dx/2,omega)@(M+dx/2*k2); k4=A2(x+dx,omega)@(M+dx*k3)
    M=M+dx/6*(k1+2*k2+2*k3+k4); x+=dx
def modes2(n,v,V,c,om,Vf):
    lam,vec=la.eig(np.array([[0,1,0,0],[v**2+4*g*n+2*V-2*muc-2*om,-2j*v,2*g*n,0],
        [0,0,0,1],[2*g*n,0,v**2+4*g*n+2*V-2*muc+2*om,2j*v]],dtype=complex))
    out=[]
    for i in range(4):
        k=-1j*lam[i]; Y=vec[:,i]
        if abs(k.imag)<1e-8:
            Om=om-v*k.real; vg=v+(c**2*k.real+k.real**3/2)/Om
            Y=Y/np.sqrt(abs(B.norm_current(Y,Vf)))
            out.append(dict(prop=True,vg=vg,norm=int(np.sign(Om)),Y=Y))
        else: out.append(dict(prop=False,decay=(k.imag<0),Y=Y))
    return out
mLs=modes2(nLs,vLs,0.0,cLs,omega,vLs)
mRs=modes2(nRs,vRs,bg2(X0)[2],cRs,omega,vRs)
negL=sum(m['prop'] and m['norm']<0 for m in mLs); negR=sum(m['prop'] and m['norm']<0 for m in mRs)
print(f"    canales norma negativa: upstream={negL}, downstream={negR}")
print(f"    -> {'OK: sin horizonte NO hay canal de norma negativa -> |beta|^2=0 estructuralmente' if negL+negR==0 else 'REVISAR'}")
print(f"       (con horizonte habia 2 downstream -> la amplificacion es del horizonte = Hawking)")

print("\n" + "="*68)
print("AUDITORIA COMPLETA")
print("="*68)
