"""
GATE P3.2-B2c — CIERRE del límite continuo de B2.3 (control pedido por ChatGPT, 2026-11-15).

NO toca el checkpoint congelado gate_P3_2B2_green.py. Lo COMPLEMENTA.

Motivo (dictamen ChatGPT sobre B2.3): la pendiente 1.07 prueba que el flujo matricial
F = dx·wᵀ[P_R,K_P]u converge al concomitante DISCRETO candidato [BF]_disc (BF con derivadas
finitas de los campos discretos). NO prueba convergencia al concomitante CONTINUO. Para validar
la correspondencia con Green en el continuo hay que:
  (1) construir una referencia CONTINUA INDEPENDIENTE del concomitante (otra discretización),
  (2) comprobar que F_matriz Y [BF]_disc convergen a ESA referencia, ajustando ambos errores
      por separado.

REFERENCIA CONTINUA (independiente, SIN MALLA — convolución con las funciones de Green de línea):
  G=1 homogéneo, κ²=2G=2, κ=√2. A=(2G+D)⁻¹=(2-∂²)⁻¹. Núcleos de Green de (κ²-∂²)ⁿ (Matérn/Yukawa):
      g1(t)=e^{-κ|t|}/(2κ)
      g2(t)=(1+κ|t|)e^{-κ|t|}/(4κ³)
      g3(t)=(3+3κ|t|+κ²t²)e^{-κ|t|}/(16κ⁵)
  Auxiliares como CONVOLUCIÓN (derivadas puestas sobre el FUENTE lorentziano, exacto por sympy, para
  que el núcleo quede SIEMPRE continuo => cuadratura sin discontinuidades, sin aliasing, sin amplificar
  frecuencias altas — el defecto del intento espectral previo):
      p  = g1*(-a'')     p'  = g1*(-a''')
      q1 = g1*a''''      q1' = g1*a⁽⁵⁾
      q2 = g2*a''''      q2' = g2*a⁽⁵⁾
      q3 = g3*a''''      q3' = g3*a⁽⁵⁾
  Derivadas 2ª/3ª por RECURSIÓN de la EDO (exacta, sin derivar el núcleo):
      (2-∂²)p=-a''   -> p''=2p+a''   , p'''=2p'+a'''
      (2-∂²)q1=a''''  -> q1''=2q1-a'''' , q1'''=2q1'-a⁽⁵⁾
      (2-∂²)q2=q1     -> q2''=2q2-q1    , q2'''=2q2'-q1'
      (2-∂²)q3=q2     -> q3''=2q3-q2    , q3'''=2q3'-q2'
  u,w y derivadas: ANALÍTICAS (lorentzianas). Cuadratura trapezoidal fina en y∈[-R,R]; convergencia
  verificada variando (R, paso h). Totalmente independiente de la malla FD del gate.

El concomitante se evalúa en x=±6 EXACTOS: la referencia continua se evalúa directamente en esos puntos
(sin malla propia), y la malla FD del gate los tiene como nodos exactos (N=401,801,1601,3201).

Expresión del concomitante (idéntica a la del checkpoint, con términos NO-locales):
  BF = ½ W(w,u)
     + c1[ W(w,p_u) + W(p_w,u) - W(p_w,p_u) ]
     + c2[ W2(w,q3u)+W2(q3w,u) - W(q1w,q3u) - W(q2w,q2u) - W(q3w,q1u) ]
  W(a,b)=a'b-ab' ; W2(a,b)=a b'''-a' b''+a'' b'-a''' b.
"""
import numpy as np
import sympy as sp
from gate_P3_2B2_green import build, chain, KP_mat, BF as BF_fd, c1, c2, G

# kappa derivada de G (no fijada a mano) para que la referencia continua use SIEMPRE el mismo operador
# que el discreto. Los nucleos de Green de linea g1,g2,g3 asumen G CONSTANTE (homogeneo): se exige.
assert np.isscalar(G) and abs(G-1.0) < 1e-15, f"B2c asume G=1 homogeneo; G={G}"
KAPPA = np.sqrt(2.0*G)   # kappa^2 = 2G

# ----------------------------------------------------------- derivadas analíticas de la lorentziana
def _lorentz_lambdas(x0):
    """Devuelve funciones d0..d5 (derivadas 0..5) de 4/((x-x0)^2+4), exactas vía sympy."""
    xs = sp.symbols('x', real=True)
    f = 4/((xs-x0)**2+4)
    return [sp.lambdify(xs, sp.diff(f, xs, m), 'numpy') for m in range(6)]

# ----------------------------------------------------- núcleos de Green de línea y sus convoluciones
def _g1(t, k=KAPPA): r=np.abs(t); return np.exp(-k*r)/(2*k)
def _g2(t, k=KAPPA): r=np.abs(t); return (1+k*r)*np.exp(-k*r)/(4*k**3)
def _g3(t, k=KAPPA): r=np.abs(t); return (3+3*k*r+(k*r)**2)*np.exp(-k*r)/(16*k**5)

def continuum_BF(R=50.0, h=2.5e-3, xa=-6.0, xb=6.0):
    """Concomitante continuo [BF]_{xa}^{xb} por convolución con las Green de línea (sin malla FD).
    h es el paso OBJETIVO; la malla es linspace exacta y la cuadratura usa el paso REAL dy."""
    n  = int(round(2*R/h))                           # nº de intervalos
    y  = np.linspace(-R, R, n+1)                     # malla exacta (último intervalo = dy exacto)
    dy = y[1]-y[0]                                    # paso real usado en la cuadratura
    du = _lorentz_lambdas(-3.0); dw = _lorentz_lambdas(5.0)

    def aux_at(d, x):
        """Devuelve (p,q1,q2,q3) cada uno como (val,d1,d2,d3) en el punto x, para la lorentziana d."""
        a2=d[2](y); a3=d[3](y); a4=d[4](y); a5=d[5](y)   # fuentes en la malla de cuadratura
        G1=_g1(x-y); G2=_g2(x-y); G3=_g3(x-y)
        trp=lambda f: np.trapezoid(f, dx=dy)
        p =trp(G1*(-a2)); p1 =trp(G1*(-a3))              # valor y 1ª derivada por convolución
        q1 =trp(G1*a4);   q11=trp(G1*a5)
        q2 =trp(G2*a4);   q21=trp(G2*a5)
        q3 =trp(G3*a4);   q31=trp(G3*a5)
        # derivadas 2ª/3ª por recursión de la EDO (2-∂²)·=fuente, con fuentes analíticas en x
        a2x=d[2](x); a3x=d[3](x); a4x=d[4](x); a5x=d[5](x)
        P =(p ,p1 , 2*p +a2x, 2*p1 +a3x)
        Q1=(q1,q11, 2*q1-a4x, 2*q11-a5x)
        Q2=(q2,q21, 2*q2-q1 , 2*q21-q11)
        Q3=(q3,q31, 2*q3-q2 , 2*q31-q21)
        return P,Q1,Q2,Q3

    def Wc(av,bv):   return av[1]*bv[0]-av[0]*bv[1]                       # a'b-ab'
    def W2c(av,bv):  return av[0]*bv[3]-av[1]*bv[2]+av[2]*bv[1]-av[3]*bv[0]

    def BF_at(x):
        u=(du[0](x),du[1](x),du[2](x),du[3](x)); w=(dw[0](x),dw[1](x),dw[2](x),dw[3](x))
        pu,q1u,q2u,q3u=aux_at(du,x); pw,q1w,q2w,q3w=aux_at(dw,x)
        r =0.5*Wc(w,u)
        r+=c1*(Wc(w,pu)+Wc(pw,u)-Wc(pw,pu))
        r+=c2*(W2c(w,q3u)+W2c(q3w,u)-Wc(q1w,q3u)-Wc(q2w,q2u)-Wc(q3w,q1u))
        return r
    return BF_at(xb)-BF_at(xa)

# ----------------------------------------------------------------------------- lado FD (del gate)
def fd_F_and_BFdisc(N, L=20.0):
    """Reproduce F_matriz (conmutador) y [BF]_disc del checkpoint, con ±6 como nodos EXACTOS."""
    x=np.linspace(-L,L,N); dx=x[1]-x[0]; D,D2,A=build(N,dx); KP=KP_mat(D,D2,A)
    w=1/(1+((x-5)/2)**2); u=1/(1+((x+3)/2)**2)
    pu,q1u,q2u,q3u=chain(u,D,D2,A); pw,q1w,q2w,q3w=chain(w,D,D2,A); aux=(pu,q1u,q2u,q3u,pw,q1w,q2w,q3w)
    ia=int(np.argmin(np.abs(x+6))); ib=int(np.argmin(np.abs(x-6)))
    assert abs(x[ia]+6)<1e-9 and abs(x[ib]-6)<1e-9, (x[ia],x[ib])   # exige nodo exacto
    bf=lambda i: BF_fd(w,u,*aux,i,dx,nonlocal_terms=True)
    BFdisc=bf(ib)-bf(ia)
    PR=np.diag((np.arange(N)>=ia)&(np.arange(N)<=ib)).astype(float)
    Fmat=dx*(w@(PR@KP-KP@PR)@u)
    return dx,Fmat,BFdisc

def run():
    print("GATE P3.2-B2c — cierre del límite continuo de B2.3\n"+"="*70)

    print("Referencia continua (convolución Green de línea): convergencia al refinar (R,h):")
    print(f"  {'R':>5} {'h':>9} {'[BF]_continuo':>18}")
    for R,h in [(40,5e-3),(50,2.5e-3),(60,1.25e-3),(70,6.25e-4)]:
        v=continuum_BF(R=R,h=h)
        print(f"  {R:>5} {h:9.2e} {v:18.10e}")

    # controles SEPARADOS (ChatGPT): aislar cuadratura (h) de truncamiento de dominio (R)
    print("\n  (a) Convergencia en h (R=70 FIJO):   -> error de CUADRATURA controlado")
    print(f"      {'h':>9} {'[BF]_continuo':>18}")
    prev=None
    for h in [5e-3,2.5e-3,1.25e-3,6.25e-4,3.125e-4]:
        v=continuum_BF(R=70,h=h)
        print(f"      {h:9.2e} {v:18.10e}"+(f"   Δ vs h previo={abs(v-prev):.2e}" if prev is not None else "")); prev=v
    print("\n  (b) Convergencia en R (h=6.25e-4 FIJO): -> truncamiento de dominio controlado")
    print(f"      {'R':>5} {'[BF]_continuo':>18}")
    prev=None
    for R in [30,40,50,60,70,90]:
        v=continuum_BF(R=R,h=6.25e-4)
        print(f"      {R:>5} {v:18.10e}"+(f"   Δ vs R previo={abs(v-prev):.2e}" if prev is not None else "")); prev=v

    BF_cont=continuum_BF(R=70,h=6.25e-4)
    print(f"\n  -> BF_continuo (R=70,h=6.25e-4) = {BF_cont:.10e}")

    print("\nConvergencia de F_matriz y [BF]_disc a la referencia CONTINUA (errores por separado):")
    print(f"  {'N':>5} {'dx':>8} {'F_matriz':>15} {'BF_disc':>15} "
          f"{'|F-BFcont|':>11} {'|BFdisc-BFcont|':>15} {'|F-BFdisc|':>11}")
    Ns=[401,801,1601,3201]
    dxs=[];eF=[];eBF=[];eFB=[]
    for N in Ns:
        dx,Fmat,BFdisc=fd_F_and_BFdisc(N)
        eF_=abs(Fmat-BF_cont); eBF_=abs(BFdisc-BF_cont); eFB_=abs(Fmat-BFdisc)
        dxs.append(dx);eF.append(eF_);eBF.append(eBF_);eFB.append(eFB_)
        print(f"  {N:>5} {dx:8.4f} {Fmat:15.6e} {BFdisc:15.6e} "
              f"{eF_:11.3e} {eBF_:15.3e} {eFB_:11.3e}")
    def slope(e): return np.polyfit(np.log(dxs),np.log(e),1)[0]
    print(f"\n  pendientes (orden en dx):")
    print(f"    F_matriz  -> BF_continuo : {slope(eF):.2f}")
    print(f"    BF_disc   -> BF_continuo : {slope(eBF):.2f}")
    print(f"    F_matriz  -> BF_disc     : {slope(eFB):.2f}   (ya en el checkpoint: ~1.07)")

    print("\nLectura:")
    print(" - Si AMBOS (F_matriz y BF_disc) -> BF_continuo con orden >0 definido => el flujo matricial")
    print("   converge al concomitante CONTINUO de Green, no sólo al discreto candidato. B2.3 cerrado.")
    print(" - BF_continuo es independiente (convolución Green de línea, sin malla FD) => comprobación válida.")
    print(" - Alcance: G=1 homogéneo, Dirichlet-fantasma, región [-6,6]. NO identidad SBP formal (B2.2).")

if __name__=='__main__':
    run()
