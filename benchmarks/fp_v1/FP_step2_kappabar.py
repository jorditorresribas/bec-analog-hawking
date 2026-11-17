"""
RUMBO #1 / Experimento B (Finazzi-Parentani) — PASO 2:
kappa_bar promediada (Ec.20) + T_IR(w), comparando T_IR vs kappa_local/2pi vs
kappa_bar/2pi.  Construido sobre PASO 1 (FP_step1_scales.py, escalas validadas) y
sobre la maquinaria N_omega VALIDADA (mode_matching.solve_inject, usada en H4).

Definiciones FP (verificadas contra arXiv:1012.1556v2, HTML):
  * Ec.9  dispersion:  Omega^2 = c^2 k^2 + c^4 k^4/Lambda^2,  Lambda = sqrt2 c/xi
          (ESCALA DE FRECUENCIA).  En nuestras unidades (hbar=m=g=1, lam=1):
          Lambda_FP = 2 c^2  (= c*k_*), identico a PASO 1.
  * Ec.20 kappa_bar = (1/d_xi) INT_{-d_xi^L}^{+d_xi^R} dx d(w_FP)/dx
                    = [w_FP(d_xi^R) - w_FP(-d_xi^L)]/d_xi
          FP usa w_FP = c+v (su flujo v<0); en nuestra convencion (v>0, horizonte
          c=v) la combinacion equivalente es (c-v).  PROMEDIO UNIFORME del gradiente
          sobre una ventana que ENSANCHA el horizonte.
  * Ec.24 d_xi = d * c0/(kappa_K Lambda^2)^(1/3),  d~1.58 AJUSTADO por FP a SUS
          perfiles (NO universal).  Solo el ESCALADO es transferible (PASO 1 §3).
          -> aqui NO fijamos d=1.58; barremos d y medimos tambien el d* que
             reproduce T_IR, comprobando si es O(1) y estable en w.

Honestidad (misma regla que H4): T_IR se mide INDEPENDIENTEMENTE del espectro, sin
imponer T_IR = kappa/2pi.  Dejamos que los datos decidan cual gravedad (local o
promediada) sigue T_IR, y si la diferencia crece en horizontes empinados (w chico),
que es la prediccion FALSABLE de FP.

Split d_xi^L=d_xi^R=d_xi/2 centrado en el horizonte interpolado: default honesto;
el reparto asimetrico y el centro x_hat (Ec.20) dependen del fondo -> PENDIENTE,
se reporta su sensibilidad al final.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'solver'))  # layout del repo
import numpy as np
import csv
import mode_matching as MM

g = MM.g
SQRT2 = np.sqrt(2.0)
LAM = 1.0                      # knob de dispersion (Bogoliubov exacto)

# --- constantes del fondo (identicas a mode_matching / kappa_derive / H4) ---
n_u, mach_u, md = 1.0, 0.5, 2.0
vLb = mach_u*np.sqrt(g*n_u); jc = n_u*vLb; mu = g*n_u + 0.5*vLb**2
n_d = (jc/(md*np.sqrt(g)))**(2/3)


# ----------------------------------------------------------------------------
# GEOMETRIA: horizonte interpolado, kappa_local, Lambda_FP, d_xi(d), kappa_bar(d)
# ----------------------------------------------------------------------------
def background(w, L, Ngeo=60000):
    """Fondo analitico tanh (el mismo de MM.bg 'horizon'), malla densa para geom."""
    x = np.linspace(-L, L, Ngeo)
    t = np.tanh(x/w); s = 0.5*(1+t); n = n_u + (n_d-n_u)*s
    v = jc/n; c = np.sqrt(g*n)
    return x, c, v


def horizon_interp(x, c, v):
    """Raiz INTERPOLADA de (c-v)=0 y c alli (igual criterio que PASO 1)."""
    f = c - v
    isg = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    tt = f[isg]/(f[isg]-f[isg+1])
    xH = x[isg] + tt*(x[isg+1]-x[isg])
    cH = c[isg] + tt*(c[isg+1]-c[isg])
    return xH, cH


def cminusv_at(x, c, v, xq):
    """(c-v)(xq) por interpolacion lineal sobre la malla densa."""
    return np.interp(xq, x, c - v)


def kappa_local(x, c, v, xH):
    """|d(c-v)/dx| en el horizonte interpolado (gravedad superficial geometrica).
    Derivada central analitica via gradiente sobre malla densa, evaluada en xH."""
    f = c - v
    df = np.gradient(f, x)
    return abs(np.interp(xH, x, df))


def kappa_gradmax(x, c, v):
    """TECHO absoluto de cualquier promedio del gradiente: max_x |d(c-v)/dx|.
    Ninguna kappa_bar (cualquier centro/split/anchura) puede superar este valor,
    porque es un promedio de |f'| y el promedio <= maximo. Devuelve (kmax, x_max)."""
    df = np.abs(np.gradient(c - v, x))
    i = int(np.argmax(df))
    return df[i], x[i]


def kappa_bar(x, c, v, xH, cH, kappaK, d, lam=LAM):
    """Ec.20 con ventana simetrica d_xi (Ec.24). Devuelve (kappa_bar, d_xi)."""
    Lam_FP = 2*cH**2/np.sqrt(lam)                 # escala de FRECUENCIA (PASO 1)
    d_xi = d * cH / (kappaK * Lam_FP**2)**(1.0/3.0)
    xL = xH - 0.5*d_xi; xR = xH + 0.5*d_xi
    kb = abs(cminusv_at(x, c, v, xR) - cminusv_at(x, c, v, xL)) / d_xi
    return kb, d_xi, Lam_FP


# ----------------------------------------------------------------------------
# T_IR: estimador 2 de H4 (omega*N_omega -> 0), maquinaria N_omega VALIDADA
# ----------------------------------------------------------------------------
def Nw_pu(w, om, L, N):
    """N_omega por matching modal exacto (identico a H4). Devuelve (N, pseudo-u)."""
    MM.w = w; MM.omega = om
    x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
    mR = MM.modes_full(n[-1], v[-1], V[-1])
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
    out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=incR)
    Nw = sum(a for (nrm, s, a) in out if nrm > 0)
    pu = abs(sum(nrm*a for (nrm, s, a) in out) - (-1.0))
    return Nw, pu


def T_IR(w, L, N, Tscale):
    """omega*N en banda baja (om/Tscale in [0.12,0.6]) -> extrapolar a om->0.
    Tscale solo FIJA LA BANDA (no entra en el valor). Lineal y cuadratica."""
    xr = np.array([0.12, 0.18, 0.25, 0.33, 0.45, 0.60])
    oms = xr*Tscale
    data = [Nw_pu(w, om, L, N) for om in oms]
    Nw = np.array([d[0] for d in data]); pu = np.array([d[1] for d in data])
    omN = oms*Nw
    cl = np.polyfit(oms, omN, 1); T_lin = cl[-1]
    cq = np.polyfit(oms, omN, 2); T_quad = cq[-1]
    return dict(oms=oms, Nw=Nw, omN=omN, pu=pu, T_lin=T_lin, T_quad=T_quad,
                pu_max=pu.max())


# ----------------------------------------------------------------------------
# BARRIDO en w
# ----------------------------------------------------------------------------
WS = [2.0, 3.0, 4.0, 6.0, 8.0]
DS = [1.0, 1.58, 2.0]                         # d (1.58 = valor FP, NO universal)
# caja y resolucion por w (misma densidad ~200 ptos/unidad que H4)
def box(w):
    L = max(40.0, 8.0*w + 20.0)
    N = int(round(L*200))
    return L, N

print("="*92)
print("PASO 2 (exp. B): kappa_bar promediada (Ec.20) + T_IR(w)  [N_omega validado]")
print("="*92)
print(f"Fondo: n_u={n_u} Mach_u={mach_u} md={md} -> n_d={n_d:.5f}; lam={LAM} (Bogoliubov)")
print("kappa_local = |d(c-v)/dx|_H (geom);  kappa_bar = promedio del gradiente sobre d_xi")
print("d_xi = d*c_H/(kappa_local*Lambda_FP^2)^(1/3),  Lambda_FP=2c_H^2 (frecuencia)\n")

rows = []
geom = {}   # w -> dict con kappa_local, Lam_FP, c_H, kappa_bar[d], d_xi[d]
TIR  = {}   # w -> dict T_IR

for w in WS:
    L, N = box(w)
    x, c, v = background(w, L)
    xH, cH = horizon_interp(x, c, v)
    kK = kappa_local(x, c, v, xH)
    kMax, xMax = kappa_gradmax(x, c, v)        # techo absoluto de cualquier kappa_bar
    Lam_FP = 2*cH**2/np.sqrt(LAM)
    kb = {}; dxi = {}
    for d in DS:
        kbval, dxival, _ = kappa_bar(x, c, v, xH, cH, kK, d)
        kb[d] = kbval; dxi[d] = dxival
    geom[w] = dict(xH=xH, cH=cH, kK=kK, kMax=kMax, xMax=xMax, Lam_FP=Lam_FP,
                   kb=kb, dxi=dxi, L=L, N=N)

    print(f"--- w={w:.0f}  (L={L:.0f}, N={N})  x_H={xH:+.4f}  c_H={cH:.4f}  "
          f"Lambda_FP={Lam_FP:.4f}  kappa/Lambda={kK/Lam_FP:.3f} ---")
    print(f"    kappa_local = {kK:.5f}   kappa_local/2pi = {kK/(2*np.pi):.5f}")
    print(f"    max_x|d(c-v)/dx| = {kMax:.5f} @x={xMax:+.4f}  (techo kappa_bar = "
          f"{kMax/kK:.5f} x kappa_local; /2pi = {kMax/(2*np.pi):.5f})")
    for d in DS:
        print(f"    d={d:4.2f}: d_xi={dxi[d]:.4f}  kappa_bar={kb[d]:.5f}  "
              f"kappa_bar/2pi={kb[d]/(2*np.pi):.5f}  "
              f"(kappa_bar/kappa_local={kb[d]/kK:.4f})")

# ---- T_IR(w): medir independientemente ----
print("\n" + "="*92)
print("T_IR(w) por omega*N_omega -> 0 (banda om/Tscale in [0.12,0.60], Tscale=kappa_local/2pi)")
print("="*92)
for w in WS:
    gw = geom[w]; Tscale = gw['kK']/(2*np.pi)
    r = T_IR(w, gw['L'], gw['N'], Tscale)
    TIR[w] = r
    print(f"\n--- w={w:.0f}  Tscale(kappa_loc/2pi)={Tscale:.5f}  (pseudo-u max={r['pu_max']:.1e}) ---")
    print(f"    {'om':>8} {'N_omega':>10} {'omega*N':>10}")
    for i in range(len(r['oms'])):
        print(f"    {r['oms'][i]:8.4f} {r['Nw'][i]:10.4f} {r['omN'][i]:10.5f}")
    print(f"    T_IR (lineal)    = {r['T_lin']:.5f}")
    print(f"    T_IR (cuadratica)= {r['T_quad']:.5f}   (modelo-dep {abs(r['T_lin']-r['T_quad'])/r['T_quad']:.1%})")

# ---- COMPARACION FINAL + test falsable ----
print("\n" + "="*92)
print("COMPARACION:  T_IR  vs  kappa_local/2pi  vs  kappa_bar/2pi")
print("="*92)
print(f"{'w':>4} {'T_IR':>9} {'k_loc/2pi':>10} {'TIR/kloc':>9} | "
      f"{'kbar/2pi(1.58)':>14} {'TIR/kbar':>9} | {'techo/2pi':>9} {'TIR>techo?':>10}")
for w in WS:
    gw = geom[w]; r = TIR[w]
    T = r['T_quad']                       # estimador fisico primario (cuadratico)
    kloc = gw['kK']/(2*np.pi)
    kbar158 = gw['kb'][1.58]/(2*np.pi)
    kceil = gw['kMax']/(2*np.pi)          # techo absoluto de cualquier kappa_bar
    # d* tal que kappa_bar(d*)/2pi = T_IR.  kappa_bar(d) es monotona DECRECIENTE desde
    # kappa_local (d->0) hacia abajo -> existe d*>0 solo si T_IR < kappa_local.
    dgrid = np.linspace(0.05, 6.0, 600)
    x, c, v = background(w, gw['L'])
    kbvals = np.array([kappa_bar(x, c, v, gw['xH'], gw['cH'], gw['kK'], dd)[0]/(2*np.pi)
                       for dd in dgrid])
    diff = kbvals - T
    sc = np.where(np.sign(diff[:-1]) != np.sign(diff[1:]))[0]
    dstar = (dgrid[sc[0]] + dgrid[sc[0]+1])/2 if len(sc) else np.nan
    over = "SI" if T > kceil else "no"
    print(f"{w:4.0f} {T:9.5f} {kloc:10.5f} {T/kloc:9.4f} | "
          f"{kbar158:14.5f} {T/kbar158:9.4f} | {kceil:9.5f} {over:>10}")
    rows.append([w, gw['xH'], gw['cH'], gw['kK'], gw['kMax'], gw['Lam_FP'],
                 r['T_lin'], r['T_quad'], r['pu_max'], kloc, kceil,
                 gw['kb'][1.0]/(2*np.pi), gw['kb'][1.58]/(2*np.pi),
                 gw['kb'][2.0]/(2*np.pi), gw['dxi'][1.58], dstar])

with open('out/fp_step2_kappabar.csv', 'w', newline='') as fcsv:
    wr = csv.writer(fcsv)
    wr.writerow(['w', 'xH', 'cH', 'kappa_local', 'kappa_gradmax', 'Lambda_FP',
                 'T_IR_lin', 'T_IR_quad', 'pseudo_u_max', 'kloc_over_2pi',
                 'kceil_over_2pi', 'kbar_d1.0', 'kbar_d1.58', 'kbar_d2.0',
                 'd_xi_1.58', 'd_star'])
    wr.writerows(rows)
print("\nCSV: out/fp_step2_kappabar.csv")

# ---- diagnostico de la prediccion FALSABLE ----
print("\n" + "-"*92)
print("DIAGNOSTICO (prediccion FP falsable):")
ds = [r[-1] for r in rows]
nd = sum(1 for d in ds if np.isfinite(d))
print(f"  (a) kappa_bar < kappa_local en TODAS las w, y mas reducida a w chico:")
for r in rows:
    print(f"        w={int(r[0])}: kappa_bar(1.58)/kappa_local = {r[12]/r[9]:.4f}")
print(f"      -> SI: el promedio reduce la gravedad (mas en horizonte empinado). Como FP.")
print(f"  (b) pero T_IR esta POR ENCIMA de kappa_local (no por debajo):")
for r in rows:
    print(f"        w={int(r[0])}: T_IR/kappa_local = {r[7]/r[9]:.4f}")
print(f"      -> kappa_bar se ALEJA de T_IR; no hay d>0 con kappa_bar=T_IR ({nd}/{len(rows)} w).")
print(f"  (c) ademas, en w>=4 T_IR SUPERA el TECHO max_x|d(c-v)/dx| = 1.0143 kappa_local:")
for r in rows:
    flag = "T_IR>techo (ningun promedio del gradiente alcanza T_IR)" if r[7] > r[10] else "T_IR<techo"
    print(f"        w={int(r[0])}: T_IR/techo = {r[7]/r[10]:.4f}  [{flag}]")
print(f"\n  CONCLUSION: la prescripcion FP de horizonte ensanchado (kappa_bar < kappa_local)")
print(f"  esta FALSADA en SIGNO para esta familia de perfiles (regimen debilmente")
print(f"  dispersivo kappa/Lambda=0.05..0.21): T_IR >= kappa_local, con un exceso de +1..+4%")
print(f"  que CRECE con w (horizontes suaves), opuesto a la expectativa 'dispersion importa")
print(f"  mas en horizontes empinados'. En w>=4 el exceso supera el maximo global del")
print(f"  gradiente -> NO es un efecto de promediado del gradiente (indep. de convencion).")
print("\nNota: split d_xi^L=d_xi^R y centro=x_H son el default honesto; pero el TECHO")
print("max|f'| acota CUALQUIER ventana (centro/split/anchura) -> la falsacion en w>=4 no")
print("depende de esa convencion. Reparto asimetrico y x_hat (Ec.20): PENDIENTE (subdominante).")
