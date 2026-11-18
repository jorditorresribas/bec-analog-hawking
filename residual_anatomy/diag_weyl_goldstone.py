"""
DIAGNOSTICO (no gate) — por que el operador Weyl no local (gate0_pseudoherm) tiene residuo de Goldstone
~1.8e-3 ya en lam = 1, y por que sus cocientes en diag_goldstone (1.45-2.38) no son 1.
Rama diag-condicionamiento (2026-11-17). Solo lectura del operador: importa gate0_pseudoherm sin tocarlo.

Algebra (Phi = (s, -s), s = sqrt n, Qp fijo): (M Phi)_{u,v} = K s + Qp s -/+ i conv s, con Qp s = s''/2.
  => r_lam = r_1 + delta_lam,  r_1 = (K_1 s + s''/2) [+ conveccion],  delta_lam = (K_lam - K_1) s.
  diag_goldstone usa como prediccion ||delta_lam||; el cociente medido es ||r_1 + delta_lam|| / ||delta_lam||,
  que solo vale 1 si r_1 = 0.
Hipotesis H-W1: en lam = 1 el simbolo de Weyl es exactamente k^2/2 (sqrt((G + k^2/2)^2) - G), asi que K_1 es
  un -1/2 d^2 espectral de banda limitada sobre una caja FINITA. Aplicado a s, que NO decae en los bordes
  (s -> 1 y s -> sqrt(n_d)), el nucleo truncado deja una cola de borde que decae con la distancia al borde
  y no con dx.
Pruebas:
  P1  descomposicion lineal: ||r_lam - (r_1 + delta_lam)|| y angulo entre r_1 y delta_lam
  P2  r_1 frente a MARGIN (distancia al borde), L, N y Nk (cuadratura en k)
  P3  control: K_1 aplicado a una gaussiana (decae en los bordes) frente a -f''/2 analitico
  P4  control: misma s y mismo dx en un dominio doble (L=60), evaluada en |x|<=20; si el residuo baja, es
      efecto del borde de la caja.
Salida: out/diag_weyl_goldstone.csv
"""
import os, sys, io, csv, time, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, HERE)
with contextlib.redirect_stdout(io.StringIO()):
    import gate0_pseudoherm as G0M
W = 2.0; t0 = time.time(); rows = []


def sdd(x):
    """s'' analitico, s = sqrt(n), n tanh."""
    n_d = (G0M.jcur/(2.0*np.sqrt(G0M.g)))**(2/3); Dn = n_d - G0M.n_u
    t = np.tanh(x/W); n = G0M.n_u + Dn*0.5*(1+t)
    n1 = Dn*(1-t**2)/(2*W); n2 = -Dn*t*(1-t**2)/W**2
    s = np.sqrt(n); return n2/(2*s) - n1**2/(4*s**3), s


def r1_vec(x, K1):
    d2, s = sdd(x); return K1@s + 0.5*d2, s


def rel(vec, s, mk): return float(np.linalg.norm(vec[mk])/np.linalg.norm(s[mk]))


# ---- P1: descomposicion en la configuracion de diag_goldstone (L=30, N=400, MARGIN=10) ----
L, N, MARG = 30.0, 400, 10.0
x = np.linspace(-L, L, N); mk = np.abs(x) <= L - MARG
K1 = G0M.weyl_kinetic(x, 1.0, W, 'horizon'); r1, s = r1_vec(x, K1)
print(f"P1 (L={L}, N={N}, MARGIN={MARG}): ||r_1||/||s|| = {rel(r1, s, mk):.3e}  (diag_goldstone C lam=1: 1.794e-03 con conveccion)")
for lam in (0.25, 0.5, 2.0, 4.0):
    Kl = G0M.weyl_kinetic(x, lam, W, 'horizon'); dl = (Kl - K1)@s
    rl = Kl@s + 0.5*sdd(x)[0]
    lin = np.linalg.norm((rl - (r1 + dl))[mk])/np.linalg.norm(rl[mk])
    cosang = float(np.dot(r1[mk], dl[mk])/(np.linalg.norm(r1[mk])*np.linalg.norm(dl[mk])))
    ratio = rel(rl, s, mk)/rel(dl, s, mk)
    pred_ratio = np.sqrt(1 + (rel(r1, s, mk)/rel(dl, s, mk))**2 + 2*cosang*rel(r1, s, mk)/rel(dl, s, mk))
    print(f"   lam={lam}: ||delta||={rel(dl, s, mk):.3e}  ||r_lam||={rel(rl, s, mk):.3e}  cociente={ratio:.4f}"
          f"  (recompuesto con r_1 y angulo: {pred_ratio:.4f})  cos(r_1,delta)={cosang:+.3f}  no-linealidad={lin:.1e}")
    rows.append(dict(prueba='P1', lam=lam, L=L, N=N, margin=MARG, Nk=3*N, r1=rel(r1, s, mk), delta=rel(dl, s, mk),
                     r_lam=rel(rl, s, mk), cociente=ratio, cos_r1_delta=cosang))
print(f"   [{time.time()-t0:.0f} s]", flush=True)

# ---- P2: r_1 frente a MARGIN, L, N, Nk ----
print("P2 ||r_1||/||s|| (lam=1) frente a distancia al borde, L, N y Nk:")
for (L_, N_, Nk_) in ((30.0, 400, None), (30.0, 800, None), (30.0, 400, 4800), (60.0, 800, None)):
    xx = np.linspace(-L_, L_, N_); Kk = G0M.weyl_kinetic(xx, 1.0, W, 'horizon', Nk_); rr, ss = r1_vec(xx, Kk)
    out = []
    for marg in (5.0, 10.0, 15.0, 20.0, 25.0) + ((40.0, 50.0) if L_ > 30 else ()):
        m_ = np.abs(xx) <= L_ - marg
        out.append(f"d={marg:.0f}:{rel(rr, ss, m_):.1e}")
        rows.append(dict(prueba='P2', lam=1.0, L=L_, N=N_, margin=marg, Nk=(Nk_ or max(800, 3*N_)), r1=rel(rr, ss, m_)))
    # residuo puntual en el centro del dominio (lejos de ambos bordes)
    ic = np.argmin(np.abs(xx)); out.append(f"|r_1(x=0)|={abs(rr[ic]):.1e}")
    print(f"   L={L_:4.0f} N={N_:4d} dx={xx[1]-xx[0]:.3f} Nk={Nk_ or max(800, 3*N_)}: " + '  '.join(out) + f"   [{time.time()-t0:.0f} s]", flush=True)

# ---- P3: control con gaussiana ----
xx = np.linspace(-30.0, 30.0, 400); Kk = G0M.weyl_kinetic(xx, 1.0, W, 'horizon')
f = np.exp(-xx**2/8); fdd = (xx**2/16 - 0.25)*f
err = np.max(np.abs(Kk@f + 0.5*fdd))/np.max(np.abs(0.5*fdd))
print(f"P3 control gaussiana (decae en los bordes): max|K_1 f + f''/2| / max|f''/2| = {err:.1e}")
rows.append(dict(prueba='P3_gaussiana', lam=1.0, L=30.0, N=400, margin=0.0, Nk=1200, r1=err))

# ---- P4: misma s, pero K_1 calculado en un dominio DOBLE (s extendida con sus colas) y evaluado en |x|<=20 ----
xx2 = np.linspace(-60.0, 60.0, 799)                     # mismo dx que L=30, N=400
Kk2 = G0M.weyl_kinetic(xx2, 1.0, W, 'horizon'); rr2, ss2 = r1_vec(xx2, Kk2)
m_in = np.abs(xx2) <= 20.0
print(f"P4 dominio doble (L=60, mismo dx), evaluado en |x|<=20: ||r_1||/||s|| = {rel(rr2, ss2, m_in):.1e}"
      f"   (frente a L=30, MARGIN=10: {rows[0]['r1']:.1e})")
rows.append(dict(prueba='P4_dominio_doble', lam=1.0, L=60.0, N=799, margin=40.0, Nk=max(800, 3*799), r1=rel(rr2, ss2, m_in)))

keys = ['prueba', 'lam', 'L', 'N', 'margin', 'Nk', 'r1', 'delta', 'r_lam', 'cociente', 'cos_r1_delta']
with open(os.path.join(OUT, 'diag_weyl_goldstone.csv'), 'w', newline='') as fh:
    wr = csv.DictWriter(fh, fieldnames=keys); wr.writeheader()
    for r in rows: wr.writerow({k: r.get(k, '') for k in keys})
print(f"fin ({time.time()-t0:.0f} s)")
