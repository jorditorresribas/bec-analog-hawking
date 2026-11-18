"""
DIAGNOSTICO (no gate, sin commit, 2026-11-18) — <dn(x) dn(x')> DIRECTO desde la ordenacion del vacio, sin la formula
del articulo. No modifica el solver.

Con dphi = sum_pos [u a + v* a^dag] + [u_d2 a_d2^dag + v_d2* a_d2] (norma - : a y a^dag intercambiados) y
dn = sqrt(n)(dphi + dphi^dag), en el vacio-in:
   <dn(x) dn(x')> = sqrt(n n') * N_w * int dw  D(x,x';w),   D = sum_pos r(x) r*(x') + r_d2*(x) r_d2(x'),  r = u + v
(incluye el termino propio n delta(x-x')). Con modos |J| = 1 (los del codigo), int dx (|u|^2-|v|^2) = 2 pi delta(w-w'),
asi que N_w = 1/(2 pi). Algebra: Re D - P_codigo = K, con P_codigo = Re(B + B^dag) (receta del codigo, sin prefactor) y
K = Re sum_I s_I [u_I u_I'* - v_I v_I'*] (completitud BdG, s_I = signo de la norma), que con la base completa es
2 pi delta(x - x').
Pruebas:
 H  homogeneo Mach 0,5: n/(2 pi) Re D frente a la formula exacta n int dk/2pi S(k) e^{ik dx} (con termino propio).
 Z1 horizonte lam=1: identidad Re D = P_codigo + K punto a punto (comprueba el algebra en la implementacion).
 Z2 K se concentra en la diagonal al subir el corte en omega (es el termino propio) y la parte antisimetrica de D
    (Im D, = <[dn, dn']>/2i, que debe ser 0) tambien decae; la parte antisimetrica de la lectura literal, no.
Salida: out/diag_recati_direct_run.log
"""
import os, sys, io, time, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C
t0 = time.time(); g = C.g

# ---------------- H: homogeneo ----------------
L, N, NC = 40.0, 8000, 241
idx = np.linspace(0, N-1, NC).astype(int); x = np.linspace(-L, L, N)[idx]
OM = np.linspace(0.004, 1.2, 300); wts = np.gradient(OM); wts[[0, -1]] *= 0.5   # trapecio
D = np.zeros((NC, NC), complex); P = np.zeros((NC, NC)); kp = []; km = []
for om, wt in zip(OM, wts):
    for ch in ('u', 'd1'):
        r = C.scatter_field(L, N, om, ch, 2.0, 1.0, 'homog')
        (kp if ch == 'u' else km).append((r['Lin'] if ch == 'u' else r['Rin_pos'])['k'])
        u = r['u'][idx]; vv = r['v_mode'][idx]; rr = u + vv
        D += np.outer(rr, np.conj(rr))*wt
        B = np.outer(np.conj(vv), rr); P += np.real(B + np.conj(B.T))*wt
n0 = C.n_u; ic = NC//2; seps = np.arange(0, 40); dxs = x[ic+seps] - x[ic]
def exact(dx_, with_self):
    out = 0.0
    for a, b in ((min(kp), max(kp)), (min(km), max(km))):
        kk = np.linspace(a, b, 20001); e = kk**2/2; E = np.sqrt(e*(e+2*g*n0))
        out += np.trapz((e/E - (0 if with_self else 1))*np.cos(kk*dx_), kk)/(2*np.pi)
    return n0*out
dD = np.array([n0*np.real(D[ic, ic+s])/(2*np.pi) for s in seps]); eD = np.array([exact(d, True) for d in dxs])
dP = np.array([n0*P[ic, ic+s]/(2*np.pi) for s in seps]); eP = np.array([exact(d, False) for d in dxs])
fD = np.sum(dD*eD)/np.sum(eD*eD); fP = np.sum(dP*eP)/np.sum(eP*eP)
print("H homogeneo Mach 0,5 (u: k>0, d1: k<0, omega en [0.004, 1.2]):")
print(f"   directo  n/(2pi) Re D  frente a  n int S e^ikdx/2pi (con termino propio):   factor {fD:.6f}, "
      f"residuo {np.max(np.abs(dD-fD*eD))/np.max(np.abs(eD)):.1e}")
print(f"   receta del codigo con 1/(2pi) en vez de 1/(4pi)  frente a  n int (S-1) e^ikdx/2pi: factor {fP:.6f}, "
      f"residuo {np.max(np.abs(dP-fP*eP))/np.max(np.abs(eP)):.1e}   [{time.time()-t0:.0f} s]", flush=True)

# ---------------- Z: horizonte lam = 1 ----------------
L, N, NC = 30.0, 3600, 420
idx = np.linspace(0, N-1, NC).astype(int); xc = np.linspace(-L, L, N)[idx]
OMS = np.concatenate([np.arange(0.004, 0.47, 0.003), np.arange(0.47, 1.0, 0.006), np.arange(1.0, 2.0, 0.01),
                      np.arange(2.0, 4.0, 0.02), np.arange(4.0, 8.0001, 0.04)])
wts = np.gradient(OMS); wts[[0, -1]] *= 0.5; CUTS = (0.47, 1.0, 2.0, 4.0, 8.0)
Dacc = np.zeros((NC, NC), complex); Pacc = np.zeros((NC, NC)); Lacc = np.zeros((NC, NC)); Kacc = np.zeros((NC, NC))
snap = {}; ci = 0; ident = 0.0
for om, wt in zip(OMS, wts):
    while ci < len(CUTS) and om > CUTS[ci] + 1e-9:
        snap[CUTS[ci]] = (Dacc.copy(), Pacc.copy(), Lacc.copy(), Kacc.copy()); ci += 1
    pres = {}
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L, N, om, ch, 2.0, 1.0, 'horizon')
        if ch == 'd1' and r['Rin_pos'] is None: continue
        if ch == 'd2' and (r['Rin_neg'] is None or r['nchan'] != 2): continue
        pres[ch] = (r['u'][idx], r['v_mode'][idx])
    for ch, (u, vv) in pres.items():
        rr = u + vv
        if ch == 'd2':
            Dw = np.outer(np.conj(rr), rr); B = np.outer(np.conj(u), rr); s = -1.0
        else:
            Dw = np.outer(rr, np.conj(rr)); B = np.outer(np.conj(vv), rr); s = +1.0
        Pw = np.real(B + np.conj(B.T)); Kw = s*np.real(np.outer(u, np.conj(u)) - np.outer(vv, np.conj(vv)))
        ident = max(ident, np.abs(np.real(Dw) - Pw - Kw).max()/max(np.abs(Dw).max(), 1e-300))
        Dacc += Dw*wt; Pacc += Pw*wt; Lacc += 2*np.real(B)*wt; Kacc += Kw*wt
snap[8.0] = (Dacc.copy(), Pacc.copy(), Lacc.copy(), Kacc.copy())
print(f"\nZ horizonte lam=1, L={L}, N={N}, {len(OMS)} omegas hasta 8   [{time.time()-t0:.0f} s]")
print(f"Z1 identidad Re D = P_codigo + K, por omega y canal: max error relativo = {ident:.1e}")
dxm = np.abs(xc[:, None] - xc[None, :]); inner = (np.abs(xc)[:, None] <= 20) & (np.abs(xc)[None, :] <= 20)
off = inner & (dxm > 3.0)
print("Z2 (|x|,|x'| <= 20; 'fuera de diagonal' = |x-x'| > 3)")
print("   corte | K: max fuera diag / max en diag | Im D antisim.: max/max|Re D| | literal: max|L-L^T| / max|P|")
for c in CUTS:
    Dm, Pm, Lm, Km = snap[c]
    kdiag = np.abs(np.diag(Km)[np.abs(xc) <= 20]).max()
    print(f"   {c:5.2f} | {np.abs(Km[off]).max()/kdiag:.3e}                        | "
          f"{np.abs(np.imag(Dm)[inner]).max()/np.abs(np.real(Dm)[inner]).max():.3e}                    | "
          f"{np.abs((Lm-Lm.T)[inner]).max()/np.abs(Pm[inner]).max():.3e}")
print(f"fin [{time.time()-t0:.0f} s]")
