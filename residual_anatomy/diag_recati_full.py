"""
DIAGNOSTICO (no gate, sin commit, 2026-11-18) — seguimiento de diag_recati_norm.py. No modifica el solver.

P1  Base de modos en la prueba homogenea (Mach 0,5): para cada omega, que canales entrantes existen, su signo de k y
    de norma, y la identidad por modo  (|U|^2-|V|^2)|v_g| = 1  y  V(U+V)|v_g| = (S(k)-1)/2. Comprueba que cada k
    aparece una sola vez (omega(k) = v k + E_k monotona en cada rama) -> no se cuenta media base ni base doble.
P2  Horizonte lam = 1 (tanh w = 2): integracion en omega INCLUYENDO omega > omega_max (alli solo existen u y d1).
    Se guarda el corchete complejo B integrado por canal y por tramo de omega, y se mide la asimetria de la lectura
    literal 2 Re B frente a la del codigo Re(B + B^dag) al aumentar el corte en omega.
P3  Pendiente de la cresta de §3 con el metodo ORIGINAL (registro de la sesion 2026-11-14: para cada x upstream,
    argmin_x' C_d2 en x' in (0.5, 14); mediana de x'/x donde C_d2 < -1e-3), en la configuracion original
    (L=30, N=3600, 140 omegas en [0.02, 0.45], Hann 15 %, ncoarse 420), para: receta del codigo, lectura literal,
    umbral escalado (la amplitud correcta es 2x), y C total con la integracion completa.
Salida: out/diag_recati_full.npz, log en out/diag_recati_full_run.log
"""
import os, sys, io, time, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C
t0 = time.time(); g = C.g


def stamp(): return f"[{time.time()-t0:.0f} s]"


# ======================= P1 =======================
print("P1 base de modos, condensado homogeneo Mach 0,5 (los mismos ajustes que diag_recati_norm B):")
print("   omega   canal  k        Om=w-vk  norma  vg       (|U|^2-|V|^2)|vg|  V(U+V)|vg|   (S-1)/2      Rin_neg")
n0 = C.n_u; v0 = C.vLb
for om in (0.004, 0.05, 0.3, 0.8, 1.2):
    for ch in ('u', 'd1'):
        r = C.scatter_field(40.0, 8000, om, ch, 2.0, 1.0, 'homog')
        m = r['Lin'] if ch == 'u' else r['Rin_pos']
        k = m['k']; U, V = m['U'], m['V']; vg = m['vg']
        eps = k*k/2; E = np.sqrt(eps*(eps + 2*g*n0)); Sm1h = (eps/E - 1)/2
        print(f"   {om:5.3f}  {ch:4s}  {k:+.4f}  {om-v0*k:+.4f}  {m['norm']:+d}    {vg:+.4f}  "
              f"{(abs(U)**2-abs(V)**2)*abs(vg):.6f}           {np.real(np.conj(V)*(U+V))*abs(vg):+.6f}    {Sm1h:+.6f}   "
              f"{'ninguno' if r['Rin_neg'] is None else 'EXISTE'}")
ks = np.linspace(0.001, 6, 4000)
mono_p = np.all(np.diff(v0*ks + np.sqrt(ks**2/2*(ks**2/2+2*g*n0))) > 0)
mono_m = np.all(np.diff(-v0*ks + np.sqrt(ks**2/2*(ks**2/2+2*g*n0))) > 0)
print(f"   omega(k) = v k + E_k monotona en k>0: {mono_p};  en k<0 (en |k|): {mono_m}  -> cada k aparece una sola vez")
print(f"   {stamp()}", flush=True)

# ======================= P2 =======================
L, N, NC = 30.0, 3600, 420
idx = np.linspace(0, N-1, NC).astype(int); xc = np.linspace(-L, L, N)[idx]
OMS = np.concatenate([np.arange(0.004, 0.47, 0.003), np.arange(0.47, 1.0, 0.006), np.arange(1.0, 2.0, 0.01),
                      np.arange(2.0, 4.0, 0.02), np.arange(4.0, 8.0001, 0.04)])
wts = np.gradient(OMS)                                   # pesos de trapecio (malla no uniforme)
CUTS = (1.0, 2.0, 4.0, 8.0)
segs = ['win'] + [f'hi{c:g}' for c in CUTS]
acc = {(s, ch): np.zeros((NC, NC), complex) for s in segs for ch in ('u', 'd1', 'd2')}
omax_seen = 0.0; nch_hi = set()
print(f"\nP2 horizonte lam=1, L={L}, N={N}, {len(OMS)} omegas en [{OMS[0]}, {OMS[-1]:.2f}] ...", flush=True)
for om, wt in zip(OMS, wts):
    present = {}
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L, N, om, ch, 2.0, 1.0, 'horizon')
        if ch == 'd1' and r['Rin_pos'] is None: continue
        if ch == 'd2' and r['Rin_neg'] is None: continue
        present[ch] = r
    has_d2 = 'd2' in present and present['d2']['nchan'] == 2
    if has_d2:
        seg = 'win'; omax_seen = om
    else:
        seg = next(f'hi{c:g}' for c in CUTS if om <= c + 1e-9)
        nch_hi.add(tuple(sorted(present)))
    for ch, r in present.items():
        if ch == 'd2' and not has_d2: continue
        u = r['u'][idx]; vv = r['v_mode'][idx]
        a = np.conj(u) if ch == 'd2' else np.conj(vv)
        acc[(seg, ch)] += np.outer(a, u+vv)*wt
print(f"   ventana de Hawking hasta omega = {omax_seen:.3f}; canales por encima: {sorted(nch_hi)}  {stamp()}")
mk = np.abs(xc) <= 20; M = np.ix_(mk, mk)
pref = np.sqrt(np.outer(*(2*[C.bg(np.linspace(-L, L, N), 2.0, 'horizon')[0][idx]])))/(4*np.pi)
def total(upto):
    S = list(segs[:1 + CUTS.index(upto) + 1]) if upto else ['win']
    return sum(acc[(s, ch)] for s in S for ch in ('u', 'd1', 'd2'))
print("   corte en omega | asim. literal max|L-L^T|/max|L| | asim. literal rms | dif. literal-codigo (rel. max)")
for upto in (None,) + CUTS:
    B = total(upto); lit = 2*np.real(B); cod = np.real(B + np.conj(B.T))
    a_max = np.abs(lit - lit.T)[M].max()/np.abs(lit[M]).max()
    a_rms = np.sqrt(np.mean((lit - lit.T)[M]**2))/np.sqrt(np.mean(lit[M]**2))
    dlc = np.abs(lit - cod)[M].max()/np.abs(cod[M]).max()
    print(f"   {('ventana %.3f' % omax_seen) if upto is None else ('%.0f' % upto):>14} | {a_max:.3e} | {a_rms:.3e} | {dlc:.3e}")
# asimetria por canal en la ventana (d2 solo existe ahi)
for ch in ('u', 'd1', 'd2'):
    B = acc[('win', ch)]; lit = 2*np.real(B)
    print(f"   canal {ch:3s} (ventana): asim. literal max = {np.abs(lit-lit.T)[M].max()/np.abs(lit[M]).max():.3e}")
Bhi = sum(acc[(s, ch)] for s in segs[1:] for ch in ('u', 'd1')); lit = 2*np.real(Bhi)
print(f"   solo omega > omega_max (u+d1): asim. literal max = {np.abs(lit-lit.T)[M].max()/np.abs(lit[M]).max():.3e}  {stamp()}", flush=True)
np.savez(os.path.join(OUT, 'diag_recati_full.npz'), xc=xc, OMS=OMS, omax=omax_seen,
         **{f'{s}_{ch}': acc[(s, ch)] for s in segs for ch in ('u', 'd1', 'd2')})

# ======================= P3 =======================
print("\nP3 pendiente de la cresta (metodo original), configuracion original de §3:")
omegas = np.linspace(0.02, 0.45, 140); dw = omegas[1]-omegas[0]; nn = len(omegas)
win = np.clip(0.5*(1-np.cos(np.pi*np.minimum(np.arange(nn), nn-1-np.arange(nn))/(0.15*nn))), 0, 1)
Bd2 = np.zeros((NC, NC), complex); Bud = np.zeros((NC, NC), complex)
for iw, om in enumerate(omegas):
    Mm = {}
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L, N, om, ch, 2.0, 1.0, 'horizon')
        if r['Rin_neg'] is None or r['nchan'] != 2: Mm = None; break
        Mm[ch] = dict(u=r['u'][idx], v=r['v_mode'][idx], r=r['u'][idx]+r['v_mode'][idx])
    if Mm is None: continue
    Bud += (np.outer(np.conj(Mm['u']['v']), Mm['u']['r']) + np.outer(np.conj(Mm['d1']['v']), Mm['d1']['r']))*dw*win[iw]
    Bd2 += np.outer(np.conj(Mm['d2']['u']), Mm['d2']['r'])*dw*win[iw]
Cd2_cod = np.real(pref*(Bd2 + np.conj(Bd2.T))); Cd2_lit = np.real(pref*2*Bd2)
ref = np.load(os.path.join(OUT, 'correlation_decomp.npz')) if os.path.exists(os.path.join(OUT, 'correlation_decomp.npz')) else None
if ref is not None and ref['Cd2'].shape == Cd2_cod.shape:
    print(f"   reproduccion de out/correlation_decomp.npz (Cd2): dif rel max = "
          f"{np.abs(Cd2_cod-ref['Cd2']).max()/np.abs(ref['Cd2']).max():.1e}")


def ridge(Cm, thr, xmin=-12, xmax=-1):
    mu = (xc < xmax) & (xc > xmin); md = (xc > 0.5) & (xc < 14)
    sub = Cm[np.ix_(mu, md)]; xu = xc[mu]; xd = xc[md]; sl = []
    for i, x in enumerate(xu):
        j = np.argmin(sub[i])
        if sub[i, j] < thr: sl.append(xd[j]/x)
    return (np.median(sl) if sl else np.nan), len(sl)


rows = []
for name, Cm in (('d2, receta del codigo (como §3)', Cd2_cod), ('d2, lectura literal 2 Re B', Cd2_lit),
                 ('d2, codigo con amplitud corregida (x2)', 2*Cd2_cod)):
    for thr in (-0.5e-3, -1e-3, -2e-3):
        s, n = ridge(Cm, thr); rows.append((name, thr, s, n))
        print(f"   {name:42s} umbral {thr:+.1e}: pendiente mediana {s:+.2f}  ({n} puntos)")
# forma de la lengua: literal frente a codigo en el cuadrante cruzado
mu = (xc < -0.5) & (xc > -15); md = (xc > 0.5) & (xc < 15); Q = np.ix_(mu, md)
corr = np.corrcoef(Cd2_cod[Q].ravel(), Cd2_lit[Q].ravel())[0, 1]
print(f"   forma de la lengua (cuadrante cruzado |x|<15): correlacion codigo-literal = {corr:.4f}; "
      f"dif rel max = {np.abs(Cd2_cod-Cd2_lit)[Q].max()/np.abs(Cd2_cod[Q]).max():.2e}")
# C total con integracion completa (u,d1 hasta omega=8, d2 en la ventana): el canal d2 no cambia
Btot = total(8.0)
Ctot_cod = np.real(pref*(Btot + np.conj(Btot.T))); Ctot_lit = np.real(pref*2*Btot)
for name, Cm in (('TOTAL integracion completa, codigo', Ctot_cod), ('TOTAL integracion completa, literal', Ctot_lit)):
    s, n = ridge(Cm, -1e-3); print(f"   {name:42s} umbral -1.0e-03: pendiente mediana {s:+.2f}  ({n} puntos)")
print(f"fin {stamp()}")
