"""
DIAGNOSTICO (no gate, sin commit) — normalizacion de la correlacion densidad-densidad de correlations.py
(2026-11-17). No modifica el solver: importa correlations.py y replica literalmente su corchete.

(A) Simetria del corchete en un caso de horizonte (lam = 1, tanh w = 2): el paper (Recati et al., arXiv:0907.4305v1,
    Ec. 25 segun el cotejo literal del revisor) suma B + c.c. (conjugado complejo en el mismo (x, x')), es decir 2 Re B.
    El codigo suma B + conj(B.T) y toma la parte real: Re B(x,x') + Re B(x',x). Son iguales si Re B es simetrica.
(B) Test ABSOLUTO de normalizacion en un condensado homogeneo subsonico (Mach 0,5), donde la correlacion de vacio es
    exacta: <dn(x) dn(x')>_conexa sin termino propio = n * int dk/2pi [S(k) - 1] e^{ik(x-x')},
    S(k) = eps_k / E_k (factor de estructura estatico de Bogoliubov; invariante galileana a tiempo igual).
    Se construye C con EXACTAMENTE la receta de density_correlation (modos |J| = 1, canales u y d1, corchete
    v* r + c.c. simetrizado, prefactor sqrt(n n')/(4 pi)) y se compara con la formula exacta restringida al mismo
    rango de k que cubren las omegas integradas. Si la amplitud difiere en un factor constante, ese factor es
    el error de normalizacion de la receta (independiente de como se normalicen los modos).
Salida: out/diag_recati_norm.csv
"""
import os, sys, io, csv, time, contextlib
import numpy as np

HERE = os.path.dirname(os.path.abspath(__file__)); ROOT = os.path.dirname(HERE); OUT = os.path.join(ROOT, 'out')
sys.path.insert(0, os.path.join(ROOT, 'solver'))
with contextlib.redirect_stdout(io.StringIO()):
    import correlations as C
t0 = time.time(); rows = []

# ---------------- (A) simetria del corchete, horizonte lam = 1 ----------------
L, N, NC = 40.0, 9600, 320
idx = np.linspace(0, N-1, NC).astype(int)
print("(A) simetria del corchete B (antes de simetrizar), horizonte tanh w=2, lam=1, L=40, N=9600:")
for om in (0.02, 0.05, 0.10):
    md = {}
    for ch in ('u', 'd1', 'd2'):
        r = C.scatter_field(L, N, om, ch, 2.0, 1.0, 'horizon')
        u = r['u'][idx]; vv = r['v_mode'][idx]; md[ch] = dict(u=u, v=vv, r=u+vv)
    B = sum(np.outer(np.conj(md[I]['v']), md[I]['r']) for I in ('u', 'd1')) + np.outer(np.conj(md['d2']['u']), md['d2']['r'])
    nB = np.abs(B).max()
    asym_re = np.abs(B.real - B.real.T).max()/nB
    asym_full = np.abs(B - B.T).max()/nB
    code = np.real(B + np.conj(B.T)); paper = 2*np.real(B)
    dif = np.abs(code - paper).max()/np.abs(paper).max()
    print(f"   omega={om}: max|B - B^T|/max|B| = {asym_full:.2e}   max|ReB - ReB^T|/max|B| = {asym_re:.2e}"
          f"   max|codigo - 2ReB|/max|2ReB| = {dif:.2e}")
    rows.append(dict(prueba='A_simetria', omega=om, asim_B=asym_full, asim_ReB=asym_re, dif_codigo_paper=dif))

# ---------------- (B) homogeneo subsonico: test absoluto ----------------
g = C.g; n0 = C.n_u; v0 = C.vLb; c0 = np.sqrt(g*n0)
L, N, NC = 40.0, 8000, 241
idx = np.linspace(0, N-1, NC).astype(int)
x = np.linspace(-L, L, N)[idx]
OM = np.linspace(0.004, 1.2, 300)       # trapecio en omega
acc = np.zeros((NC, NC), complex); kpos = []; kneg = []
for iw, om in enumerate(OM):
    B = np.zeros((NC, NC), complex)
    for ch in ('u', 'd1'):
        r = C.scatter_field(L, N, om, ch, 2.0, 1.0, 'homog')
        if ch == 'u':
            kpos.append(r['Lin']['k'])
        else:
            kneg.append(r['Rin_pos']['k'])
        u = r['u'][idx]; vv = r['v_mode'][idx]
        B += np.outer(np.conj(vv), u+vv)
    B = B + np.conj(B.T)
    wgt = (OM[1]-OM[0])*(0.5 if iw in (0, len(OM)-1) else 1.0)
    acc += B*wgt
Ccode = np.real(np.sqrt(n0*n0)/(4*np.pi)*acc)
kpos = np.array(kpos); kneg = np.array(kneg)
print(f"\n(B) homogeneo Mach {v0/c0:.2f}: k cubierto (u) [{kpos.min():.4f}, {kpos.max():.4f}], (d1) [{kneg.max():.4f}, {kneg.min():.4f}]")


def exact(dx_):
    """n * int dk/2pi (S-1) cos(k dx) sobre los mismos intervalos de k que cubren las omegas integradas."""
    out = 0.0
    for a, b in ((kpos.min(), kpos.max()), (kneg.min(), kneg.max())):
        kk = np.linspace(a, b, 20001); eps = kk**2/2; E = np.sqrt(eps*(eps + 2*g*n0))
        out += np.trapz((eps/E - 1)*np.cos(kk*dx_), kk)/(2*np.pi)
    return n0*out


ic = NC//2
seps = np.arange(0, 40)
cc = np.array([Ccode[ic, ic+s] for s in seps]); dxs = x[ic+seps] - x[ic]
ee = np.array([exact(d) for d in dxs])
ratio = np.sum(cc*ee)/np.sum(ee*ee)                     # factor de mejor ajuste codigo/exacto
resid = np.max(np.abs(cc - ratio*ee))/np.max(np.abs(ee))
print(f"   fila central, {len(seps)} separaciones 0..{dxs[-1]:.2f}:")
for s in (0, 1, 2, 4, 8, 16, 32):
    print(f"     dx={dxs[s]:6.3f}: codigo {cc[s]:+.6e}   exacto {ee[s]:+.6e}   cociente {cc[s]/ee[s]:+.5f}")
print(f"   FACTOR codigo/exacto (minimos cuadrados) = {ratio:.6f}   residuo tras el factor = {resid:.1e}")
print(f"   referencias: 1/2 = 0.5, 1/(2 pi) = {1/(2*np.pi):.5f}, 2 = 2, 4 pi = {4*np.pi:.4f}")
for s in range(len(seps)):
    rows.append(dict(prueba='B_homogeneo', omega='', dx=dxs[s], C_codigo=cc[s], C_exacto=ee[s]))
rows.append(dict(prueba='B_factor', omega='', factor=ratio, residuo=resid))
keys = ['prueba', 'omega', 'asim_B', 'asim_ReB', 'dif_codigo_paper', 'dx', 'C_codigo', 'C_exacto', 'factor', 'residuo']
with open(os.path.join(OUT, 'diag_recati_norm.csv'), 'w', newline='') as f:
    w = csv.DictWriter(f, fieldnames=keys); w.writeheader()
    for r in rows: w.writerow({k: r.get(k, '') for k in keys})
print(f"fin ({time.time()-t0:.0f} s)")
