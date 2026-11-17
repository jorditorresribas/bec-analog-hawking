"""
RUMBO #1 (b) — PROFUNDIZAR ROBUSTEZ, paso 1: CONVERGENCIA del propio sweep.

La robustez ~+-10% de la lengua de Hawking frente a lam se midio a (N=3000, L=30)
sin estudio de convergencia (caveat 2 de CORRELATIONS.md). Aqui re-medimos la
fuerza de la lengua (misma metrica integral que run_lambda_sweep) variando:
  (A) N   (resolucion)  a L fijo
  (B) L   (tamano de caja) a espaciado grueso fijo
para lam = {0.5, 1.0, 2.0}. Verificamos que el cociente fuerza/lam1 es estable.

Metrica (identica a run_lambda_sweep): integral de la parte NEGATIVA del canal d2
(Hawking puro) en el cuadrante cruzado upstream x downstream, |x|<14, con dx^2 del
grid grueso. Para comparar entre L mantenemos el espaciado grueso dc FIJO
(ncoarse = round(2L/dc)) -> la metrica es comparable.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import correlations as C

w = 2.0
lams = [0.5, 1.0, 2.0]
DC = 2*30.0/360          # espaciado grueso de referencia (run_lambda_sweep: L=30,ncoarse=360)

def omega_max(lam):
    lo, hi = 0.05, 1.3
    for _ in range(26):
        mid = 0.5*(lo+hi)
        r = C.scatter_field(L=24, N=2000, omega=mid, channel='d2', lam=lam,
                            mode='horizon', op='lambda')
        lo, hi = (mid, hi) if (r['Rin_neg'] is not None and r['nchan'] == 2) else (lo, mid)
    return lo

def tongue_strength(L, N, lam, nomega=90, ncoarse=None):
    """Fuerza integral de la lengua (parte negativa de C_d2 en cuadrante cruzado)."""
    if ncoarse is None:
        ncoarse = int(round(2*L/DC))
    omax = omega_max(lam)
    omegas = np.linspace(0.02, 0.9*omax, nomega)
    xc, Cd2, _ = C.density_correlation(L, N, omegas, w=w, lam=lam, mode='horizon',
                                       ncoarse=ncoarse, terms='d2', op='lambda')
    dx = xc[1]-xc[0]
    iu = (xc < -0.5) & (xc > -14); idn = (xc > 0.5) & (xc < 14)
    blk = Cd2[np.ix_(iu, idn)]
    return blk[blk < 0].sum()*dx**2, omax

def ratios(vals):
    return {l: vals[l]/vals[1.0] for l in lams}

print("="*70)
print("CONVERGENCIA DEL SWEEP DE ROBUSTEZ  (metrica integral, canal d2)")
print("="*70)

# --- (A) convergencia en N a L=30 fijo ---
print("\n(A) Convergencia en N  (L=30, ncoarse=360):")
print(f"{'N':>6} | " + " ".join(f"lam={l:<4} str" for l in lams) + " |  fuerza/lam1")
rowsN = {}
for N in [2000, 3000, 4500]:
    s = {}
    for l in lams:
        s[l], _ = tongue_strength(30.0, N, l, ncoarse=360)
    r = ratios(s); rowsN[N] = (s, r)
    print(f"{N:>6} | " + " ".join(f"{s[l]:+.3e}" for l in lams) +
          " |  " + " ".join(f"{r[l]:.3f}" for l in lams))

# --- (B) convergencia en L a espaciado grueso fijo ---
print("\n(B) Convergencia en L  (dx fino ~ 2L/N fijo, ncoarse=round(2L/dc)):")
print(f"{'L':>5} {'N':>6} | " + " ".join(f"lam={l:<4} str" for l in lams) + " |  fuerza/lam1")
rowsL = {}
dxfine = 2*30.0/3000     # mantener resolucion fina ~0.02
for L in [24.0, 30.0, 40.0]:
    N = int(round(2*L/dxfine))
    s = {}
    for l in lams:
        s[l], _ = tongue_strength(L, N, l)
    r = ratios(s); rowsL[L] = (s, r, N)
    print(f"{L:>5.0f} {N:>6} | " + " ".join(f"{s[l]:+.3e}" for l in lams) +
          " |  " + " ".join(f"{r[l]:.3f}" for l in lams))

print("\n" + "="*70)
print("LECTURA: el cociente fuerza/lam1 (robustez) debe ser estable en N y L.")
print("  Referencia (sesion previa, N=3000/L=30): 1.02 / 1.00 / 0.91")
rN = rowsN[4500][1]; rL = rowsL[40.0][1]
print(f"  Mejor N (4500,L=30): {rN[0.5]:.3f} / 1.000 / {rN[2.0]:.3f}")
print(f"  Mejor L (40,N={rowsL[40.0][2]}): {rL[0.5]:.3f} / 1.000 / {rL[2.0]:.3f}")
