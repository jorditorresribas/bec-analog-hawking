"""
AUDITORIA de PADE_SPEC.md (4 pruebas de ChatGPT) ANTES de escribir el solver.
TOLERANCIA FIJADA: error de dispersion Omega_P^2 < 1% en la banda fisica (para no
contaminar observables al nivel del % que se busca).
Weyl (gate0) se usa como referencia. Operadores en DENSO (auditoria, N moderado).
"""
import numpy as np, contextlib, os
src = open('residual_anatomy/gate0_pseudoherm.py').read().split('# ===========')[0]
ns = {}; exec(src, ns)
weyl_kinetic = ns['weyl_kinetic']; bg = ns['bg']; nval = ns['nval']; g = ns['g']; mu = ns['mu']

G0 = g*1.0  # homogeneo n=1
def eps_exact(k, lam): return -G0 + np.sqrt(G0**2 + G0*k**2 + lam*k**4/4)
def eps_pade(k, lam):  s = k**2; return s/2 + (lam-1)*s**2/(4*(2*G0+s))

print("="*70); print("AUDITORIA PADE (4 pruebas de ChatGPT). Tol dispersion Om^2 < 1% en banda"); print("="*70)

# --- PRUEBA 1: error de eps_P y Omega_P^2 en TODA la banda, cada lam ---
print("\n[1] Error relativo en la banda fisica k in (0, 2.5], malla densa")
kb = np.linspace(0.05, 2.5, 200)
print(f"    {'lam':>5} {'max err eps_P':>13} {'max err Om^2':>13} {'k donde Om^2=1%':>16}")
for lam in (0.5, 1.0, 2.0, 4.0):
    ee = eps_exact(kb, lam); ep = eps_pade(kb, lam)
    Om2e = ee*(ee+2*G0); Om2p = ep*(ep+2*G0)
    erre = np.abs(ep-ee)/np.abs(ee); errO = np.abs(Om2p-Om2e)/np.abs(Om2e)
    over = kb[errO > 0.01]
    kcut = f"{over[0]:.2f}" if len(over) else ">2.5 (ok banda)"
    print(f"    {lam:5.1f} {erre.max():13.2e} {errO.max():13.2e} {kcut:>16}")
print("    -> Om^2<1% se cumple hasta el k indicado; por encima, [2/1] insuficiente (subir orden).")

# --- operador Pade cinetico inhomogeneo denso:  D/2 + (lam-1)/4 D (2G(x)+D)^-1 D ---
def Dmat(x):
    N = len(x); dx = x[1]-x[0]; e = np.ones(N-1)
    D2 = (np.diag(e, 1) - 2*np.diag(np.ones(N)) + np.diag(e, -1))/dx**2
    return -D2   # D = -d^2/dx^2, autovalores ~ +k^2 >=0
def Kpade(x, lam, w=2.0, mode='horizon'):
    n, v, V = bg(x, w, mode); G = g*n
    D = Dmat(x); A = np.linalg.inv(np.diag(2*G) + D)       # (2G(x)+D)^-1
    return 0.5*D + (lam-1)/4.0 * (D @ A @ D)

# --- PRUEBA 2: positividad/invertibilidad de 2G(x)+D (homog y horizonte) ---
print("\n[2] Positividad e invertibilidad de 2G(x)+D (min autovalor > 0)")
for mode in ('homog', 'horizon'):
    x = np.linspace(-40, 40, 400); n, v, V = bg(x, 2.0, mode); G = g*n
    ev = np.linalg.eigvalsh(np.diag(2*G) + Dmat(x))
    print(f"    {mode:>8}: min autovalor(2G+D) = {ev.min():.4f}  (2G_min={2*G.min():.3f}) -> "
          f"{'POSITIVO DEFINIDO, invertible' if ev.min() > 0 else 'NO POSITIVO (revisar)'}")

# --- PRUEBA 3: espectro del cinetico Pade: ramas/autovalores espurios ---
print("\n[3] Espectro de Kpade (cinetico): debe ser real>=0, reproducir eps_P(k), sin espurios")
for mode in ('homog', 'horizon'):
    x = np.linspace(-40, 40, 400)
    Kp = Kpade(x, 2.0, mode=mode)
    ev = np.linalg.eigvals(Kp)
    im = np.abs(ev.imag).max(); neg = ev.real.min()
    print(f"    {mode:>8} lam=2: max|Im(ev)|={im:.2e}  min Re(ev)={neg:.3f}  "
          f"max Re(ev)={ev.real.max():.1f}  (eps_P es real>=0; sin ramas complejas espurias)")

# --- PRUEBA 4: Pade inhomogeneo vs Weyl inhomogeneo (accion), horizonte vs asintotico ---
print("\n[4] Kpade vs Kweyl en fondo con HORIZONTE (lam): accion sobre funciones de prueba")
print("    (lam=1 mide la diferencia FD-vs-espectral = linea base de discretizacion)")
x = np.linspace(-60, 60, 600)
xH = 0.0  # horizonte ~ x=0.46 en el perfil; region central
ihor = np.abs(x) < 6.0          # region del horizonte
iasy = np.abs(x) > 30.0         # region asintotica
for lam in (1.0, 2.0):
    Kp = Kpade(x, lam, mode='horizon')
    Kw = weyl_kinetic(x, lam, mode='horizon')
    # funcion de prueba suave y paquete localizado cerca del horizonte
    fs = np.exp(-(x/15.0)**2)*np.cos(0.5*x)
    fl = np.exp(-(x/3.0)**2)*np.cos(1.0*x)          # localizado cerca del horizonte
    for name, f in (("suave", fs), ("loc.horizonte", fl)):
        d = (Kp@f) - (Kw@f)
        rel_h = np.linalg.norm(d[ihor])/max(np.linalg.norm((Kw@f)[ihor]), 1e-30)
        rel_a = np.linalg.norm(d[iasy])/max(np.linalg.norm((Kw@f)[iasy]), 1e-30)
        print(f"    lam={lam} {name:>13}: |Kp-Kw| rel  horizonte={rel_h:.2e}  asintotico={rel_a:.2e}")
print("    (lam=1: dif = discretizacion FD/espectral; lam=2: + aprox Pade + ordenamiento inhom)")
