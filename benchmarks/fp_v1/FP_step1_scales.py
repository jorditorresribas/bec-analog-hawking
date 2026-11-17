"""
RUMBO #1 / Experimento B (Finazzi-Parentani) — PASO 1 (v2, auditado por ChatGPT):
fijar las escalas en NUESTRAS unidades y validarlas con DERIVACION + ASERCIONES.
Correcciones de la auditoria:
  (1) distinguir longitud de curacion EFECTIVA xi_lam = sqrt(lam)/(sqrt2 c) de la de
      REFERENCIA xi_ref = 1/(sqrt2 c) (= xi_lam en lam=1). xi_lam/ell=sqrt2 (indep lam);
      xi_ref/ell=sqrt(2/lam).
  (2) Lambda_FP es una ESCALA DE FRECUENCIA: de la Ec.9 de FP,
      Omega^2 = c^2 k^2 + c^4 k^4/Lambda_FP^2, igualando a c^2k^2+(lam/4)k^4 ->
      Lambda_FP = 2 c^2/sqrt(lam) = c * k_*.  NO es k_* (eso es un wavenumber).
      Forma FP: Omega^2 = c^2 k^2 (1 + c^2 k^2/Lambda_FP^2).
  (3) horizonte = raiz INTERPOLADA (c=v), no el punto de malla mas cercano.
  (4) extremos = bordes de la CAJA numerica (no x->+-inf matematico).
  (5) aserciones + residuos acumulados (reproducibilidad).
No se implementa kappa_bar ni barridos (paso 2).
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'solver'))  # layout del repo
import numpy as np
import correlations as C

g = C.g
SQRT2 = np.sqrt(2.0)

def scales(cval, lam):
    """Devuelve dict con las escalas en nuestras unidades."""
    kstar  = 2*cval/np.sqrt(lam)        # numero de onda de cruce dispersivo
    ell    = 1.0/kstar                  # = sqrt(lam)/(2c)
    xi_ref = 1.0/(SQRT2*cval)           # curacion de REFERENCIA (lam=1), fija
    xi_lam = np.sqrt(lam)/(SQRT2*cval)  # curacion EFECTIVA del modelo modificado
    Lam_FP = 2*cval**2/np.sqrt(lam)     # escala de FRECUENCIA de FP (= c*k_*)
    return dict(c=cval, kstar=kstar, ell=ell, xi_ref=xi_ref, xi_lam=xi_lam, Lam_FP=Lam_FP)

# ---- fondo (w=2) y HORIZONTE INTERPOLADO ----
x = np.linspace(-60, 60, 6000)
n, v, V = C.bg(x, w=2.0, mode='horizon')
c = np.sqrt(g*n)
f = c - v
# raiz interpolada de c-v=0
isg = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
t = f[isg]/(f[isg]-f[isg+1])                     # fraccion lineal
xH = x[isg] + t*(x[isg+1]-x[isg])
nH = n[isg] + t*(n[isg+1]-n[isg])
vH = v[isg] + t*(v[isg+1]-v[isg])
cH = np.sqrt(g*nH)
print("="*78)
print("PASO 1 v2 (exp. B): escalas dispersivas + validacion con aserciones")
print("="*78)
print(f"Horizonte INTERPOLADO: x_H={xH:.5f}  c_H={cH:.5f}  v_H={vH:.5f}  "
      f"|c-v|={abs(cH-vH):.2e}")
assert abs(cH - vH) < 1e-4, "horizonte interpolado no cumple c=v"

# ---- tabla de escalas (lam=1) en borde-izq, horizonte, borde-der ----
rows = [("borde izq. caja (x=-60)", c[0], n[0], v[0]),
        ("horizonte interp. (c=v)", cH, nH, vH),
        ("borde der. caja (x=+60)", c[-1], n[-1], v[-1])]
print(f"\nEscalas (lambda=1):")
print(f"{'punto':>26} | {'c':>7} {'v':>7} | {'xi_ref':>7} {'xi_lam':>7} {'ell':>7} "
      f"{'k_*':>7} {'Lam_FP':>7}")
for name, cv, nv, vv in rows:
    s = scales(cv, 1.0)
    print(f"{name:>26} | {cv:>7.4f} {vv:>7.4f} | {s['xi_ref']:>7.4f} {s['xi_lam']:>7.4f} "
          f"{s['ell']:>7.4f} {s['kstar']:>7.4f} {s['Lam_FP']:>7.4f}")

# ---- ASERCIONES de las relaciones derivadas (varios lam) ----
print("\nASERCIONES de las relaciones (lam=0.5,1,2) en c=c_H:")
for lam in [0.5, 1.0, 2.0]:
    s = scales(cH, lam)
    r_lam = s['xi_lam']/s['ell']           # debe ser sqrt2
    r_ref = s['xi_ref']/s['ell']           # debe ser sqrt(2/lam)
    lam_ckstar = s['Lam_FP']/(cH*s['kstar'])   # debe ser 1
    # FP-form reproduce la dispersion?
    kk = 0.7
    fp  = cH**2*kk**2*(1 + cH**2*kk**2/s['Lam_FP']**2)
    our = cH**2*kk**2 + lam*kk**4/4
    # xi_FP = sqrt2 c0/Lam_FP con c0=c local debe ser xi_lam (NO xi_ref salvo lam=1)
    xi_FP = SQRT2*cH/s['Lam_FP']
    assert abs(xi_FP - s['xi_lam']) < 1e-12
    assert abs(r_lam - SQRT2) < 1e-12
    assert abs(r_ref - np.sqrt(2/lam)) < 1e-12
    assert abs(lam_ckstar - 1) < 1e-12
    assert abs(fp - our) < 1e-12
    print(f"  lam={lam}: xi_lam/ell={r_lam:.6f}(=sqrt2) xi_ref/ell={r_ref:.6f}"
          f"(=sqrt(2/lam)={np.sqrt(2/lam):.6f}) Lam_FP/(c k_*)={lam_ckstar:.3e} "
          f"xi_FP=xi_lam={s['xi_lam']:.4f}(1/(sqrt2 c)={s['xi_ref']:.4f}) "
          f"FP-form err={abs(fp-our):.1e}  OK")

# ---- VALIDACION dispersion: cuartico del solver vs analitica (AMBOS bordes + HORIZONTE) ----
print("\nVALIDACION dispersion  |Omega^2 - (c^2k^2+lam k^4/4)|  (bordes + horizonte interp.):")
pts = [('borde izq', c[0], v[0]), ('horizonte', cH, vH), ('borde der', c[-1], v[-1])]
for lam in [0.5, 1.0, 2.0]:
    maxerr = 0.0; nprop = 0
    for _, cv, vv in pts:
        for omega in [0.05, 0.1, 0.2]:
            for k in C.quartic_roots_lam(omega, vv, cv, lam):
                if abs(k.imag) < 1e-7:
                    kr = k.real; Om = omega - vv*kr
                    maxerr = max(maxerr, abs(Om**2 - (cv**2*kr**2 + lam*kr**4/4)))
                    nprop += 1
    assert maxerr < 1e-10, f"dispersion no valida lam={lam}: {maxerr}"
    print(f"  lam={lam}: max err = {maxerr:.2e}  ({nprop} raices propagantes)  OK")

# ---- VALIDACION 2: simbolo del operador build_op_lambda, rama fisica ----
print("\nVALIDACION 2: simbolo build_op_lambda, autovalor fisico ~ omega (bordes + horizonte):")
for lam in [0.5, 1.0, 2.0]:
    maxerr = 0.0
    for _, cv, vv in pts:
        nv = cv**2/g
        for omega in [0.1, 0.2]:
            for k in C.quartic_roots_lam(omega, vv, cv, lam):
                if abs(k.imag) < 1e-7:
                    M = C._symbol_lam(k.real, nv, vv, lam, omega)
                    maxerr = max(maxerr, np.min(np.abs(np.linalg.eigvals(M) - omega)))
    assert maxerr < 1e-10, f"simbolo no valida lam={lam}: {maxerr}"
    print(f"  lam={lam}: max|autovalor_fisico - omega| = {maxerr:.2e}  OK")

print("\n" + "="*78)
print("MAPEO A FINAZZI-PARENTANI (CORREGIDO)")
print("="*78)
print("DETERMINADO sin ambiguedad por nuestra dispersion:")
print("  k_*(x)   = 2 c(x)/sqrt(lam)        [NUMERO DE ONDA de cruce]")
print("  Lam_FP(x)= 2 c(x)^2/sqrt(lam) = c k_*  [ESCALA DE FRECUENCIA de FP, Ec.9]")
print("  xi_lam(x)= sqrt(lam)/(sqrt2 c)    [curacion EFECTIVA; xi_lam/ell=sqrt2]")
print("  xi_ref(x)= 1/(sqrt2 c)            [curacion de REFERENCIA lam=1; =xi_lam en lam=1]")
print(f"  horizonte lam=1: k_*={2*cH:.4f}  Lam_FP={2*cH**2:.4f}  xi={1/(SQRT2*cH):.4f}")
print("\nFP: xi_FP = sqrt2 c0/Lam_FP. Con c0=c local -> xi_FP = sqrt(lam)/(sqrt2 c) = xi_lam")
print("    (curacion EFECTIVA; = 1/(sqrt2 c) SOLO en lam=1, no para lam!=1).")
print("CONVENCION / PENDIENTE (no fijar sin justificar):")
print("  - c0 de referencia de FP = velocidad de sonido de referencia del PERFIL FP,")
print("    no automaticamente la local en cada punto -> elegir y justificar.")
print("  - d~1.58 (Ec.24): FIT a los perfiles de FP, NO universal. Transferible SOLO el")
print("    escalado d_xi ~ c0/(kappa_K Lam_FP^2)^(1/3), con Lam_FP de FRECUENCIA (no k_*).")
print("  - reparto d_xi^L vs d_xi^R y centro x_hor (Ec.20): dependen del fondo, pendientes.")
print("\nTODAS LAS ASERCIONES PASARON.")
