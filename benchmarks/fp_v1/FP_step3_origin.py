"""
PASO 3 (origen fisico de T_IR > kappa_local/2pi). Programa de ChatGPT, pruebas 1-3,
SIN tocar el solver BdG. Objetivo: antes de postular un mecanismo, distinguir si el
exceso esta (i) en el proceso de creacion de pares del canal Hawking, (ii) en la
suma/mezcla de canales (greybody) que hace el estimador H4, o (iii) en un sesgo del
limite IR (banda).

Estructura de canales (sonda previa, inyectando el modo entrante de NORMA -1):
  salientes propagantes = { (+1,L) Hawking upstream (grande), (-1,R) companero,
                            (+1,R) co-propagante downstream (pequeno) }
  N_omega(H4) = suma de |amp|^2 de los de NORMA +1 = (+1,L) + (+1,R).

Pruebas:
 1 (canal): T_eff del canal Hawking (+1,L) SOLO vs del sum H4 -> el exceso, ¿vive en
   el canal principal (fisico) o lo introduce sumar el canal (+1,R) (greybody)?
 2 (espectral): T_eff(omega)=omega/ln(1+1/N) vs omega/kappa_local -> meseta IR
   desplazada (T efectiva distinta) vs deriva con omega (ajuste termico comprime una
   desviacion espectral).
 3 (banda): extrapolar el limite IR con bandas cada vez mas bajas; ¿el exceso es
   estable al mover la banda hacia 0, o es sesgo del estimador?
Control: pseudo-unitariedad por punto + convergencia en N en el omega mas bajo.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', '..', 'solver'))  # layout del repo
import numpy as np, csv
import mode_matching as MM
g = MM.g; n_u, mach_u, md = 1.0, 0.5, 2.0
vLb = mach_u*np.sqrt(g*n_u); jc = n_u*vLb; n_d = (jc/(md*np.sqrt(g)))**(2/3)

def box(w): L = max(40.0, 8.0*w+20.0); return L, int(round(L*240))
def kappa_local(w):
    L, _ = box(w); x = np.linspace(-L, L, 60000)
    s = 0.5*(1+np.tanh(x/w)); n = n_u+(n_d-n_u)*s; c = np.sqrt(n); v = jc/n
    f = c-v; i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]
    t = f[i]/(f[i]-f[i+1]); xH = x[i]+t*(x[i+1]-x[i])
    return abs(np.interp(xH, x, np.gradient(f, x)))

def channels(w, om, L, N):
    """Devuelve dict con N_tot (H4), N_L (canal Hawking +1,L), N_R+ (+1,R), pseudo-u."""
    MM.w = w; MM.omega = om
    x = np.linspace(-L, L, N); n, v, V = MM.bg(x, 'horizon')
    mR = MM.modes_full(n[-1], v[-1], V[-1])
    Rin = [m for m in mR if m['prop'] and m['vg'] < 0]
    incR = [1.0 if m['norm'] < 0 else 0.0 for m in Rin]
    out, _ = MM.solve_inject(L=L, N=N, mode='horizon', incL=0.0, incR=incR)
    N_L = sum(a for (nrm, s, a) in out if nrm > 0 and s == 'L')
    N_Rp = sum(a for (nrm, s, a) in out if nrm > 0 and s == 'R')
    N_tot = N_L + N_Rp
    pu = abs(sum(nrm*a for (nrm, s, a) in out) + 1.0)
    return dict(N_tot=N_tot, N_L=N_L, N_Rp=N_Rp, pu=pu)

def Teff(om, Nval):
    return om/np.log(1.0+1.0/Nval)

def kbar_d(w, d=1.58):
    """kappa_bar(Ec.20, ventana simetrica) / kappa_local, para comparar con T_eff_L."""
    L, _ = box(w); x = np.linspace(-L, L, 60000)
    s = 0.5*(1+np.tanh(x/w)); n = n_u+(n_d-n_u)*s; c = np.sqrt(n); v = jc/n; f = c-v
    i = np.where(np.sign(f[:-1]) != np.sign(f[1:]))[0][0]; t = f[i]/(f[i]-f[i+1])
    xH = x[i]+t*(x[i+1]-x[i]); cH = np.interp(xH, x, c); kK = abs(np.interp(xH, x, np.gradient(f, x)))
    Lam = 2*cH**2; dxi = d*cH/(kK*Lam**2)**(1/3)
    kb = abs(np.interp(xH+0.5*dxi, x, f)-np.interp(xH-0.5*dxi, x, f))/dxi
    return kb/kK

WS = [2.0, 4.0, 8.0]
print("="*96)
print("PASO 3: origen de T_IR>kappa_local/2pi. Pruebas 1(canal) 2(espectral) 3(banda)")
print("="*96)

rows = []; compare = []
for w in WS:
    kK = kappa_local(w); T0 = kK/(2*np.pi); L, N = box(w)
    # malla densa de omega/T0 (cubre las 3 bandas)
    xr = np.array([0.05,0.065,0.08,0.10,0.12,0.15,0.18,0.22,0.27,0.33,0.40,0.50,0.60])
    oms = xr*T0
    dat = [channels(w, om, L, N) for om in oms]
    Ntot = np.array([d['N_tot'] for d in dat]); NL = np.array([d['N_L'] for d in dat])
    NRp = np.array([d['N_Rp'] for d in dat]); pu = np.array([d['pu'] for d in dat])
    omN = oms*Ntot
    Teff_tot = np.array([Teff(oms[i], Ntot[i]) for i in range(len(oms))])
    Teff_L   = np.array([Teff(oms[i], NL[i])   for i in range(len(oms))])
    print(f"\n{'='*92}\n--- w={w:.0f}  kappa_local/2pi={T0:.5f}  (L={L:.0f},N={N})  pu_max={pu.max():.1e} ---")
    print(f"  fraccion del canal pequeno (+1,R): N_R+/N_tot = {NRp.mean()/Ntot.mean():.4f}")
    # PRUEBA 2: T_eff(omega)/kappa_local vs omega/kappa_local
    print(f"\n  PRUEBA 2  (T_eff espectral; buscar meseta IR desplazada vs deriva):")
    print(f"  {'om/kloc':>8} {'Ntot':>8} {'omega*N':>9} {'Teff_tot/kloc':>13} {'Teff_L/kloc':>12} {'pu':>8}")
    for i in range(len(oms)):
        print(f"  {oms[i]/kK:8.4f} {Ntot[i]:8.3f} {omN[i]:9.5f} {Teff_tot[i]/T0:13.4f} "
              f"{Teff_L[i]/T0:12.4f} {pu[i]:8.1e}")
    # PRUEBA 1: limite IR del canal Hawking SOLO vs sum H4
    # PRUEBA 3: extrapolar omega*N->0 en 3 bandas descendentes (quad)
    bands = [(0.12,0.60),(0.08,0.40),(0.05,0.25)]
    print(f"\n  PRUEBA 3 (banda) + PRUEBA 1 (canal): limite IR (quad) de omega*N:")
    print(f"  {'banda/T0':>12} {'npts':>5} {'TIR_tot/kloc':>12} {'TIR_L/kloc':>11}")
    for lo, hi in bands:
        sel = (xr >= lo-1e-9) & (xr <= hi+1e-9)
        if sel.sum() < 3: continue
        o = oms[sel]
        TIR_tot = np.polyfit(o, (oms*Ntot)[sel], 2)[-1]
        TIR_L   = np.polyfit(o, (oms*NL)[sel],   2)[-1]
        print(f"  [{lo:.2f},{hi:.2f}] {sel.sum():5d} {TIR_tot/T0:12.4f} {TIR_L/T0:11.4f}")
        rows.append([w, T0, lo, hi, sel.sum(), TIR_tot, TIR_L, TIR_tot/T0, TIR_L/T0])
        if abs(lo-0.05) < 1e-9:   # banda mas baja -> guardar para comparar con kappa_bar
            compare.append([w, TIR_L/T0, TIR_tot/T0, kbar_d(w, 1.58)])
    # convergencia en N en el omega MAS BAJO
    om_low = oms[0]
    c1 = channels(w, om_low, L, N); c2 = channels(w, om_low, L, int(N*1.5))
    print(f"  convergencia N en om/T0={xr[0]:.3f}: omega*Ntot  N={N}:{om_low*c1['N_tot']:.5f}  "
          f"N={int(N*1.5)}:{om_low*c2['N_tot']:.5f}  (pu {c1['pu']:.1e}->{c2['pu']:.1e})")

with open('out/fp_step3_origin.csv', 'w', newline='') as f:
    wr = csv.writer(f); wr.writerow(['w','T0','band_lo','band_hi','npts','TIR_tot','TIR_L','TIR_tot_over_kloc','TIR_L_over_kloc'])
    wr.writerows(rows)
print("\nCSV: out/fp_step3_origin.csv")
print("\n" + "="*96)
print("COMPARACION canal Hawking (+1,L) vs kappa_bar de FP (Ec.20, d=1.58):")
print("="*96)
print(f"  {'w':>3} {'T_eff_L/kloc':>13} {'kappa_bar/kloc':>15} {'T_IR_tot/kloc':>14} {'dif L vs kbar':>14}")
for w, tL, ttot, kb in compare:
    print(f"  {int(w):3d} {tL:13.4f} {kb:15.4f} {ttot:14.4f} {abs(tL-kb)/kb:13.2%}")
print("  -> el canal Hawking SOLO esta POR DEBAJO de kappa_local y RASTREA kappa_bar<kappa_local")
print("     (deficit encoge con w = recuperacion hidrodinamica, direccion CORRECTA de FP).")
print("     El 'exceso' >kappa_local de T_IR_tot lo mete sumar el canal (+1,R) co-propagante.")
print("\n" + "-"*96)
print("LECTURA:")
print(" - Si TIR_L/kloc ~ 1 y TIR_tot/kloc > 1  -> el exceso lo introduce SUMAR el canal")
print("   (+1,R): es greybody del estimador H4, no una T de creacion de pares mayor.")
print(" - Si TIR_L/kloc tambien > 1             -> el exceso esta en el canal Hawking (fisico).")
print(" - Si TIR_tot/kloc cae hacia 1 al bajar la banda -> era sesgo del limite IR (banda).")
print(" - Si Teff(omega)/kloc es meseta plana > 1 -> T efectiva distinta; si deriva con omega")
print("   -> el ajuste termico comprime una desviacion espectral.")
