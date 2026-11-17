"""
RUMBO #1 / Opcion A — comprobacion pedida por ChatGPT: la ASIMETRIA entre ramas
(lambda<1 vs lambda>1) con w, es real o la introduce normalizar cada config a su
propio lambda=1? Representamos R_lambda(w) normalizado Y los str SIN normalizar, y
las desviaciones ABSOLUTAS Delta respecto a lambda=1.

Datos de w_dependence_v3.py (log; deterministas). |str(1)| CRECE con w (banda termica
captura mas / ventana espacial ~w): el denominador de R cambia con w.
"""
import numpy as np
import matplotlib
matplotlib.use('Agg')
import matplotlib.pyplot as plt

ws = np.array([2.0, 4.0, 6.0])
# str(w, lam) sin normalizar (del log v3)
STR = {
    0.5: np.array([-0.2316, -0.2716, -0.3044]),
    1.0: np.array([-0.1890, -0.2279, -0.2656]),
    2.0: np.array([-0.1698, -0.1943, -0.2263]),
}
R = {l: STR[l]/STR[1.0] for l in (0.5, 1.0, 2.0)}
# desviaciones ABSOLUTAS respecto a lambda=1
D_lo = STR[0.5] - STR[1.0]     # rama hidrodinamica (lam<1)
D_hi = STR[2.0] - STR[1.0]     # rama dispersiva  (lam>1)

print("Comprobacion de la asimetria (absoluto vs normalizado):")
print(f"{'w':>4} | {'|str(1)|':>9} | {'R0.5':>6} {'R2':>6} | "
      f"{'D_lo(abs)':>9} {'D_hi(abs)':>9} | {'R0.5-1':>7} {'R2-1':>7}")
for i, w in enumerate(ws):
    print(f"{w:>4} | {abs(STR[1.0][i]):>9.4f} | {R[0.5][i]:>6.3f} {R[2.0][i]:>6.3f} | "
          f"{D_lo[i]:>+9.4f} {D_hi[i]:>+9.4f} | {R[0.5][i]-1:>+7.3f} {R[2.0][i]-1:>+7.3f}")

print("\nLECTURA:")
print(f"  |str(1)| CRECE con w: {abs(STR[1.0][0]):.4f} -> {abs(STR[1.0][2]):.4f} "
      f"(x{abs(STR[1.0][2])/abs(STR[1.0][0]):.2f})")
print(f"  D_lo (lam<1, ABS) ~ CONSTANTE: {D_lo[0]:+.4f} -> {D_lo[2]:+.4f}")
print(f"  D_hi (lam>1, ABS) CRECE:       {D_hi[0]:+.4f} -> {D_hi[2]:+.4f} "
      f"(x{D_hi[2]/D_hi[0]:.2f})")
print("  => El 'R0.5 -> 1 al crecer w' es casi todo DENOMINADOR creciente (D_lo plano),")
print("     NO una mejora real de robustez en la rama hidrodinamica.")
print("  => En ABSOLUTO, la rama dispersiva (lam>1) se DESVIA MAS con w (D_hi crece).")
print("  => La asimetria de R esta dominada por la normalizacion (confirma sospecha ChatGPT).")

fig, ax = plt.subplots(1, 3, figsize=(16, 4.8), constrained_layout=True)
# (1) str sin normalizar
for l, c, m in [(0.5, 'tab:blue', 'o'), (1.0, 'k', 's'), (2.0, 'tab:red', '^')]:
    ax[0].plot(ws, -STR[l], m+'-', color=c, label=f"lam={l}")
ax[0].set_xlabel("w"); ax[0].set_ylabel("|str| (sin normalizar)")
ax[0].set_title("Fuerza de la lengua SIN normalizar\n(|str(1)| crece con w)")
ax[0].legend(); ax[0].grid(alpha=0.3)
# (2) R normalizado
ax[1].plot(ws, R[0.5], 'o-', color='tab:blue', label='R(0.5)')
ax[1].plot(ws, R[2.0], '^-', color='tab:red', label='R(2)')
ax[1].axhline(1, color='k', lw=0.5, alpha=0.4)
ax[1].set_xlabel("w"); ax[1].set_ylabel("R = str/str(1)")
ax[1].set_title("Normalizado a lam=1\n(R0.5->1 por denominador, no por fisica)")
ax[1].legend(); ax[1].grid(alpha=0.3)
# (3) desviaciones absolutas
ax[2].plot(ws, D_lo, 'o-', color='tab:blue', label='D_lo = str(0.5)-str(1)')
ax[2].plot(ws, D_hi, '^-', color='tab:red', label='D_hi = str(2)-str(1)')
ax[2].axhline(0, color='k', lw=0.5, alpha=0.4)
ax[2].set_xlabel("w"); ax[2].set_ylabel("desviacion ABSOLUTA vs lam=1")
ax[2].set_title("Desviaciones absolutas\n(lam<1 plano, lam>1 crece)")
ax[2].legend(); ax[2].grid(alpha=0.3)
fig.suptitle("Asimetria de ramas: absoluto vs normalizado (comprobacion ChatGPT)", fontsize=12)
fig.savefig("out/asymmetry_check.png", dpi=120)
print("\nFig -> out/asymmetry_check.png")
