"""
Validacion del solver GPE1D contra dos soluciones analiticas exactas.
Genera PNGs en out/ e imprime errores cuantitativos.
"""
import os as _os, sys as _sys; _sys.path.insert(0, _os.path.join(_os.path.dirname(_os.path.abspath(__file__)), '..', 'solver'))  # layout del repo
import numpy as np
import matplotlib
matplotlib.use("Agg")
import matplotlib.pyplot as plt
from gpe import GPE1D

OUT = "out"


def test_bright_soliton():
    """g=-1, sin trampa. Soliton exacto:
       psi(x,t) = eta*sech(eta*(x-v t)) * exp(i[v x - (v^2-eta^2) t/2])
       Debe trasladarse a velocidad v SIN cambiar de forma.
    """
    print("\n=== TEST 1: Soliton brillante (g=-1) ===")
    eta, v = 1.0, 2.0
    sim = GPE1D(L=80.0, N=4096, g=-1.0)
    x = sim.x
    psi0 = eta / np.cosh(eta * x) * np.exp(1j * v * x)

    T = 10.0
    dt = 1e-3
    nsteps = int(T / dt)
    psi = sim.evolve_real(psi0.astype(complex), dt, nsteps)

    # solucion analitica en t=T (centro desplazado a x=v*T, con periodicidad L)
    xc = (v * T + sim.L/2) % sim.L - sim.L/2
    dx_per = (x - xc + sim.L/2) % sim.L - sim.L/2
    psi_exact = eta / np.cosh(eta * dx_per) * np.exp(
        1j * (v * dx_per - (v**2 - eta**2) * T / 2))

    dens = np.abs(psi)**2
    dens_exact = np.abs(psi_exact)**2
    err = np.sqrt(np.sum((dens - dens_exact)**2) * sim.dx /
                  (np.sum(dens_exact**2) * sim.dx))
    norm_err = abs(sim.norm(psi) - sim.norm(psi0)) / sim.norm(psi0)
    peak = dens.max()
    print(f"  Norma conservada:      error = {norm_err:.2e}  (deberia ~0)")
    print(f"  Pico densidad:         {peak:.5f}  (exacto = {eta**2:.5f})")
    print(f"  Error RMS de forma:    {err:.2e}  (deberia ser pequeno)")

    plt.figure(figsize=(9, 4))
    plt.plot(x, dens_exact, 'k-', lw=3, alpha=0.4, label='analitico (t=10)')
    plt.plot(x, dens, 'r--', lw=1.2, label='numerico (t=10)')
    plt.plot(x, np.abs(psi0)**2, 'b:', lw=1, label='inicial (t=0)')
    plt.xlabel('x'); plt.ylabel(r'$|\psi|^2$'); plt.legend()
    plt.title(f'Soliton brillante: traslacion v={v}  (err RMS={err:.1e})')
    plt.xlim(-40, 40); plt.tight_layout()
    plt.savefig(f"{OUT}/test1_soliton.png", dpi=110)
    plt.close()
    return err < 1e-2


def test_thomas_fermi():
    """g>0, trampa armonica V=1/2 w^2 x^2. Estado fundamental via tiempo
       imaginario vs perfil Thomas-Fermi:
         n_TF(x) = max(mu - 1/2 w^2 x^2, 0)/g,  R=(3Ng/2w^2)^(1/3), mu=1/2 w^2 R^2
    """
    print("\n=== TEST 2: Estado fundamental vs Thomas-Fermi (g>0) ===")
    w, g, N_atoms = 1.0, 50.0, 1.0
    V = lambda x: 0.5 * w**2 * x**2
    sim = GPE1D(L=40.0, N=2048, g=g)
    sim.V = V(sim.x)

    psi0 = np.exp(-sim.x**2 / 2)  # gaussiana semilla
    psi = sim.ground_state(psi0, dtau=0.002, nsteps=40000, N_target=N_atoms)

    # Thomas-Fermi analitico
    R = (3 * N_atoms * g / (2 * w**2))**(1/3)
    mu_TF = 0.5 * w**2 * R**2
    n_TF = np.maximum(mu_TF - 0.5 * w**2 * sim.x**2, 0) / g

    n_num = np.abs(psi)**2
    mu_num = sim.chem_potential(psi)
    # error en el interior (|x| < 0.8 R, donde TF es valido)
    mask = np.abs(sim.x) < 0.8 * R
    err = np.sqrt(np.sum((n_num[mask] - n_TF[mask])**2) /
                  np.sum(n_TF[mask]**2))
    print(f"  Radio TF:              R = {R:.4f}")
    print(f"  Potencial quimico:     mu_num = {mu_num:.5f}  vs  mu_TF = {mu_TF:.5f}"
          f"   (err {abs(mu_num-mu_TF)/mu_TF:.2%})")
    print(f"  Error RMS perfil:      {err:.2e}  (interior |x|<0.8R)")

    plt.figure(figsize=(9, 4))
    plt.plot(sim.x, n_TF, 'k-', lw=3, alpha=0.4, label='Thomas-Fermi (analitico)')
    plt.plot(sim.x, n_num, 'r--', lw=1.3, label='numerico (tiempo imag.)')
    plt.xlabel('x'); plt.ylabel(r'$n(x)=|\psi|^2$'); plt.legend()
    plt.title(f'Estado fundamental vs Thomas-Fermi  (err={err:.1e}, mu err={abs(mu_num-mu_TF)/mu_TF:.1%})')
    plt.xlim(-R*1.5, R*1.5); plt.tight_layout()
    plt.savefig(f"{OUT}/test2_thomas_fermi.png", dpi=110)
    plt.close()
    return err < 0.05


if __name__ == "__main__":
    print("Validando solver GPE1D (split-step de Fourier)...")
    ok1 = test_bright_soliton()
    ok2 = test_thomas_fermi()
    print("\n" + "="*50)
    print(f"  TEST 1 (soliton):       {'PASA' if ok1 else 'FALLA'}")
    print(f"  TEST 2 (Thomas-Fermi):  {'PASA' if ok2 else 'FALLA'}")
    print("="*50)
    print("PNGs en out/")
    # Codigo de salida (auditoria fase 1, 2026-11-17): un test que falla debe fallar.
    # Tolerancias: las que ya usaba este archivo (soliton err RMS < 1e-2; TF err < 0.05).
    _sys.exit(0 if (ok1 and ok2) else 1)
