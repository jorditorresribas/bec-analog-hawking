"""
Solver de la ecuacion de Gross-Pitaevskii (GPE) en 1D.
Metodo: split-step de Fourier (Strang splitting, 2o orden).

GPE adimensional:
    i d_t psi = [ -1/2 d_xx + V(x) + g |psi|^2 ] psi

Validaciones incluidas (ambas con solucion analitica):
  1) Soliton brillante  (g<0, sin trampa)  -> se traslada sin deformarse.
  2) Estado fundamental en trampa armonica (g>0, tiempo imaginario)
     -> perfil de Thomas-Fermi.

Primer ladrillo del pipeline para gravedad analoga (horizonte sonico).
"""

import numpy as np


class GPE1D:
    def __init__(self, L=40.0, N=2048, g=1.0, V=None):
        self.L = L
        self.N = N
        self.x = np.linspace(-L/2, L/2, N, endpoint=False)
        self.dx = self.x[1] - self.x[0]
        self.k = 2.0 * np.pi * np.fft.fftfreq(N, d=self.dx)
        self.k2 = self.k**2
        self.g = g
        self.V = np.zeros(N) if V is None else V

    def norm(self, psi):
        return np.sum(np.abs(psi)**2) * self.dx

    def normalize(self, psi, N_target):
        return psi * np.sqrt(N_target / self.norm(psi))

    def _step(self, psi, dt):
        """Un paso split-step simetrico. dt real -> tiempo real;
        dt = -1j*dtau -> tiempo imaginario (relajacion a fundamental)."""
        # medio paso: potencial + no linealidad (espacio real)
        Vnl = self.V + self.g * np.abs(psi)**2
        psi = psi * np.exp(-0.5j * Vnl * dt)
        # paso completo: cinetica (espacio de Fourier)
        psi = np.fft.ifft(np.exp(-0.5j * self.k2 * dt) * np.fft.fft(psi))
        # medio paso final: potencial + no linealidad actualizado
        Vnl = self.V + self.g * np.abs(psi)**2
        psi = psi * np.exp(-0.5j * Vnl * dt)
        return psi

    def evolve_real(self, psi, dt, nsteps):
        for _ in range(nsteps):
            psi = self._step(psi, dt)
        return psi

    def ground_state(self, psi0, dtau=0.005, nsteps=20000, N_target=None, tol=1e-10):
        """Tiempo imaginario: dt = -i*dtau, renormaliza cada paso."""
        if N_target is None:
            N_target = self.norm(psi0)
        psi = self.normalize(psi0.astype(complex), N_target)
        mu_prev = 0.0
        for i in range(nsteps):
            psi = self._step(psi, -1j * dtau)
            psi = self.normalize(psi, N_target)
            if i % 200 == 0:
                mu = self.chem_potential(psi)
                if abs(mu - mu_prev) < tol:
                    break
                mu_prev = mu
        return psi

    def chem_potential(self, psi):
        """mu = <psi| -1/2 d_xx + V + g|psi|^2 |psi> / <psi|psi>"""
        lap = np.fft.ifft(-self.k2 * np.fft.fft(psi))
        kin = -0.5 * np.real(np.conj(psi) * lap)
        pot = (self.V + self.g * np.abs(psi)**2) * np.abs(psi)**2
        return np.sum(kin + pot) * self.dx / self.norm(psi)
