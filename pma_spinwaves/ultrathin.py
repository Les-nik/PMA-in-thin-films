"""Analytical spin-wave dispersion of an ultrathin film with perpendicular anisotropy.

The film lies in the (x, y) plane and is magnetized in-plane by the field H0.
For a film much thinner than the exchange length the magnetization is uniform
across the thickness, and linearizing the Landau-Lifshitz equation gives two
stiffness frequencies, one for each precession component:

    omega_oop = C [h + (k l_ex)^2 + (1 - P) - Q]       out-of-plane component
    omega_ip  = C [h + (k l_ex)^2 + P sin^2(phi)]      in-plane component

    omega = sqrt(omega_oop * omega_ip)

with

    C    = gamma mu0 Ms
    h    = H0 / Ms
    Q    = 2 Ku / (mu0 Ms^2)       (reduced perpendicular anisotropy)
    l_ex = sqrt(2 A / (mu0 Ms^2))  (exchange length)
    P    = 1 - (1 - exp(-k d)) / (k d)
    phi  = angle between k and the static magnetization

phi = 0 is the backward-volume geometry (k parallel to M) and phi = pi/2 the
Damon-Eshbach geometry (k perpendicular to M).

What PMA does: Q enters only omega_oop and lowers it. The term (1 - P) is the
out-of-plane demagnetizing field, which also decreases with k. When Q is
close to 1 + h, omega_oop is small at k = 0 and the exchange term pulls it up
again at large k, so omega(k) has a minimum at finite k: the "sombrero"
dispersion of the paper.
"""

from __future__ import annotations

from dataclasses import dataclass

import numpy as np

MU0 = 4e-7 * np.pi


@dataclass(frozen=True)
class Film:
    """Material and thickness. Defaults: the 2 nm CoFeB film of the paper."""

    Ms: float = 1344e3        # saturation magnetization, A/m
    Aex: float = 13.6e-12     # exchange stiffness, J/m
    thickness: float = 2e-9   # m
    gamma: float = 176e9      # gyromagnetic ratio, rad / (s T)

    @property
    def exchange_length(self) -> float:
        return np.sqrt(2 * self.Aex / (MU0 * self.Ms**2))

    @property
    def C(self) -> float:
        """gamma mu0 Ms in rad/s."""
        return self.gamma * MU0 * self.Ms

    def h(self, B0: float) -> float:
        """Reduced field H0 / Ms for an applied field B0 = mu0 H0 in tesla."""
        return B0 / (MU0 * self.Ms)


def _P(kd: np.ndarray) -> np.ndarray:
    """Thin-film dipolar factor P(kd) = 1 - (1 - exp(-|kd|)) / |kd|, with P(0) = 0."""
    kd = np.abs(np.asarray(kd, dtype=float))
    out = np.zeros_like(kd)
    nz = kd > 0
    out[nz] = 1 - (1 - np.exp(-kd[nz])) / kd[nz]
    return out


def omega_components(k, B0: float, Q: float, phi, film: Film = Film()):
    """The two stiffness frequencies (rad/s): out-of-plane and in-plane.

    ``k`` in rad/m, ``phi`` in radians (broadcast against ``k``).
    """
    k = np.asarray(k, dtype=float)
    P = _P(k * film.thickness)
    h = film.h(B0)
    exch = (k * film.exchange_length) ** 2
    omega_oop = film.C * (h + exch + (1 - P) - Q)
    omega_ip = film.C * (h + exch + P * np.sin(phi) ** 2)
    return omega_oop, omega_ip


def omega(k, B0: float, Q: float, phi=np.pi / 2, film: Film = Film()) -> np.ndarray:
    """Angular frequency (rad/s). NaN where omega_oop * omega_ip < 0 (unstable state)."""
    w_oop, w_ip = omega_components(k, B0, Q, phi, film)
    product = w_oop * w_ip
    with np.errstate(invalid="ignore"):
        return np.where(product >= 0, np.sqrt(np.abs(product)), np.nan)


def frequency_surface(kmax: float, B0: float, Q: float, n: int = 301, film: Film = Film()):
    """f(kx, ky) in Hz on a square grid, with the field (and M) along y.

    With M along y, phi is the angle between k and the y axis, so
    sin(phi) = kx / |k|.
    """
    k1 = np.linspace(-kmax, kmax, n)
    kx, ky = np.meshgrid(k1, k1)
    km = np.hypot(kx, ky)
    phi = np.arctan2(np.abs(kx), np.abs(ky))  # angle from the y axis, in radians
    return kx, ky, omega(km, B0, Q, phi, film) / (2 * np.pi)


def dispersion_minimum(B0: float, Q: float, phi=np.pi / 2, film: Film = Film(),
                       kmax: float = 300e6, n: int = 20001):
    """Wavenumber (rad/m) and frequency (Hz) of the lowest point of f(k).

    Returns k = 0 when the dispersion has no minimum at finite k.
    """
    k = np.linspace(0, kmax, n)
    f = omega(k, B0, Q, phi, film) / (2 * np.pi)
    i = int(np.nanargmin(f))
    return k[i], f[i]
