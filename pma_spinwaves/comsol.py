"""Readers for the COMSOL text exports in ``data/comsol``.

COMSOL writes a few header lines starting with ``%`` and then two columns,
``x`` and ``Height``. What the columns mean depends on the export:

* dispersion: k (rad/um) and f (GHz). Several bands are written for each k,
  one after another, so row ``i`` belongs to band ``i % n_bands``.
* angular sweep: the same, but the repeated rows are different directions of
  k (0 to pi/2 in steps of pi/32) instead of different bands.
* profiles and torques: z (nm, film centred at 0) and the value.
"""

from __future__ import annotations

from pathlib import Path

import numpy as np

DATA = Path(__file__).resolve().parents[1] / "data" / "comsol"


def load_txt(path: str | Path) -> np.ndarray:
    """Two-column COMSOL export as an (n, 2) array."""
    return np.loadtxt(path, comments="%")


def load_bands(path: str | Path, n_bands: int) -> tuple[np.ndarray, np.ndarray]:
    """Split an interleaved dispersion export into bands.

    Returns ``k`` of shape (nk,) and ``f`` of shape (n_bands, nk).
    """
    data = load_txt(path)
    if len(data) % n_bands:
        raise ValueError(f"{path}: {len(data)} rows is not a multiple of {n_bands}")
    k = data[::n_bands, 0]
    f = data[:, 1].reshape(-1, n_bands).T
    return k, f


def load_angular_sweep(path: str | Path, n_angles: int = 17,
                       angle_step: float = np.pi / 32):
    """Frequency on scattered (kx, ky) points, mirrored to all four quadrants.

    The export holds one band for directions 0, pi/32, ..., pi/2 measured from
    the kx axis. The film is symmetric under kx -> -kx and ky -> -ky, so the
    other quadrants are copies. The result is meant for
    ``matplotlib.pyplot.tricontour(kx, ky, f)``.
    """
    k, f = load_bands(path, n_angles)  # here "bands" are directions
    angles = np.arange(n_angles) * angle_step
    kx, ky, ff = [], [], []
    for a, fa in zip(angles, f):
        for sx, sy in [(1, 1), (-1, 1), (1, -1), (-1, -1)]:
            kx.append(sx * k * np.cos(a))
            ky.append(sy * k * np.sin(a))
            ff.append(fa)
    return np.concatenate(kx), np.concatenate(ky), np.concatenate(ff)


def load_profile(component: str, B_mT: int, k: int, folder: Path = DATA / "thick_20nm" / "profiles"):
    """Complex thickness profile of m_x or m_z: returns z (nm) and m(z)."""
    re = load_txt(folder / f"{component}_B{B_mT}mT_k{k}.txt")
    im = load_txt(folder / f"{component}_imag_B{B_mT}mT_k{k}.txt")
    return re[:, 0], re[:, 1] + 1j * im[:, 1]


TORQUE_TERMS = ("exch", "pma", "zeeman", "demag", "eff")


def load_torques(B_mT: int, k: int, component: str = "tx",
                 folder: Path = DATA / "thick_20nm" / "torques") -> tuple[np.ndarray, dict]:
    """Torque contributions across the thickness: z (nm) and {term: values}."""
    z, terms = None, {}
    for term in TORQUE_TERMS:
        data = load_txt(folder / f"{component}_{term}_B{B_mT}_k{k}.txt")
        z, terms[term] = data[:, 0], data[:, 1]
    return z, terms
