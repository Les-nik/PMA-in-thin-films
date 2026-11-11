import numpy as np
import pytest

from pma_spinwaves import comsol, precession
from pma_spinwaves.fourier import isofrequency_map
from pma_spinwaves.ultrathin import MU0, Film, dispersion_minimum, frequency_surface, omega

FILM = Film()


def test_kittel_limit():
    """At k = 0 and Q = 0: omega = gamma mu0 Ms sqrt(h (1 + h))."""
    B0 = 0.4
    h = B0 / (MU0 * FILM.Ms)
    expected = FILM.C * np.sqrt(h * (1 + h))
    assert omega(np.array([0.0]), B0, 0.0, film=FILM)[0] == pytest.approx(expected)


def test_sombrero_minimum_only_with_pma():
    k_min, _ = dispersion_minimum(0.4, 1.2, film=FILM)
    assert 30e6 < k_min < 45e6  # minimum at finite k for Q = 1.2
    k0, _ = dispersion_minimum(0.4, 0.0, film=FILM)
    assert k0 == 0.0  # no minimum without PMA


def test_damon_eshbach_above_backward_volume():
    k = np.array([30e6])
    assert omega(k, 0.4, 1.2, np.pi / 2, FILM)[0] > omega(k, 0.4, 1.2, 0.0, FILM)[0]


def test_surface_orientation():
    """Field along y: waves along x are Damon-Eshbach, so f(kx) > f(ky) at equal |k|."""
    kx, ky, f = frequency_surface(60e6, 0.4, 1.2, n=121, film=FILM)
    i0, iend = 60, 120
    assert f[i0, iend] > f[iend, i0]
    np.testing.assert_allclose(f, f[::-1, :])
    np.testing.assert_allclose(f, f[:, ::-1])


def test_comsol_bands():
    k, f = comsol.load_bands(comsol.DATA / "thick_20nm" / "dispersion" / "dispersion_B250mT_Q0.6_3bands.txt", 3)
    assert f.shape == (3, k.size)
    assert k[0] == 0 and f[0, 0] == pytest.approx(14.0)
    assert np.all(np.diff(k) > 0)


def test_angular_sweep_mirrored():
    kx, ky, f = comsol.load_angular_sweep(comsol.DATA / "angular" / "thick_20nm_B250mT_angles_0_to_pi2.txt")
    assert kx.size == ky.size == f.size == 527 * 4
    assert kx.min() == pytest.approx(-kx.max())


def test_handedness():
    assert precession.handedness(np.array([1.0]), np.array([1j]))[0] == 1
    assert precession.handedness(np.array([1.0]), np.array([-1j]))[0] == -1
    # k = 0 at 250 mT: the whole film precesses in the same (Larmor) sense
    _, mx = comsol.load_profile("mx", 250, 0)
    _, mz = comsol.load_profile("mz", 250, 0)
    s = precession.handedness(mx, mz)
    assert np.all(s == s[0])


def test_isofrequency_sign():
    """A wave moving along +x gives exp(-i k0 x) in a forward-FFT field; the map must peak at +k0."""
    n, dx = 128, 5e-9
    x = np.arange(n) * dx
    k0 = 2 * np.pi * 10 / (n * dx)
    field = np.exp(-1j * k0 * x)[None, :] * np.ones((64, 1))
    kx, ky, amp, _ = isofrequency_map(field, dx, dx)
    iy, ix = np.unravel_index(np.argmax(amp), amp.shape)
    assert kx[ix] == pytest.approx(k0)
    assert ky[iy] == pytest.approx(0, abs=1e-6)
