import numpy as np
import pytest
from dummies import _dummy_structure, _dummy_kpts

from visualizeBGW.analysis.berkeleygw_processing import (
    p_in_kpts,
    p2p,
    get_kpts_in_plane,
    gaussian,
    lorentzian,
    get_absorption_spectra,
    get_grid_and_levels,
    calc_rho,
    get_eps_head,
    get_eigenvectors_components,
)

# ---------- p_in_kpts ----------


def test_p_in_kpts_finds_point():
    kpts = _dummy_kpts()
    idx = p_in_kpts(kpts, [0.5, 0.5, 0.0])
    assert idx == 6  # see simple_kpts fixture in conftest.py


def test_p_in_kpts_raises_if_not_found():
    kpts = _dummy_kpts()
    with pytest.raises(ValueError):
        p_in_kpts(kpts, [0.1, 0.2, 0.3])


# ---------- p2p ----------


def test_p2p_picks_points_on_segment_and_sorts():
    # Points along x-axis between 0 and 1, plus an off-line point
    kpts = _dummy_kpts()
    p1 = [0.0, 0.0, 0.0]
    p2 = [1.0, 1.0, 1.0]

    indices = p2p(kpts, p1, p2)

    assert np.array_equal(indices, np.array([0, 7]))


def test_p2p_raises_for_identical_points():
    kpts = np.zeros((3, 3))
    with pytest.raises(ValueError):
        p2p(kpts, [0.0, 0.0, 0.0], [0.0, 0.0, 0.0])


# ---------- get_kpts_in_plane ----------


def test_get_kpts_in_plane_aggregates_quantities():
    # Two distinct projections in kx-ky plane: (0,0) and (1,0)
    kpts = _dummy_kpts()
    quantity = np.array([1.0, 2.0, 3.0, 4.0, 5.0, 6.0, 7.0, 8.0])
    dims = [0, 1]

    k2d, q2d = get_kpts_in_plane(kpts, quantity, dims)

    assert k2d.shape == (4, 2)
    assert np.allclose(k2d, np.array([[0.0, 0.0], [0.0, 0.5], [0.5, 0.0], [0.5, 0.5]]))
    assert q2d.shape == (4,)
    assert np.all(q2d == np.array([1.0 + 2.0, 3.0 + 4.0, 5.0 + 6.0, 7.0 + 8.0]))


# ---------- gaussian / lorentzian ----------


def test_gaussian_is_normalized_numerically():
    eta = 0.1
    x = np.linspace(-10 * eta, 10 * eta, 50001)
    dx = x[1] - x[0]
    vals = gaussian(x, eta)
    integral = np.sum(vals) * dx
    assert np.isclose(integral, 1.0, atol=1e-3)


def test_lorentzian_is_normalized_numerically():
    eta = 0.1
    # Need a wider range for Lorentzian tails
    x = np.linspace(-100 * eta, 100 * eta, 50001)
    dx = x[1] - x[0]
    vals = lorentzian(x, eta)
    integral = np.sum(vals) * dx
    assert np.isclose(integral, 1.0, atol=5e-3)


# ---------- get_absorption_spectra ----------


@pytest.mark.parametrize("function", ["lorentzian", "gaussian", "voigt"])
def test_get_absorption_spectra_peak_near_transition_and_normalized(function):
    eigs = np.array([1.0])  # one excitation at 1 eV
    dipole = np.array([2.0 + 0.0j])
    broadening = 0.1

    energy, spectra = get_absorption_spectra(
        eigs=eigs,
        dipole=dipole,
        broadening=broadening,
        function=function,
        de=0.001,
    )

    # Normalized
    assert np.isclose(spectra.max(), 1.0, atol=1e-6)

    # Peak should be close to the excitation energy
    peak_E = energy[np.argmax(spectra)]
    assert np.isclose(peak_E, 1.0, atol=5 * broadening)


# ---------- get_grid_and_levels ----------


def test_get_grid_and_levels_returns_valid_grid_and_levels():
    rho = np.zeros((4, 4, 4), dtype=float)
    rho[1, 1, 1] = 1.0
    rho[2, 2, 2] = 0.5

    structure = _dummy_structure()
    grid, levels = get_grid_and_levels(
        rho, structure, n_levels=3, p_low=0.5, p_high=0.9
    )

    # Grid has same dimensions as rho
    assert list(grid.dimensions) == list(rho.shape)
    # Data array exists and has the right length
    assert "rho" in grid.array_names
    assert grid["rho"].shape[0] == rho.size

    # Levels are in ascending order and positive (since rho >= 0)
    assert np.all(np.diff(levels) >= 0)
    assert np.all(levels >= 0.0)


def test_get_grid_and_levels_handles_degenerate_case():
    # All zeros -> amp.size == 0 -> should fall back to a single level (0.0)
    rho = np.ones((2, 2, 2), dtype=float)
    structure = _dummy_structure()
    grid, levels = get_grid_and_levels(rho, structure, n_levels=3)

    assert list(grid.dimensions) == list(rho.shape)
    assert levels.shape == (1,)
    assert levels[0] == 1.0


# ---------- calc_rho / get_eps_head / get_eigenvectors_components (error paths) ----------


def test_calc_rho_raises_for_missing_file(tmp_path):
    missing = tmp_path / "does_not_exist.h5"
    with pytest.raises(FileNotFoundError):
        calc_rho(str(missing), bnd=0, kpt=[0.0, 0.0, 0.0])


def test_get_eps_head_raises_for_missing_file(tmp_path):
    missing = tmp_path / "epsmat_missing.h5"
    with pytest.raises(FileNotFoundError):
        get_eps_head(str(missing))


def test_get_eigenvectors_components_raises_for_missing_file(tmp_path):
    missing = tmp_path / "eigenvectors_missing.h5"
    with pytest.raises(FileNotFoundError):
        get_eigenvectors_components(str(missing), index=0)
