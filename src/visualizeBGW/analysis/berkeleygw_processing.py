import os
import numpy as np
from scipy.special import voigt_profile
import h5py
import pyvista as pv
from ..io.berkeleygw_readers import read_h5_group


def p_in_kpts(kpts, p, tol=1e-6):
    """
    Input:
        kpts (np.ndarray): of shape (N, 3), list of k-points
        p (list or np.ndarray): of shape (3,), point to find in kpts
        tol (float): tolerance for matching points - default 1e-6

    Output:
        idx (int): index of point p in kpts
    """
    p = np.array(p)
    diffs = np.linalg.norm(kpts - p, axis=1)
    idx = np.where(diffs <= tol)[0]
    if len(idx) == 0:
        raise ValueError(f"Point {p} not found in k-points")
    return idx[0]


def p2p(kpts, p1, p2, tol=1e-6):
    """
    Input:
        kpts (np.ndarray): of shape (N, 3), list of k-points
        p1 (list or np.ndarray): of shape (3,), starting point of the line segment
        p2 (list or np.ndarray): of shape (3,), ending point of the line segment
        tol (float): tolerance for distance to the line segment - default 1e-6

    Output:
        online (np.ndarray): of shape (M,), indices of k-points that lie on the line segment from p1 to p2, ordered along the segment
    """
    p1 = np.array(p1)
    p2 = np.array(p2)
    d = p2 - p1
    d2 = np.dot(d, d.T)
    if d2 == 0:
        raise ValueError("p1 and p2 must be different points")

    proj = np.einsum("ij,j->i", kpts - p1, d) / d2
    dist = np.linalg.norm(kpts - (p1 + proj[:, None] * d), axis=1)
    mask = (dist <= tol) & (proj >= -tol) & (proj <= 1 + tol)

    idx = np.where(mask)[0]
    order = np.argsort(proj[idx])
    online = idx[order]

    return online


def get_kpts_in_plane(kpts, quantity, dims):
    """
    Input:
        kpts (np.ndarray): of shape (N, 3), list of k-points
        quantity (np.ndarray): of shape (N,), quantity defined at each k-point
        dims (list or np.ndarray): of shape (2,), indices of the two dimensions defining the plane (e.g., [0, 1] for kx-ky plane)

    Output:
        kpts_in_plane (np.ndarray): of shape (M, 2), unique k-points projected onto the specified 2D plane
        quantity_in_plane (np.ndarray): of shape (M,), summed quantity at each unique k-point in the plane

    Description:
        Projects k-points onto the specified 2D plane and sums the quantity for k-points that map to the same point in that plane.
    """
    kpts_in_plane = np.unique(kpts[:, dims], axis=0)
    quantity_in_plane = np.array([], dtype=quantity.dtype)
    for unique_kpt in kpts_in_plane:
        indices = np.argwhere(np.all(np.isclose(kpts[:, dims], unique_kpt), axis=1))
        quantity_in_plane = np.append(quantity_in_plane, np.sum(quantity[indices]))

    return kpts_in_plane, quantity_in_plane


def calc_rho(wfn_file, bnd, kpt):
    """
    Input:
        wfn_file (str): path to the wavefunction HDF5 file
        bnd (int): band index
        kpt (list or np.ndarray): of shape (3,), k-point at which to calculate the charge density

    Output:
        rho (np.ndarray): charge density grid
        mf_header (dict): metadata from the wavefunction file

    Description:
        Computes the mean-field wavefunction density on a 3D grid in real space.
    """
    if not os.path.isfile(wfn_file):
        raise FileNotFoundError(f"File not found: {wfn_file}")

    with h5py.File(wfn_file, "r") as f:
        mf_header = read_h5_group(f["mf_header"])
        kpts = mf_header["kpoints"]["rk"]
        ns = mf_header["kpoints"]["nspin"] * mf_header["kpoints"]["nspinor"]
        fftgrid = mf_header["gspace"]["FFTgrid"]
        ngk = mf_header["kpoints"]["ngk"]
        occ = np.sum(mf_header["kpoints"]["occ"], axis=-1)
        if not np.allclose(occ, occ[0, 0]):
            raise ValueError("Occupations vary with k-point, cannot proceed")
        try:
            occ = int(occ[0, 0])
        except Exception as exc:
            raise ValueError(
                "Occupation number is not an integer, cannot proceed"
            ) from exc
        k = p_in_kpts(kpts, kpt)
        ngk_offset = np.sum(ngk[:k])
        gvecs = f["wfns/gvecs"][ngk_offset : ngk_offset + ngk[k], :]
        coeffs = f["wfns/coeffs"][bnd + occ, :, ngk_offset : ngk_offset + ngk[k], :]
    coeffs = coeffs[..., 0] + coeffs[..., 1] * 1j

    rho = np.zeros((fftgrid[0], fftgrid[1], fftgrid[2]), dtype=np.float64)
    for s in range(ns):
        fftbox = np.zeros((fftgrid[0], fftgrid[1], fftgrid[2]), dtype=np.complex128)
        fftbox[tuple(gvecs.T)] = coeffs[s, :]
        phi = np.fft.ifftn(fftbox, norm=None)
        rho += np.square(np.abs(phi)) * np.prod(fftgrid)

    return rho, mf_header


def get_grid_and_levels(rho, structure, n_levels=1, p_low=0.98, p_high=0.995):
    """
    Input:
        rho (np.ndarray): 3D charge density grid
        structure (dict): structure information containing lattice vectors
        n_levels (int): number of isosurface levels to generate - default 1
        p_low (float): lower percentile for isosurface levels - default 0.98
        p_high (float): upper percentile for isosurface levels - default 0.995

    Output:
        grid (pyvista.StructuredGrid): structured grid with charge density data
        levels (np.ndarray): array of isosurface levels

    Description:
        Creates a structured grid from the charge density and computes isosurface levels based on specified percentiles.
    """
    Gx, Gy, Gz = rho.shape
    a1, a2, a3 = structure["lattice_vectors"]
    FX, FY, FZ = np.meshgrid(
        np.arange(Gx) / Gx, np.arange(Gy) / Gy, np.arange(Gz) / Gz, indexing="ij"
    )
    grid = pv.StructuredGrid(
        FX * a1[0] + FY * a2[0] + FZ * a3[0],
        FX * a1[1] + FY * a2[1] + FZ * a3[1],
        FX * a1[2] + FY * a2[2] + FZ * a3[2],
    )
    grid["rho"] = rho.ravel(order="F").astype(np.float64)

    amp = rho.ravel()
    amp = amp[amp > 0]
    v_low = np.quantile(amp, p_low)
    v_high = np.quantile(amp, p_high)
    if not np.isfinite(v_low) or not np.isfinite(v_high) or v_high <= v_low:
        v_high = float(np.max(amp)) if amp.size else 0.0
        return grid, np.array([v_high])

    return grid, np.linspace(v_low, v_high, n_levels)


def get_eps_head(epsmat_file):
    """
    Input:
        epsmat_file (str): path to the epsilon matrix HDF5 file

    Output:
        epshead (list): head of the dielectric function at each q-point
        qlen (list): length of each q-point in the bdot metric

    Description:
        Extracts the head of the dielectric function and the lengths of q-points from the epsilon
    """
    if not os.path.isfile(epsmat_file):
        raise FileNotFoundError(f"File not found: {epsmat_file}")

    with h5py.File(epsmat_file, "r") as f:
        mf_header = read_h5_group(f["mf_header"])
        bdot = mf_header["crystal"]["bdot"]
        gvec = mf_header["gspace"]["components"]
        get_heaader = read_h5_group(f["eps_header"])
        qpts = get_heaader["qpoints"]["qpts"]
        gind_eps2rho = get_heaader["gspace"]["gind_eps2rho"] - 1
        nmtx = get_heaader["gspace"]["nmtx"]
        matrix = f["mats/matrix"]
        epshead = np.zeros((len(qpts),), dtype=float)
        qlen = np.zeros((len(qpts),), dtype=float)
        for iq, qpt in enumerate(qpts):
            qlen[iq] = float(np.sqrt(qpt @ (bdot @ qpt)))
            for ig_eps in range(nmtx[iq]):
                ig_rho = gind_eps2rho[iq, ig_eps]
                if ig_rho < 0 or ig_rho >= gvec.shape[0]:
                    continue
                if np.all(gvec[ig_rho, :] == 0):
                    epshead[iq] = matrix[iq, 0, 0, ig_eps, ig_eps, 0]
                    break

    qlen, epshead = zip(*sorted(zip(qlen, epshead)))

    return epshead, qlen


def gaussian(x, eta):
    return np.exp(-np.square(x) / (2 * np.square(eta))) / (np.sqrt(2 * np.pi) * eta)


def lorentzian(x, eta):
    return (eta / 2) / (np.pi * (np.square(x) + np.square(eta / 2)))


def get_absorption_spectra(eigs, dipole, broadening, function, de=0.001):
    """
    Input:
        eigs (np.ndarray): of shape (N,), excitation energies
        dipole (np.ndarray): of shape (N,), dipole matrix elements
        broadening (float): broadening parameter
        function (str): type of broadening function ('lorentzian', 'voigt', or 'gaussian')
        de (float): energy step for the spectra - default 0.001 eV

    Output:
        energy (np.ndarray): energy grid
        spectra (np.ndarray): normalized absorption spectra

    Description:
        Computes the absorption spectra using the specified broadening function.
    """
    energy = np.arange(0, np.max(eigs) + 3 * broadening, de)
    spectra = np.zeros_like(energy)
    for i, omega in enumerate(energy):
        if function == "lorentzian":
            spectra[i] = np.sum(
                lorentzian(omega - eigs, broadening) * np.abs(dipole) ** 2
            )
            spectra[i] -= np.sum(
                lorentzian(-omega - eigs, broadening) * np.abs(dipole) ** 2
            )
        elif function == "voigt":
            spectra[i] = np.sum(
                voigt_profile(omega - eigs, broadening, broadening)
                * np.abs(dipole) ** 2
            )
            spectra[i] -= np.sum(
                voigt_profile(-omega - eigs, broadening, broadening)
                * np.abs(dipole) ** 2
            )
        else:
            spectra[i] = np.sum(
                gaussian(omega - eigs, broadening) * np.abs(dipole) ** 2
            )
            spectra[i] -= np.sum(
                gaussian(-omega - eigs, broadening) * np.abs(dipole) ** 2
            )

    spectra /= np.max(spectra)
    return energy, spectra


def get_eigenvectors_components(eigenvectors_file, index):
    """
    Input:
        eigenvectors_file (str): path to the eigenvectors HDF5 file
        index (int): index of the exciton state to analyze

    Output:
        Ak (np.ndarray): of shape (Nk,), k-point contribution to the exciton state
        Ac (np.ndarray): of shape (Nc,), conduction band contribution to the exciton state
        Av (np.ndarray): of shape (Nv,), valence band contribution to the exciton state
        exciton_header (dict): metadata from the exciton file
        mf_header (dict): metadata from the mean-field file

    Description:
        Decomposes the exciton eigenvector into contributions from k-points, conduction bands, and valence bands.
    """
    if not os.path.isfile(eigenvectors_file):
        raise FileNotFoundError(f"File not found: {eigenvectors_file}")

    with h5py.File(eigenvectors_file, "r") as f:
        exciton_header = read_h5_group(f["exciton_header"])
        mf_header = read_h5_group(f["mf_header"])
        A = f["exciton_data/eigenvectors"][0, index, :, :, :, :, :]
    A = A[:, :, :, :, :, 0] + 1j * A[:, :, :, :, :, 1]
    Ak = np.sum(np.square(np.abs(A)), axis=(2, 3, 4))
    Ac = np.sum(np.square(np.abs(A)), axis=(1, 3, 4))
    Av = np.sum(np.square(np.abs(A)), axis=(1, 2, 4))

    return Ak, Ac, Av, exciton_header, mf_header
