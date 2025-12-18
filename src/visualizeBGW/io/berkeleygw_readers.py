import os
from collections import defaultdict
import numpy as np
import h5py


def read_h5_group(group):
    """
    Input:
        group (h5py.Group) — HDF5 group object

    Output:
        Nested dictionary of NumPy arrays

    Description:
        Recursively extracts a full HDF5 group into a Python dictionary.
        ! Note: Use only for small to medium-sized datasets! Useful for reading metadata.
    """
    result = {}
    for key, item in group.items():
        if isinstance(item, h5py.Dataset):
            result[key] = item[()]
        elif isinstance(item, h5py.Group):
            result[key] = read_h5_group(item)
    return result


def read_structure(wfn_file):
    """
    Input:
        wfn_file (str) — path to WFN.h5

    Output:
        structure (dict) — alat, lattice vectors, atomic species, fractional positions
        mf_header (dict) — parsed metadata from the HDF5 file

    Description:
        Reads structural and metadata information from BerkeleyGW’s WFN.h5.
    """
    if not os.path.isfile(wfn_file):
        raise FileNotFoundError(f"File not found: {wfn_file}")

    with h5py.File(wfn_file, "r") as f:
        mf_header = read_h5_group(f["mf_header"])
        alat = mf_header["crystal"]["alat"]
        lattice_vectors = mf_header["crystal"]["avec"]
        atom_types = mf_header["crystal"]["atyp"]
        atom_positions = mf_header["crystal"]["apos"]
        structure = {
            "alat": alat,
            "lattice_vectors": lattice_vectors,
            "atom_types": atom_types,
            "atom_positions": atom_positions,
        }

    return structure


def read_chi_converge(chi_converge_file):
    """
    Input:
        path (str) — path to chi_converge.dat

    Output:
        dict of convergence data per q-point

    Description:
        Parses BerkeleyGW chi_converge.dat and returns convergence information.
    """
    if not os.path.isfile(chi_converge_file):
        raise FileNotFoundError(f"File not found: {chi_converge_file}")

    chi_converge_data = defaultdict(
        lambda: {
            "ncbands": [],
            "chi(0,0)": [],
            "chi(0,0)-extrap": 0.0,
            "chi(G_max,G_max)": [],
            "chi(G_max,G_max)-extrap": 0.0,
        }
    )
    with open(chi_converge_file, "r", encoding="utf-8") as f:
        for line in f:
            if "q=" in line:
                qpt = tuple(float(i) for i in line.split()[2:5])
            elif "#" in line or line.strip() == "":
                continue
            else:
                data = line.split()
                chi_converge_data[qpt]["ncbands"].append(int(data[0]))
                chi_converge_data[qpt]["chi(0,0)"].append(float(data[1]))
                chi_converge_data[qpt]["chi(G_max,G_max)"].append(float(data[3]))
                if chi_converge_data[qpt]["chi(0,0)-extrap"] == 0.0:
                    chi_converge_data[qpt]["chi(0,0)-extrap"] = float(data[2])
                if chi_converge_data[qpt]["chi(G_max,G_max)-extrap"] == 0.0:
                    chi_converge_data[qpt]["chi(G_max,G_max)-extrap"] = float(data[4])

    return chi_converge_data


def read_ch_converge(ch_converge_file):
    """
    Input:
        path (str) — path to ch_converge.dat

    Output:
        dict with Coulomb-hole term vs number of bands

    Description:
        Reads the COHSEX Coulomb-hole convergence output file.
    """
    if not os.path.isfile(ch_converge_file):
        raise FileNotFoundError(f"File not found: {ch_converge_file}")

    ch_converge_data = defaultdict(
        lambda: {
            "nbands": [],
            "CH(vbm)": [],
            "CH(cbm)": [],
            "diff": [],
            "CH(vbm)-extrap": 0.0,
            "CH(cbm)-extrap": 0.0,
            "diff-extrap": 0.0,
            "CH(vbm)-error": 0.0,
            "CH(cbm)-error": 0.0,
            "diff-error": 0.0,
        }
    )
    with open(ch_converge_file, "r", encoding="utf-8") as f:
        for line in f:
            if "k =" in line:
                qpt = tuple(float(i) for i in line.split()[3:6])
            elif "1/N extrap" in line:
                data = line.split()
                ch_converge_data[qpt]["CH(vbm)-extrap"] = float(data[8])
                ch_converge_data[qpt]["CH(cbm)-extrap"] = float(data[9])
                ch_converge_data[qpt]["diff-extrap"] = float(data[10])
            elif "Error est." in line:
                data = line.split()
                ch_converge_data[qpt]["CH(vbm)-error"] = float(data[8])
                ch_converge_data[qpt]["CH(cbm)-error"] = float(data[9])
                ch_converge_data[qpt]["diff-error"] = float(data[10])
            elif "#" in line or line.strip() == "":
                continue
            else:
                data = line.split()
                ch_converge_data[qpt]["nbands"].append(int(data[0]))
                ch_converge_data[qpt]["CH(vbm)"].append(float(data[1]))
                ch_converge_data[qpt]["CH(cbm)"].append(float(data[2]))
                ch_converge_data[qpt]["diff"].append(float(data[3]))

    return ch_converge_data


def read_eqp(eqp_file):
    """
    Input:
        path (str) — path to eqp.dat

    Output:
        kpts (np.ndarray) — array of k-points
        emf (np.ndarray) — array of mean-field energies
        eqp (np.ndarray) — array of quasiparticle energies

    Description:
        Reads BerkeleyGW eqp.dat file and extracts k-points, mean-field energies, and quasiparticle energies.
    """
    if not os.path.isfile(eqp_file):
        raise FileNotFoundError(f"File not found: {eqp_file}")

    with open(eqp_file, "r", encoding="utf-8") as f:
        lines = f.readlines()
        nbands = int(lines[0].split()[-1])
        nk = len(lines) // (nbands + 1)
        kpts = np.zeros((nk, 3), dtype=float)
        emf = np.zeros((nk, nbands), dtype=float)
        eqp = np.zeros((nk, nbands), dtype=float)
        for k in range(nk):
            kpts[k, :] = [float(i) for i in lines[k * (nbands + 1)].split()[:3]]
            for b in range(nbands):
                data_line = lines[k * (nbands + 1) + b + 1].split()
                emf[k, b] = float(data_line[2])
                eqp[k, b] = float(data_line[3])

    return kpts, emf, eqp


def read_eigenvalues_noeh(eigenvalues_noeh_file):
    """
    Input:
        path (str) — path to eigenvalues_noeh.dat

    Output:
        eigs (np.ndarray) — array of eigenvalues without electron-hole interaction
        dipole (np.ndarray) — array of complex dipole matrix elements

    Description:
        Reads BerkeleyGW eigenvalues_noeh.dat file and extracts eigenvalues and dipole matrix elements.
    """
    if not os.path.isfile(eigenvalues_noeh_file):
        raise FileNotFoundError(f"File not found: {eigenvalues_noeh_file}")

    eigs = np.loadtxt(eigenvalues_noeh_file, skiprows=4, usecols=6, dtype=float)
    raw_data = np.loadtxt(
        eigenvalues_noeh_file, skiprows=4, usecols=(8, 9), dtype=float
    )
    dipole = raw_data[:, 0] + 1j * raw_data[:, 1]

    return eigs, dipole


def read_eigenvalues(eigenvalues_file):
    """
    Input:
        path (str) — path to eigenvalues.dat

    Output:
        eigs (np.ndarray) — array of eigenvalues with electron-hole interaction
        dipole (np.ndarray) — array of complex dipole matrix elements

    Description:
        Reads BerkeleyGW eigenvalues.dat file and extracts eigenvalues and dipole matrix elements.
    """
    if not os.path.isfile(eigenvalues_file):
        raise FileNotFoundError(f"File not found: {eigenvalues_file}")

    eigs = np.loadtxt(eigenvalues_file, skiprows=4, usecols=0, dtype=float)
    raw_data = np.loadtxt(eigenvalues_file, skiprows=4, usecols=(2, 3), dtype=float)
    dipole = raw_data[:, 0] + 1j * raw_data[:, 1]

    return eigs, dipole
