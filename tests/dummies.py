import numpy as np


def _dummy_structure():
    return {
        "alat": 10.0,
        "avec": np.eye(3),  # primitive lattice in Bohr (will be scaled)
        "lattice_vectors": np.eye(3),
        "atom_types": np.array([14, 14]),  # Silicon
        "atom_positions": np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
    }


def _dummy_mf_header():
    # Minimal header structure used by export_xsf
    return {
        "crystal": {
            "alat": 10.0,
            "avec": np.eye(3),
            "nat": 2,
            "atyp": np.array([14, 14]),
            "apos": np.array([[0.0, 0.0, 0.0], [0.25, 0.25, 0.25]]),
        },
        "gspace": {
            "FFTgrid": np.array([4, 4, 4]),
        },
    }


def _dummy_kpts():
    return np.array(
        [
            [0.0, 0.0, 0.0],
            [0.0, 0.0, 0.5],
            [0.0, 0.5, 0.0],
            [0.0, 0.5, 0.5],
            [0.5, 0.0, 0.0],
            [0.5, 0.0, 0.5],
            [0.5, 0.5, 0.0],
            [0.5, 0.5, 0.5],
        ]
    )
