import os
import numpy as np
from ase.data import chemical_symbols
from ..analysis.units import bohr2angstrom


def export_structure(path, structure):
    """
    Input:
        path (str) — directory
        structure (dict) — lattice & atoms

    Output:
        Writes structure.vasp in POSCAR format

    Description:
        Exports structural geometry to a VASP-compatible file.
    """
    uniques = []
    counts = []
    last_elem = None
    for Z in structure["atom_types"]:
        elem = chemical_symbols[Z]
        if elem == last_elem:
            counts[-1] += 1
        else:
            uniques.append(elem)
            counts.append(1)
        last_elem = elem

    prefix = " ".join(f"{e}{n}" for e, n in zip(uniques, counts))

    with open(os.path.join(path, "structure.vasp"), "w", encoding="utf-8") as f:
        f.write(prefix + "\n")
        f.write(" 1.0\n")

        for vec in structure["lattice_vectors"]:
            x, y, z = vec * structure["alat"] * bohr2angstrom
            f.write(f"  {x:.10f}  {y:.10f}  {z:.10f}\n")

        unique_elems = []
        counts = []
        for Z in structure["atom_types"]:
            elem = chemical_symbols[Z]
            if elem not in unique_elems:
                unique_elems.append(elem)
                counts.append(1)
            else:
                counts[unique_elems.index(elem)] += 1

        f.write("  " + "  ".join(unique_elems) + "\n")
        f.write("  " + "  ".join(map(str, counts)) + "\n")

        f.write("Direct\n")

        inv_lat = np.linalg.inv(structure["lattice_vectors"])
        for pos in structure["atom_positions"]:
            fx, fy, fz = pos @ inv_lat
            f.write(f"  {fx:.10f}  {fy:.10f}  {fz:.10f}\n")


def export_xsf(xsf_dir, rho, mf_header):
    """
    Input:
        xsf_dir (str) — directory
        rho (ndarray) — scalar field
        mf_header (dict) — structural metadata

    Output:
        Writes projwfn.xsf file

    Description:
        Saves the 3D real-space wavefunction density into XSF format for XCrySDen.
    """
    alat = mf_header["crystal"]["alat"]
    avec = mf_header["crystal"]["avec"]
    nat = mf_header["crystal"]["nat"]
    atyp = mf_header["crystal"]["atyp"]
    apos = mf_header["crystal"]["apos"]
    fftgrid = mf_header["gspace"]["FFTgrid"]

    with open(os.path.join(xsf_dir, "projwfn.xsf"), "w", encoding="utf-8") as f:
        f.write(" CRYSTAL\n")
        f.write(" PRIMVEC\n")
        for vec in avec:
            vec_ang = vec * alat * bohr2angstrom
            f.write(
                f"    {vec_ang[0]:1.9f}    {vec_ang[1]:1.9f}    {vec_ang[2]:1.9f}\n"
            )
        f.write(" PRIMCOORD\n")
        f.write(f"          {nat}          1\n")
        for i, pos in enumerate(apos):
            pos_ang = pos * alat * bohr2angstrom
            f.write(
                f"{chemical_symbols[atyp[i]]}         {pos_ang[0]:1.9f}    {pos_ang[1]:1.9f}    {pos_ang[2]:1.9f}\n"
            )
        f.write("BEGIN_BLOCK_DATAGRID_3D\n")
        f.write("3D_PWSCF\n")
        f.write("BEGIN_DATAGRID_3D_UNKNOWN\n")
        f.write(
            f"          {fftgrid[0] + 1}          {fftgrid[1] + 1}          {fftgrid[2] + 1}\n"
        )
        f.write("  0.000000  0.000000  0.000000\n")
        for vec in avec:
            vec_ang = vec * alat * bohr2angstrom
            f.write(
                f"    {vec_ang[0]:1.9f}    {vec_ang[1]:1.9f}    {vec_ang[2]:1.9f}\n"
            )
        i = 0
        for z in list(range(fftgrid[2])) + [0]:
            for y in list(range(fftgrid[1])) + [0]:
                for x in list(range(fftgrid[0])) + [0]:
                    f.write(f"  {rho[x, y, z]:0.6e}")
                    i += 1
                    if i % 6 == 0:
                        f.write("\n")
        f.write("\nEND_DATAGRID_3D\n")
        f.write("END_BLOCK_DATAGRID_3D")
