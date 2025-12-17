import numpy as np
from ase.data import chemical_symbols
from ..analysis.units import bohr2angstrom

def export_structure(path, structure):
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

    with open(path + "/structure.vasp", "w") as f:
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
    alat = mf_header['crystal']['alat']
    avec = mf_header['crystal']['avec']
    nat = mf_header['crystal']['nat']
    atyp = mf_header['crystal']['atyp']
    apos = mf_header['crystal']['apos']
    fftgrid = mf_header['gspace']['FFTgrid']

    with open(xsf_dir + "/projwfn.xsf", 'w') as f:
        f.write(' CRYSTAL\n')
        f.write(' PRIMVEC\n')
        for vec in avec:
            vec *= alat * bohr2angstrom
            f.write('    %1.9f    %1.9f    %1.9f\n' % (vec[0], vec[1], vec[2]))
        f.write(' PRIMCOORD\n')
        f.write('          %i          1\n' % nat)
        for i, pos in enumerate(apos):
            pos *= alat * bohr2angstrom
            f.write('%s         %1.9f    %1.9f    %1.9f\n' % (chemical_symbols[atyp[i]], pos[0], pos[1], pos[2]))
        f.write('BEGIN_BLOCK_DATAGRID_3D\n')
        f.write('3D_PWSCF\n')
        f.write('BEGIN_DATAGRID_3D_UNKNOWN\n')
        f.write('          %i          %i          %i\n' % (fftgrid[0] + 1, fftgrid[1] + 1, fftgrid[2] + 1))
        f.write('  0.000000  0.000000  0.000000\n')
        for vec in avec:
            f.write('    %1.9f    %1.9f    %1.9f\n' % (vec[0], vec[1], vec[2]))
        i = 0
        for z in list(range(fftgrid[2])) + [0]:
            for y in list(range(fftgrid[1])) + [0]:
                for x in list(range(fftgrid[0])) + [0]:
                    f.write('  %0.6e' % rho[x, y, z])
                    i += 1
                    if i % 6 == 0:
                        f.write('\n')
        f.write('\nEND_DATAGRID_3D\n')
        f.write('END_BLOCK_DATAGRID_3D')
