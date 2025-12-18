import numpy as np
from dummies import _dummy_mf_header, _dummy_structure
from visualizeBGW.io.berkeleygw_exports import export_structure, export_xsf


def test_export_structure_writes_poscar(tmp_path):
    structure = _dummy_structure()
    export_structure(tmp_path, structure)

    out_file = tmp_path / "structure.vasp"
    assert out_file.exists()
    content = out_file.read_text()

    # Very light sanity checks
    assert "Generated from visualizeBGW" in content or len(content.splitlines()) > 5
    assert "Si" in content or "14" in content  # depending on how you map Z→symbol


def test_export_xsf_writes_file(tmp_path):
    mf_header = _dummy_mf_header()
    # Simple 4x4x4 density
    rho = np.ones((4, 4, 4), dtype=float)

    export_xsf(tmp_path, rho, mf_header)

    out_file = tmp_path / "projwfn.xsf"
    assert out_file.exists()
    content = out_file.read_text()

    # Basic format checks
    assert "CRYSTAL" in content
    assert "PRIMVEC" in content
    assert "PRIMCOORD" in content
    assert "BEGIN_BLOCK_DATAGRID_3D" in content
