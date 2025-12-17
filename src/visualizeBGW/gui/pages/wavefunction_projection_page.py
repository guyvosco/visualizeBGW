"""
WavefunctionProjectionPage: GUI page for plotting wavefunction projections
and exporting them as .xsf files.

Data flow
---------
Inputs:
- WFN.h5 -> structure via read_structure(wfn_file)
- WFN.h5, band index, k-point -> (rho, mf_header) via calc_rho(wfn_file, bnd, kpt)

Options:
- Band index to plot (integer, 0 = first unoccupied band) with up/down arrows.
- k-point components: kx, ky, kz (floats in open text boxes).
- Replication counts (non-negative ints) in each direction:
    (-x, +x), (-y, +y), (-z, +z)

3D plotting (PyVista, new window, separate process):
- times_X = [n_minus_x, n_plus_x]
- times_Y = [n_minus_y, n_plus_y]
- times_Z = [n_minus_z, n_plus_z]
- plot_proj_wfn(pl, structure, rho, times_X, times_Y, times_Z)

Export:
- export_xsf(xsf_dir, rho, mf_header)

To avoid Qt/VTK event-loop issues and crashes, the 3D viewer
runs in a separate process that recomputes structure & rho.
"""

from __future__ import annotations

from typing import Callable, Optional, Tuple
import multiprocessing as mp

import numpy as np
import pyvista as pv

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from ..matplotlib_widget import MatplotlibWidget  # placeholder on the right
from ...io.berkeleygw_readers import read_structure
from ...analysis.berkeleygw_processing import calc_rho
from ...io.berkeleygw_exports import export_xsf
from ...plotting.berkeleygw_plots import plot_proj_wfn


# ----------------------------------------------------------------------
# Worker: PyVista in a separate process
# ----------------------------------------------------------------------


def _pv_worker_wavefunction_projection(
    wfn_path: str,
    band_index: int,
    kpt: Tuple[float, float, float],
    times_x: Tuple[int, int],
    times_y: Tuple[int, int],
    times_z: Tuple[int, int],
) -> None:
    """
    Worker function for plotting the wavefunction projection in a separate process.

    This isolates VTK/PyVista from the Qt main process and avoids crashes.
    """
    try:
        structure = read_structure(wfn_path)
        kpt_array = np.array(kpt, dtype=float)
        rho, mf_header = calc_rho(wfn_path, band_index, kpt_array)

        pl = pv.Plotter()
        # times_X, times_Y, times_Z are lists in the API, but we pass as lists
        plot_proj_wfn(
            pl,
            structure,
            rho,
            [times_x[0], times_x[1]],
            [times_y[0], times_y[1]],
            [times_z[0], times_z[1]],
        )
        pl.show(auto_close=True)
    except Exception as exc:  # noqa: BLE001
        print("[_pv_worker_wavefunction_projection] Error:", exc)


class WavefunctionProjectionPage(QWidget):
    """
    Page for wavefunction projection visualization and .xsf export.

    Parameters
    ----------
    on_back : callable, optional
        Callback invoked when "Back to menu" is pressed.
    """

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        on_back: Optional[Callable[[], None]] = None,
    ) -> None:
        super().__init__(parent)

        self._on_back = on_back

        # Cached rho/mf_header for export
        self._rho = None
        self._mf_header = None
        self._current_params: Optional[Tuple[str, int, Tuple[float, float, float]]] = None

        self._build_ui()

    # -------------------------
    # UI construction
    # -------------------------

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout()
        main_layout.setContentsMargins(16, 16, 16, 16)
        main_layout.setSpacing(10)

        # --- Top row: back button + title ---
        top_row = QHBoxLayout()
        top_row.setSpacing(8)

        back_button = QPushButton("Back to menu")
        back_button.setFixedWidth(120)
        back_button.clicked.connect(self._handle_back_clicked)

        title_label = QLabel("Wavefunction projection")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 18px; font-weight: bold;")

        top_row.addWidget(back_button, alignment=Qt.AlignLeft)
        top_row.addStretch(1)
        top_row.addWidget(title_label, stretch=0, alignment=Qt.AlignCenter)
        top_row.addStretch(1)

        main_layout.addLayout(top_row)

        # --- WFN file selector + status ---
        file_row = QHBoxLayout()
        file_row.setSpacing(8)

        file_label = QLabel("WFN.h5:")
        self.wfn_edit = QLineEdit()
        self.wfn_edit.setPlaceholderText("Select WFN.h5 …")
        self.wfn_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.wfn_edit.setMinimumWidth(450)

        browse_button = QPushButton("Browse…")
        browse_button.clicked.connect(self._browse_wfn)

        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.status_label.setMinimumWidth(200)

        file_row.addWidget(file_label)
        file_row.addWidget(self.wfn_edit)
        file_row.addWidget(browse_button)
        file_row.addStretch(1)
        file_row.addWidget(self.status_label)

        main_layout.addLayout(file_row)

        # --- Main area: left controls, right placeholder panel ---
        center_row = QHBoxLayout()
        center_row.setSpacing(10)

        # Left panel: options + buttons
        left_layout = QVBoxLayout()
        left_layout.setSpacing(10)

        # Band index (0 = first unoccupied band)
        band_row = QHBoxLayout()
        band_row.setSpacing(6)

        band_label = QLabel("band index (0 = first unoccupied):")
        self.band_down = QPushButton("▼")
        self.band_down.setFixedWidth(24)
        self.band_down.clicked.connect(self._decrement_band)

        self.band_edit = QLineEdit()
        self.band_edit.setText("0")
        self.band_edit.setMaximumWidth(60)

        self.band_up = QPushButton("▲")
        self.band_up.setFixedWidth(24)
        self.band_up.clicked.connect(self._increment_band)

        band_row.addWidget(band_label)
        band_row.addWidget(self.band_down)
        band_row.addWidget(self.band_edit)
        band_row.addWidget(self.band_up)
        band_row.addStretch(1)

        left_layout.addLayout(band_row)

        # k-point row: kx ky kz
        kpt_row = QHBoxLayout()
        kpt_row.setSpacing(6)

        kpt_label = QLabel("k-point (kx, ky, kz):")
        self.kx_edit = QLineEdit()
        self.kx_edit.setPlaceholderText("kx")
        self.kx_edit.setMaximumWidth(80)

        self.ky_edit = QLineEdit()
        self.ky_edit.setPlaceholderText("ky")
        self.ky_edit.setMaximumWidth(80)

        self.kz_edit = QLineEdit()
        self.kz_edit.setPlaceholderText("kz")
        self.kz_edit.setMaximumWidth(80)

        kpt_row.addWidget(kpt_label)
        kpt_row.addWidget(self.kx_edit)
        kpt_row.addWidget(self.ky_edit)
        kpt_row.addWidget(self.kz_edit)
        kpt_row.addStretch(1)

        left_layout.addLayout(kpt_row)

        # Replication controls
        repl_title = QLabel("Replication (number of times in each direction):")
        repl_title.setStyleSheet("font-weight: bold;")
        left_layout.addWidget(repl_title)

        x_row, self._minus_x_edit, self._plus_x_edit = self._create_axis_row("x")
        left_layout.addLayout(x_row)

        y_row, self._minus_y_edit, self._plus_y_edit = self._create_axis_row("y")
        left_layout.addLayout(y_row)

        z_row, self._minus_z_edit, self._plus_z_edit = self._create_axis_row("z")
        left_layout.addLayout(z_row)

        # Plot button
        self.plot_button = QPushButton("Plot wavefunction projection (PyVista)")
        self.plot_button.clicked.connect(self.plot_wavefunction_3d)
        left_layout.addWidget(self.plot_button)

        left_layout.addSpacing(8)

        # Export section
        export_title = QLabel("Export wavefunction to .xsf:")
        export_title.setStyleSheet("font-weight: bold;")
        left_layout.addWidget(export_title)

        export_row = QHBoxLayout()
        export_row.setSpacing(8)

        export_label = QLabel("Output directory:")
        self.export_dir_edit = QLineEdit()
        self.export_dir_edit.setPlaceholderText("Choose directory to save *.xsf …")
        self.export_dir_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.export_dir_edit.setMinimumWidth(350)

        export_browse = QPushButton("Browse…")
        export_browse.clicked.connect(self._browse_export_dir)

        export_row.addWidget(export_label)
        export_row.addWidget(self.export_dir_edit)
        export_row.addWidget(export_browse)

        left_layout.addLayout(export_row)

        self.export_button = QPushButton("Export .xsf")
        self.export_button.clicked.connect(self.export_to_xsf)
        left_layout.addWidget(self.export_button)

        left_layout.addStretch(1)

        center_row.addLayout(left_layout, stretch=0)

        # Right side: placeholder Matplotlib widget (keeps layout consistent)
        self.placeholder_plot = MatplotlibWidget(self)
        self.placeholder_plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        center_row.addWidget(self.placeholder_plot, stretch=1)

        main_layout.addLayout(center_row)

        self.setLayout(main_layout)

    # -------------------------
    # Axis row helpers
    # -------------------------

    def _create_axis_row(
        self, axis_label: str
    ) -> Tuple[QHBoxLayout, QLineEdit, QLineEdit]:
        """
        Create a row like:

            x:  (-) ▼ [0] ▲    (+) ▼ [0] ▲

        Returns the layout and the two QLineEdits (minus, plus).
        """
        row = QHBoxLayout()
        row.setSpacing(6)

        label = QLabel(f"{axis_label}:")
        row.addWidget(label)

        # Negative side
        minus_label = QLabel("−")
        minus_edit = QLineEdit()
        minus_edit.setText("0")
        minus_edit.setMaximumWidth(50)

        minus_down = QPushButton("▼")
        minus_down.setFixedWidth(24)
        minus_down.clicked.connect(
            lambda checked=False, le=minus_edit: self._decrement_int(le)
        )

        minus_up = QPushButton("▲")
        minus_up.setFixedWidth(24)
        minus_up.clicked.connect(
            lambda checked=False, le=minus_edit: self._increment_int(le)
        )

        row.addWidget(minus_label)
        row.addWidget(minus_down)
        row.addWidget(minus_edit)
        row.addWidget(minus_up)

        row.addSpacing(10)

        # Positive side
        plus_label = QLabel("+")
        plus_edit = QLineEdit()
        plus_edit.setText("0")
        plus_edit.setMaximumWidth(50)

        plus_down = QPushButton("▼")
        plus_down.setFixedWidth(24)
        plus_down.clicked.connect(
            lambda checked=False, le=plus_edit: self._decrement_int(le)
        )

        plus_up = QPushButton("▲")
        plus_up.setFixedWidth(24)
        plus_up.clicked.connect(
            lambda checked=False, le=plus_edit: self._increment_int(le)
        )

        row.addWidget(plus_label)
        row.addWidget(plus_down)
        row.addWidget(plus_edit)
        row.addWidget(plus_up)

        row.addStretch(1)

        return row, minus_edit, plus_edit

    # -------------------------
    # Integer helpers (replication)
    # -------------------------

    def _get_int(self, line_edit: QLineEdit) -> Optional[int]:
        txt = line_edit.text().strip()
        if not txt:
            return 0
        try:
            val = int(txt)
        except ValueError:
            self._set_status("Replication must be an integer")
            return None
        if val < 0:
            self._set_status("Replication must be ≥ 0")
            return None
        return val

    def _set_int(self, line_edit: QLineEdit, value: int) -> None:
        if value < 0:
            value = 0
        line_edit.setText(str(value))

    def _increment_int(self, line_edit: QLineEdit) -> None:
        val = self._get_int(line_edit)
        if val is None:
            return
        self._set_int(line_edit, val + 1)

    def _decrement_int(self, line_edit: QLineEdit) -> None:
        val = self._get_int(line_edit)
        if val is None:
            return
        if val > 0:
            self._set_int(line_edit, val - 1)

    # -------------------------
    # Band index helpers
    # -------------------------

    def _get_band_index(self) -> Optional[int]:
        txt = self.band_edit.text().strip()
        if not txt:
            return 0
        try:
            idx = int(txt)
        except ValueError:
            self._set_status("Band index must be an integer")
            return None
        if idx < 0:
            self._set_status("Band index must be ≥ 0")
            return None
        return idx

    def _set_band_index(self, idx: int) -> None:
        if idx < 0:
            idx = 0
        self.band_edit.setText(str(idx))

    def _increment_band(self) -> None:
        idx = self._get_band_index()
        if idx is None:
            return
        self._set_band_index(idx + 1)

    def _decrement_band(self) -> None:
        idx = self._get_band_index()
        if idx is None:
            return
        if idx > 0:
            self._set_band_index(idx - 1)

    # -------------------------
    # General helpers
    # -------------------------

    def _set_status(self, text: str) -> None:
        """Update the status label and process events so it appears immediately."""
        from PySide6.QtWidgets import QApplication

        self.status_label.setText(text)
        QApplication.processEvents()

    def _handle_back_clicked(self) -> None:
        if self._on_back is not None:
            self._on_back()

    def _browse_wfn(self) -> None:
        """Select WFN.h5 file."""
        dlg = QFileDialog(self, "Select WFN.h5 file")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setNameFilter("HDF5 files (*.h5 *.hdf5);;All files (*.*)")
        if dlg.exec():
            files = dlg.selectedFiles()
            if files:
                self.wfn_edit.setText(files[0])
                # Invalidate cached rho/mf_header
                self._rho = None
                self._mf_header = None
                self._current_params = None

    def _browse_export_dir(self) -> None:
        """Select directory for exporting .xsf."""
        dir_path = QFileDialog.getExistingDirectory(self, "Select export directory")
        if dir_path:
            self.export_dir_edit.setText(dir_path)

    # -------------------------
    # Common parameter parsing
    # -------------------------

    def _get_kpt(self) -> Optional[Tuple[float, float, float]]:
        kx_str = self.kx_edit.text().strip()
        ky_str = self.ky_edit.text().strip()
        kz_str = self.kz_edit.text().strip()

        if not (kx_str and ky_str and kz_str):
            self._set_status("All k-point components (kx, ky, kz) are required")
            return None

        try:
            kx = float(kx_str)
            ky = float(ky_str)
            kz = float(kz_str)
        except ValueError:
            self._set_status("k-point components must be floats")
            return None

        return (kx, ky, kz)

    def _get_replication_vectors(
        self,
    ) -> Optional[Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]]:
        nx_minus = self._get_int(self._minus_x_edit)
        nx_plus = self._get_int(self._plus_x_edit)
        ny_minus = self._get_int(self._minus_y_edit)
        ny_plus = self._get_int(self._plus_y_edit)
        nz_minus = self._get_int(self._minus_z_edit)
        nz_plus = self._get_int(self._plus_z_edit)

        if None in (nx_minus, nx_plus, ny_minus, ny_plus, nz_minus, nz_plus):
            return None

        return (nx_minus, nx_plus), (ny_minus, ny_plus), (nz_minus, nz_plus)

    # -------------------------
    # rho/mf_header loading (for export)
    # -------------------------

    def _load_rho_for_export(self, force: bool = False) -> bool:
        """
        Compute rho and mf_header for the current (wfn_path, band_index, kpt),
        caching the result for export_xsf.
        """
        wfn_path = self.wfn_edit.text().strip()
        if not wfn_path:
            self._set_status("WFN.h5 file is required")
            return False

        band_index = self._get_band_index()
        if band_index is None:
            return False

        kpt = self._get_kpt()
        if kpt is None:
            return False

        params = (wfn_path, band_index, kpt)

        if not force and self._current_params == params and self._rho is not None:
            self._set_status("rho already computed")
            return True

        try:
            self._set_status("computing rho...")
            rho, mf_header = calc_rho(wfn_path, band_index, np.array(kpt, dtype=float))
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"rho error: {exc}")
            print("[WavefunctionProjectionPage] Error computing rho:", exc)
            self._rho = None
            self._mf_header = None
            self._current_params = None
            return False

        self._rho = rho
        self._mf_header = mf_header
        self._current_params = params
        self._set_status("rho computed")
        return True

    # -------------------------
    # Actions
    # -------------------------

    def plot_wavefunction_3d(self) -> None:
        """
        Launch a separate PyVista process to show the wavefunction projection.
        """
        wfn_path = self.wfn_edit.text().strip()
        if not wfn_path:
            self._set_status("WFN.h5 file is required")
            return

        band_index = self._get_band_index()
        if band_index is None:
            return

        kpt = self._get_kpt()
        if kpt is None:
            return

        repl = self._get_replication_vectors()
        if repl is None:
            return
        times_x, times_y, times_z = repl

        try:
            self._set_status("launching PyVista viewer...")
            p = mp.Process(
                target=_pv_worker_wavefunction_projection,
                args=(wfn_path, band_index, kpt, times_x, times_y, times_z),
            )
            p.start()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[WavefunctionProjectionPage] Error spawning PyVista process:", exc)
            return

        self._set_status("3D viewer running")

    def export_to_xsf(self) -> None:
        """
        Export the current wavefunction projection to .xsf using export_xsf(xsf_dir, rho, mf_header).
        """
        if not self._load_rho_for_export(force=False):
            return

        if self._rho is None or self._mf_header is None:
            self._set_status("No rho data")
            return

        xsf_dir = self.export_dir_edit.text().strip()
        if not xsf_dir:
            self._set_status("Export directory is required")
            return

        try:
            self._set_status("exporting .xsf...")
            export_xsf(xsf_dir, self._rho, self._mf_header)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Export error: {exc}")
            print("[WavefunctionProjectionPage] Error exporting .xsf:", exc)
            return

        self._set_status("Exported .xsf")
