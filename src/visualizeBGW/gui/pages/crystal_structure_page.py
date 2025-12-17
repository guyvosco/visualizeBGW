"""
CrystalStructurePage: GUI page for plotting and exporting the crystal structure.

Data flow
---------
Input:
- WFN.h5 -> structure via read_structure(wfn_file)

Options:
- Replication counts (non-negative ints) in each direction:
    (-x, +x), (-y, +y), (-z, +z)
  Each pair is on one line, with text boxes and up/down arrows.

Plotting:
- times_x = [n_minus_x, n_plus_x]
- times_y = [n_minus_y, n_plus_y]
- times_z = [n_minus_z, n_plus_z]

plot_structure(pl, structure, times_x, times_y, times_z)
    where pl = pv.Plotter()

Export:
- export_structure(path, structure)
  where `path` is a directory where a POSCAR file named `structure.vasp`
  will be written (handled inside export_structure).

PyVista safety
-------------
To avoid Qt/VTK event-loop conflicts and segmentation faults, the
3D viewer runs in a separate process, which reloads the structure
from WFN.h5 in that process.
"""

from __future__ import annotations

from typing import Callable, Optional, Tuple
import multiprocessing as mp

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

# from ..matplotlib_widget import MatplotlibWidget  # just a placeholder panel on the right
from ...io.berkeleygw_readers import read_structure
from ...io.berkeleygw_exports import export_structure
from ...plotting.berkeleygw_plots import plot_structure


# ----------------------------------------------------------------------
# Worker function for PyVista in a separate process
# ----------------------------------------------------------------------


def _pv_worker_crystal_structure(
    wfn_path: str,
    times_x: Tuple[int, int],
    times_y: Tuple[int, int],
    times_z: Tuple[int, int],
) -> None:
    """
    Worker function for plotting the 3D crystal structure in a separate process.

    This prevents Qt/VTK issues from crashing the main GUI.
    """
    try:
        structure = read_structure(wfn_path)
        pl = pv.Plotter()
        plot_structure(pl, structure, list(times_x), list(times_y), list(times_z))
        pl.show(auto_close=True)
    except Exception as exc:  # noqa: BLE001
        print("[_pv_worker_crystal_structure] Error:", exc)


class CrystalStructurePage(QWidget):
    """
    Page for crystal structure visualization and export.

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

        # Cached structure (for export; plotting uses separate process)
        self._structure = None
        self._current_wfn_file: Optional[str] = None

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

        title_label = QLabel("Crystal structure")
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
        self.status_label.setMinimumWidth(180)

        file_row.addWidget(file_label)
        file_row.addWidget(self.wfn_edit)
        file_row.addWidget(browse_button)
        file_row.addStretch(1)
        file_row.addWidget(self.status_label)

        main_layout.addLayout(file_row)

        # --- Main area: left controls, right placeholder panel ---
        center_row = QHBoxLayout()
        center_row.setSpacing(10)

        # Left panel: replication options + buttons
        left_layout = QVBoxLayout()
        left_layout.setSpacing(10)

        # Replication controls title
        repl_title = QLabel("Replication (number of times in each direction):")
        repl_title.setStyleSheet("font-weight: bold;")
        left_layout.addWidget(repl_title)

        # x line: (-x and +x)
        x_row, self._minus_x_edit, self._plus_x_edit = self._create_axis_row("x")
        left_layout.addLayout(x_row)

        # y line: (-y and +y)
        y_row, self._minus_y_edit, self._plus_y_edit = self._create_axis_row("y")
        left_layout.addLayout(y_row)

        # z line: (-z and +z)
        z_row, self._minus_z_edit, self._plus_z_edit = self._create_axis_row("z")
        left_layout.addLayout(z_row)

        # Plot button
        self.plot_button = QPushButton("Plot structure (PyVista)")
        self.plot_button.clicked.connect(self.plot_structure_3d)
        left_layout.addWidget(self.plot_button)

        left_layout.addSpacing(8)

        # Export section
        export_title = QLabel("Export structure as POSCAR (structure.vasp):")
        export_title.setStyleSheet("font-weight: bold;")
        left_layout.addWidget(export_title)

        export_row = QHBoxLayout()
        export_row.setSpacing(8)

        export_label = QLabel("Output directory:")
        self.export_dir_edit = QLineEdit()
        self.export_dir_edit.setPlaceholderText("Choose directory to save structure.vasp …")
        self.export_dir_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.export_dir_edit.setMinimumWidth(350)

        export_browse = QPushButton("Browse…")
        export_browse.clicked.connect(self._browse_export_dir)

        export_row.addWidget(export_label)
        export_row.addWidget(self.export_dir_edit)
        export_row.addWidget(export_browse)

        left_layout.addLayout(export_row)

        self.export_button = QPushButton("Export structure")
        self.export_button.clicked.connect(self.export_structure_to_poscar)
        left_layout.addWidget(self.export_button)

        left_layout.addStretch(1)

        center_row.addLayout(left_layout, stretch=0)

        # # Right side: placeholder Matplotlib widget just to keep layout consistent
        # self.placeholder_plot = MatplotlibWidget(self)
        # self.placeholder_plot.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)
        # center_row.addWidget(self.placeholder_plot, stretch=1)

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

        # Negative side group
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

        # Positive side group
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
    # Integer helpers
    # -------------------------

    def _get_int(self, line_edit: QLineEdit) -> Optional[int]:
        """Parse a non-negative integer from a line edit."""
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
                # Invalidate cached structure
                self._structure = None
                self._current_wfn_file = None

    def _browse_export_dir(self) -> None:
        """Select directory for exporting structure.vasp."""
        dir_path = QFileDialog.getExistingDirectory(self, "Select export directory")
        if dir_path:
            self.export_dir_edit.setText(dir_path)

    # -------------------------
    # Data loading & actions
    # -------------------------

    def _load_structure(self, force: bool = False) -> bool:
        """Load structure from WFN.h5 if needed (for export / sanity checks)."""
        wfn_path = self.wfn_edit.text().strip()
        if not wfn_path:
            self._set_status("WFN.h5 file is required")
            return False

        if (
            not force
            and self._current_wfn_file == wfn_path
            and self._structure is not None
        ):
            self._set_status("Structure loaded")
            return True

        try:
            self._set_status("reading structure...")
            structure = read_structure(wfn_path)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Error: {exc}")
            print("[CrystalStructurePage] Error reading structure:", exc)
            self._structure = None
            self._current_wfn_file = None
            return False

        self._structure = structure
        self._current_wfn_file = wfn_path
        self._set_status("Structure loaded")
        return True

    def _get_replication_vectors(
        self,
    ) -> Optional[Tuple[Tuple[int, int], Tuple[int, int], Tuple[int, int]]]:
        """Read replication counts from widgets and return (times_x, times_y, times_z)."""
        nx_minus = self._get_int(self._minus_x_edit)
        nx_plus = self._get_int(self._plus_x_edit)
        ny_minus = self._get_int(self._minus_y_edit)
        ny_plus = self._get_int(self._plus_y_edit)
        nz_minus = self._get_int(self._minus_z_edit)
        nz_plus = self._get_int(self._plus_z_edit)

        if None in (nx_minus, nx_plus, ny_minus, ny_plus, nz_minus, nz_plus):
            # Status already set by _get_int
            return None

        times_x = (nx_minus, nx_plus)
        times_y = (ny_minus, ny_plus)
        times_z = (nz_minus, nz_plus)
        return times_x, times_y, times_z

    def plot_structure_3d(self) -> None:
        """Launch a separate process with a PyVista window showing the replicated structure."""
        # We only use _load_structure here to validate the file and make sure it exists.
        if not self._load_structure(force=False):
            return
        if self._structure is None:
            self._set_status("No structure data")
            return

        repl = self._get_replication_vectors()
        if repl is None:
            return
        times_x, times_y, times_z = repl

        wfn_path = self.wfn_edit.text().strip()
        try:
            self._set_status("launching PyVista viewer...")
            p = mp.Process(
                target=_pv_worker_crystal_structure,
                args=(wfn_path, times_x, times_y, times_z),
            )
            p.start()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[CrystalStructurePage] Error while spawning structure plot process:", exc)
            return

        self._set_status("3D viewer running")

    def export_structure_to_poscar(self) -> None:
        """Export structure as a POSCAR file (structure.vasp) in the chosen directory."""
        if not self._load_structure(force=False):
            return
        if self._structure is None:
            self._set_status("No structure data")
            return

        dir_path = self.export_dir_edit.text().strip()
        if not dir_path:
            self._set_status("Export directory is required")
            return

        try:
            self._set_status("exporting structure...")
            export_structure(dir_path, self._structure)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Export error: {exc}")
