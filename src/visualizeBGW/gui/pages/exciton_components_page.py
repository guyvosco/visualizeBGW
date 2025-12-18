"""
ExcitonComponentsPage: GUI page for plotting exciton components and 3D k-space.

Data flow
---------
Inputs:
- eigenvalues.dat   -> eigs, dipole
- eigenvectors.h5   -> Ak, Ac, Av, exciton_header, mf_header (for a given index)

Functions:
- read_eigenvalues(eigenvalues_file) -> eigs, dipole
- get_eigenvectors_components(eigenvectors_file, np.array([index]))
    -> Ak, Ac, Av, exciton_header, mf_header

Plotting:
- plot_eigenvector_components(
      ax,
      index,
      Ak, Ac, Av,
      exciton_header, mf_header,
      eigs[index], dipole[index],
  )

3D Brillouin zone:
- plot_eigenvector_k_3d(pl, Ak, exciton_header, mf_header)
  where pl = pv.Plotter()

To avoid PyVista/Qt conflicts and segfaults, the 3D plot is done
in a separate process.
"""

from __future__ import annotations

from typing import Callable, Optional
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

from ..matplotlib_widget import MatplotlibWidget
from ...io.berkeleygw_readers import read_eigenvalues
from ...analysis.berkeleygw_processing import get_eigenvectors_components
from ...plotting.berkeleygw_plots import (
    plot_eigenvector_components,
    plot_eigenvector_k_3d,
)


def _pv_worker_plot_eigenvector_k_3d(Ak, exciton_header, mf_header) -> None:
    """
    Worker function for plotting the 3D Brillouin zone in a separate process.

    This isolates VTK/PyVista from the Qt main process and avoids crashes
    or event-loop conflicts bringing down the GUI.
    """
    try:
        pl = pv.Plotter()
        plot_eigenvector_k_3d(pl, Ak, exciton_header, mf_header)
        pl.show(auto_close=True)
    except Exception as exc:  # noqa: BLE001
        print("[_pv_worker_plot_eigenvector_k_3d] Error:", exc)


class ExcitonComponentsPage(QWidget):
    """
    Page for plotting exciton components and 3D Brillouin-zone distribution.

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

        # Cached eigenvalues
        self._eigs: Optional[np.ndarray] = None
        self._dipole: Optional[np.ndarray] = None
        self._current_eig_file: Optional[str] = None

        # We only store the path for eigenvectors; we read them per index
        self._current_vec_file: Optional[str] = None

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

        title_label = QLabel("Exciton components")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 18px; font-weight: bold;")

        top_row.addWidget(back_button, alignment=Qt.AlignLeft)
        top_row.addStretch(1)
        top_row.addWidget(title_label, stretch=0, alignment=Qt.AlignCenter)
        top_row.addStretch(1)

        main_layout.addLayout(top_row)

        # --- File selectors + status ---

        # eigenvectors.h5
        vec_row = QHBoxLayout()
        vec_row.setSpacing(8)

        vec_label = QLabel("eigenvectors.h5:")
        self.vec_file_edit = QLineEdit()
        self.vec_file_edit.setPlaceholderText("Select eigenvectors.h5 …")
        self.vec_file_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.vec_file_edit.setMinimumWidth(450)

        vec_browse = QPushButton("Browse…")
        vec_browse.clicked.connect(
            lambda: self._browse_file(self.vec_file_edit, is_eig=False)
        )

        vec_row.addWidget(vec_label)
        vec_row.addWidget(self.vec_file_edit)
        vec_row.addWidget(vec_browse)
        vec_row.addStretch(1)

        main_layout.addLayout(vec_row)

        # eigenvalues.dat + status
        eig_row = QHBoxLayout()
        eig_row.setSpacing(8)

        eig_label = QLabel("eigenvalues.dat:")
        self.eig_file_edit = QLineEdit()
        self.eig_file_edit.setPlaceholderText("Select eigenvalues.dat …")
        self.eig_file_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.eig_file_edit.setMinimumWidth(450)

        eig_browse = QPushButton("Browse…")
        eig_browse.clicked.connect(
            lambda: self._browse_file(self.eig_file_edit, is_eig=True)
        )

        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.status_label.setMinimumWidth(180)

        eig_row.addWidget(eig_label)
        eig_row.addWidget(self.eig_file_edit)
        eig_row.addWidget(eig_browse)
        eig_row.addStretch(1)
        eig_row.addWidget(self.status_label)

        main_layout.addLayout(eig_row)

        # --- Main area: left options, right plot ---
        center_row = QHBoxLayout()
        center_row.setSpacing(10)

        # Left panel: exciton index + buttons
        left_layout = QVBoxLayout()
        left_layout.setSpacing(10)

        # Exciton index control: ▼ [index] ▲
        idx_row = QHBoxLayout()
        idx_row.setSpacing(4)

        idx_label = QLabel("exciton index:")
        self.idx_down_button = QPushButton("▼")
        self.idx_down_button.setFixedWidth(28)
        self.idx_down_button.clicked.connect(self._decrement_index)

        self.idx_edit = QLineEdit()
        self.idx_edit.setPlaceholderText("0")
        self.idx_edit.setMaximumWidth(60)
        self.idx_edit.setText("0")

        self.idx_up_button = QPushButton("▲")
        self.idx_up_button.setFixedWidth(28)
        self.idx_up_button.clicked.connect(self._increment_index)

        idx_row.addWidget(idx_label)
        idx_row.addWidget(self.idx_down_button)
        idx_row.addWidget(self.idx_edit)
        idx_row.addWidget(self.idx_up_button)
        idx_row.addStretch(1)

        left_layout.addLayout(idx_row)

        # Plot buttons – vertical: Plot components above Plot 3D
        self.plot_button = QPushButton("Plot components")
        self.plot_button.clicked.connect(self.plot_components)

        self.plot_3d_button = QPushButton("Plot 3D Brillouin zone")
        self.plot_3d_button.clicked.connect(self.plot_3d_brillouin)

        left_layout.addWidget(self.plot_button)
        left_layout.addWidget(self.plot_3d_button)

        left_layout.addStretch(1)

        # Right: Matplotlib plot
        self.plot_widget = MatplotlibWidget(self)
        self.plot_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        center_row.addLayout(left_layout, stretch=0)
        center_row.addWidget(self.plot_widget, stretch=1)

        main_layout.addLayout(center_row)

        self.setLayout(main_layout)

    # -------------------------
    # Helper methods
    # -------------------------

    def _set_status(self, text: str) -> None:
        """Update the status label and process events so it appears immediately."""
        from PySide6.QtWidgets import QApplication

        self.status_label.setText(text)
        QApplication.processEvents()

    def _handle_back_clicked(self) -> None:
        if self._on_back is not None:
            self._on_back()

    def _browse_file(self, line_edit: QLineEdit, is_eig: bool) -> None:
        """Open a file dialog and store the result in line_edit."""
        dlg = QFileDialog(self, "Select file")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setNameFilter("All files (*.*)")
        if dlg.exec():
            files = dlg.selectedFiles()
            if files:
                line_edit.setText(files[0])
                if is_eig:
                    self._eigs = self._dipole = None
                    self._current_eig_file = None
                else:
                    self._current_vec_file = files[0]

    # ---- exciton index controls ----

    def _get_index(self) -> Optional[int]:
        """Parse exciton index from the text box; must be integer >= 0."""
        idx_str = self.idx_edit.text().strip()
        if not idx_str:
            return 0
        try:
            idx = int(idx_str)
        except ValueError:
            self._set_status("Index must be an integer")
            return None
        if idx < 0:
            self._set_status("Index must be >= 0")
            return None
        return idx

    def _set_index(self, idx: int) -> None:
        if idx < 0:
            idx = 0
        self.idx_edit.setText(str(idx))

    def _increment_index(self) -> None:
        idx = self._get_index()
        if idx is None:
            return
        self._set_index(idx + 1)

    def _decrement_index(self) -> None:
        idx = self._get_index()
        if idx is None:
            return
        if idx > 0:
            self._set_index(idx - 1)

    # -------------------------
    # Data loading
    # -------------------------

    def _load_eigenvalues(self, force: bool = False) -> bool:
        file_path = self.eig_file_edit.text().strip()
        if not file_path:
            self._set_status("eigenvalues.dat file is required")
            return False

        if (
            not force
            and self._current_eig_file == file_path
            and self._eigs is not None
            and self._dipole is not None
        ):
            self._set_status("Eigenvalues loaded")
            return True

        try:
            self._set_status("reading eigenvalues...")
            eigs, dipole = read_eigenvalues(file_path)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Error: {exc}")
            print("[ExcitonComponentsPage] Error reading eigenvalues:", exc)
            self._eigs = self._dipole = None
            self._current_eig_file = None
            return False

        self._eigs, self._dipole = eigs, dipole
        self._current_eig_file = file_path
        self._set_status("Eigenvalues loaded")
        return True

    def _load_eigenvectors_for_index(self, index: int):
        """
        Load eigenvector components for a given exciton index.

        Returns
        -------
        Ak, Ac, Av, exciton_header, mf_header
        or None on failure.
        """
        vec_path = self.vec_file_edit.text().strip()
        if not vec_path:
            self._set_status("eigenvectors.h5 file is required")
            return None

        try:
            self._set_status("reading eigenvectors...")
            Ak, Ac, Av, exciton_header, mf_header = get_eigenvectors_components(
                vec_path, np.array([index])
            )
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Eigenvectors error: {exc}")
            print("[ExcitonComponentsPage] Error reading eigenvectors:", exc)
            return None

        return Ak, Ac, Av, exciton_header, mf_header

    # -------------------------
    # Plotting
    # -------------------------

    def plot_components(self) -> None:
        """Plot exciton components in the Matplotlib axes."""
        # Parse index
        index = self._get_index()
        if index is None:
            return

        # Ensure eigenvalues loaded
        if not self._load_eigenvalues(force=False):
            return

        if self._eigs is None or self._dipole is None:
            self._set_status("No eigenvalues data")
            return

        if index >= len(self._eigs):
            self._set_status("Index out of range for eigenvalues")
            return

        # Load eigenvectors for this index
        ev_data = self._load_eigenvectors_for_index(index)
        if ev_data is None:
            return
        Ak, Ac, Av, exciton_header, mf_header = ev_data

        # Clear the entire figure and create a fresh Axes
        fig = self.plot_widget.canvas.figure
        fig.clear()
        ax = fig.add_subplot(111)
        self.plot_widget.canvas.axes = ax  # keep attribute consistent

        try:
            self._set_status("plotting components...")
            plot_eigenvector_components(
                ax,
                index,
                Ak,
                Ac,
                Av,
                exciton_header,
                mf_header,
                self._eigs[index],
                self._dipole[index],
            )
            self.plot_widget.canvas.draw()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[ExcitonComponentsPage] Error while plotting components:", exc)
            return

        self._set_status("Done")

    def plot_3d_brillouin(self) -> None:
        """
        Open a PyVista window with the k-space distribution for the current exciton.

        This is done in a separate process to avoid Qt/VTK event-loop issues
        and potential segmentation faults crashing the GUI.
        """
        # Parse index
        index = self._get_index()
        if index is None:
            return

        # Load eigenvectors for this index
        ev_data = self._load_eigenvectors_for_index(index)
        if ev_data is None:
            return
        Ak, _Ac, _Av, exciton_header, mf_header = ev_data

        try:
            self._set_status("launching 3D k-space viewer...")
            p = mp.Process(
                target=_pv_worker_plot_eigenvector_k_3d,
                args=(Ak, exciton_header, mf_header),
            )
            p.start()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"3D plot error: {exc}")
            print("[ExcitonComponentsPage] Error while spawning 3D plot process:", exc)
            return

        self._set_status("3D viewer running")
