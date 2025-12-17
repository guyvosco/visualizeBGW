"""
DielectricFunctionPage: GUI page for plotting the head of the dielectric matrix.

Data workflow
-------------
We read TWO files:

1. eps0mat.h5  (unscreened / bare or first q?)
2. epsmat.h5   (screened / second q?)

For each file we:
    - open it with h5py.File
    - call get_eps_head(epsmat_file)  ->  epshead_i, qlen_i

Then we append (concatenate) results:
    - epshead = concat([epshead_1, epshead_2], axis=0)
    - qlen   = concat([qlen_1,   qlen_2],   axis=0)

Plotting options
----------------
We have two modes:

1. "plot the screening function"
2. "plot inverse dielectric matrix"

Both call:

    plot_eps_head(ax, qlen, epshead, plot_inverse)

where:
    plot_inverse = False for option 1
    plot_inverse = True  for option 2
"""

from __future__ import annotations

from typing import Callable, Optional, Tuple

import numpy as np
import h5py

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QButtonGroup,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QRadioButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from ..matplotlib_widget import MatplotlibWidget
from ...analysis.berkeleygw_processing import get_eps_head
from ...plotting.berkeleygw_plots import plot_eps_head


class DielectricFunctionPage(QWidget):
    """
    Page for plotting the dielectric function (screening function or inverse ε).

    Parameters
    ----------
    on_back : callable, optional
        Callback invoked when "Back to menu" is pressed.
        Typical usage: on_back() → main_window.show_page("menu").
    """

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        on_back: Optional[Callable[[], None]] = None,
    ) -> None:
        super().__init__(parent)

        self._on_back = on_back

        # Cached data
        self._epshead: Optional[np.ndarray] = None
        self._qlen: Optional[np.ndarray] = None
        self._current_files: Optional[Tuple[str, str]] = None

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

        title_label = QLabel("Dielectric function")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 18px; font-weight: bold;")

        top_row.addWidget(back_button, alignment=Qt.AlignLeft)
        top_row.addStretch(1)
        top_row.addWidget(title_label, stretch=0, alignment=Qt.AlignCenter)
        top_row.addStretch(1)

        main_layout.addLayout(top_row)

        # --- File selectors + status ---
        # First row: eps0mat.h5
        eps0_row = QHBoxLayout()
        eps0_row.setSpacing(8)

        eps0_label = QLabel("eps0mat.h5 file:")
        self.eps0_edit = QLineEdit()
        self.eps0_edit.setPlaceholderText("Select eps0mat.h5 …")
        self.eps0_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.eps0_edit.setMinimumWidth(450)

        eps0_browse = QPushButton("Browse…")
        eps0_browse.clicked.connect(lambda: self._browse_file(self.eps0_edit))

        eps0_row.addWidget(eps0_label)
        eps0_row.addWidget(self.eps0_edit)
        eps0_row.addWidget(eps0_browse)
        eps0_row.addStretch(1)

        main_layout.addLayout(eps0_row)

        # Second row: epsmat.h5 + status on the right
        eps_row = QHBoxLayout()
        eps_row.setSpacing(8)

        eps_label = QLabel("epsmat.h5 file:")
        self.eps_edit = QLineEdit()
        self.eps_edit.setPlaceholderText("Select epsmat.h5 …")
        self.eps_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.eps_edit.setMinimumWidth(450)

        eps_browse = QPushButton("Browse…")
        eps_browse.clicked.connect(lambda: self._browse_file(self.eps_edit))

        # Status
        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.status_label.setMinimumWidth(160)

        eps_row.addWidget(eps_label)
        eps_row.addWidget(self.eps_edit)
        eps_row.addWidget(eps_browse)
        eps_row.addStretch(1)
        eps_row.addWidget(self.status_label)

        main_layout.addLayout(eps_row)

        # --- Main area: left options, right plot ---
        center_row = QHBoxLayout()
        center_row.setSpacing(10)

        # Left: options + plot button
        left_layout = QVBoxLayout()
        left_layout.setSpacing(10)

        # Plot type radio buttons
        plot_type_label = QLabel("Plot type:")
        plot_type_label.setStyleSheet("font-weight: bold;")

        left_layout.addWidget(plot_type_label)

        self.radio_group = QButtonGroup(self)

        self.radio_screening = QRadioButton("plot the screening function")
        self.radio_inverse = QRadioButton("plot inverse dielectric matrix")

        self.radio_group.addButton(self.radio_screening, 0)
        self.radio_group.addButton(self.radio_inverse, 1)

        # Default: screening function
        self.radio_screening.setChecked(True)

        left_layout.addWidget(self.radio_screening)
        left_layout.addWidget(self.radio_inverse)

        # Plot button
        self.plot_button = QPushButton("Plot")
        self.plot_button.clicked.connect(self.plot)
        left_layout.addWidget(self.plot_button)

        left_layout.addStretch(1)

        # Right: Matplotlib plot
        self.plot_widget = MatplotlibWidget(self)
        self.plot_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        center_row.addLayout(left_layout, stretch=0)
        center_row.addWidget(self.plot_widget, stretch=1)

        main_layout.addLayout(center_row)

        self.setLayout(main_layout)

    # -------------------------
    # Helpers
    # -------------------------

    def _set_status(self, text: str) -> None:
        """Update the status label and process events so it appears immediately."""
        from PySide6.QtWidgets import QApplication

        self.status_label.setText(text)
        QApplication.processEvents()

    def _handle_back_clicked(self) -> None:
        if self._on_back is not None:
            self._on_back()

    def _browse_file(self, line_edit: QLineEdit) -> None:
        """Open a file dialog to select an HDF5 file and put the path into line_edit."""
        dlg = QFileDialog(self, "Select HDF5 file")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setNameFilter("HDF5 files (*.h5 *.hdf5);;All files (*.*)")
        if dlg.exec():
            files = dlg.selectedFiles()
            if files:
                line_edit.setText(files[0])
                # Force reload next time
                self._epshead = None
                self._qlen = None
                self._current_files = None

    # -------------------------
    # Data loading & plotting
    # -------------------------

    def _load_eps_head(self, force: bool = False) -> bool:
        """Load epshead and qlen from both eps0mat.h5 and epsmat.h5. Returns True on success."""
        eps0_path = self.eps0_edit.text().strip()
        eps_path = self.eps_edit.text().strip()

        if not eps0_path or not eps_path:
            self._set_status("Both eps0mat.h5 and epsmat.h5 are required")
            return False

        paths_tuple = (eps0_path, eps_path)

        if (
            not force
            and self._current_files == paths_tuple
            and self._epshead is not None
            and self._qlen is not None
        ):
            self._set_status("Data loaded")
            return True

        try:
            self._set_status("reading eps0mat.h5...")
            epshead0, qlen0 = get_eps_head(eps0_path)

            self._set_status("reading epsmat.h5...")
            epshead1, qlen1 = get_eps_head(eps_path)

        except Exception as exc:  # noqa: BLE001
            self._epshead = None
            self._qlen = None
            self._current_files = None
            self._set_status(f"Error: {exc}")
            print("[DielectricFunctionPage] Error reading eps files:", exc)
            return False

        try:
            # Concatenate along first axis (adjust if your functions return different shapes)
            self._epshead = np.concatenate([epshead0, epshead1], axis=0)
            self._qlen = np.concatenate([qlen0, qlen1], axis=0)
        except Exception as exc:  # noqa: BLE001
            self._epshead = None
            self._qlen = None
            self._current_files = None
            self._set_status(f"Concat error: {exc}")
            print("[DielectricFunctionPage] Error concatenating epshead/qlen:", exc)
            return False

        self._current_files = paths_tuple
        self._set_status("Data loaded")
        return True

    def plot(self) -> None:
        """Main entry point for plotting eps head."""
        if not self._load_eps_head(force=False):
            return

        if self._epshead is None or self._qlen is None:
            self._set_status("No data")
            return

        plot_inverse = self.radio_inverse.isChecked()

        ax = self.plot_widget.canvas.axes
        ax.cla()

        try:
            self._set_status("plotting...")
            plot_eps_head(ax, self._qlen, self._epshead, plot_inverse)
            self.plot_widget.canvas.draw()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[DielectricFunctionPage] Error while plotting eps head:", exc)
            return

        self._set_status("Done")
