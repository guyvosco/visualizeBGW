"""
AbsorptionSpectraPage: GUI page for plotting BSE absorption spectra.

Workflow
--------
- Read main eigenvalues file (eigenvalues.dat) with:
    eigs, dipole = read_eigenvalues(eigenvalues_file)

- Optionally read a second file (e.g. "no e-h" case) with:
    eigs_noeh, dipole_noeh = read_eigenvalues_noeh(eigenvalues_noeh_file)

- Options:
    * broadening (float)
    * broadening function: ["Gaussian", "Lorentzian", "Voigt"]
    * x-axis limits: xmin, xmax (floats)

- Build spectra with:
    energy, spectra = get_absorption_spectra(eigs, dipole, broadening, function)

- Optional no-e–h spectrum:
    energy_noeh, spectra_noeh = get_absorption_spectra(eigs_noeh, dipole_noeh, broadening, function)

- Plot with:
    plot_absorption_spectra(ax, energy, spectra, xmin, xmax, noeh_spectra=None or spectra_noeh)
"""

from __future__ import annotations

from typing import Callable, Optional, Tuple

import numpy as np

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QFileDialog,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QComboBox,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from ..matplotlib_widget import MatplotlibWidget
from ...io.berkeleygw_readers import read_eigenvalues, read_eigenvalues_noeh
from ...analysis.berkeleygw_processing import get_absorption_spectra
from ...plotting.berkeleygw_plots import plot_absorption_spectra


class AbsorptionSpectraPage(QWidget):
    """
    Page for absorption spectra plots.

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

        # Cached main data
        self._eigs: Optional[np.ndarray] = None
        self._dipole: Optional[np.ndarray] = None
        self._current_file: Optional[str] = None

        # Cached optional no-e–h data
        self._eigs_noeh: Optional[np.ndarray] = None
        self._dipole_noeh: Optional[np.ndarray] = None
        self._current_file_noeh: Optional[str] = None

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

        title_label = QLabel("Absorption spectra")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 18px; font-weight: bold;")

        top_row.addWidget(back_button, alignment=Qt.AlignLeft)
        top_row.addStretch(1)
        top_row.addWidget(title_label, stretch=0, alignment=Qt.AlignCenter)
        top_row.addStretch(1)

        main_layout.addLayout(top_row)

        # --- File selectors + status ---

        # Main eigenvalues.dat
        main_file_row = QHBoxLayout()
        main_file_row.setSpacing(8)

        main_file_label = QLabel("eigenvalues.dat file:")
        self.main_file_edit = QLineEdit()
        self.main_file_edit.setPlaceholderText("Select eigenvalues.dat …")
        self.main_file_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.main_file_edit.setMinimumWidth(450)

        main_browse = QPushButton("Browse…")
        main_browse.clicked.connect(
            lambda: self._browse_file(self.main_file_edit, main=True)
        )

        main_file_row.addWidget(main_file_label)
        main_file_row.addWidget(self.main_file_edit)
        main_file_row.addWidget(main_browse)
        main_file_row.addStretch(1)

        main_layout.addLayout(main_file_row)

        # Toggle for optional second file
        toggle_row = QHBoxLayout()
        toggle_row.setSpacing(8)

        self.noeh_toggle = QCheckBox("Add optional no e–h (noeh) spectrum")
        self.noeh_toggle.stateChanged.connect(self._update_noeh_visibility)

        toggle_row.addWidget(self.noeh_toggle)
        toggle_row.addStretch(1)

        main_layout.addLayout(toggle_row)

        # Optional second eigenvalues file + status
        noeh_file_row = QHBoxLayout()
        noeh_file_row.setSpacing(8)

        noeh_file_label = QLabel("eigenvalues_noeh.dat file:")
        self.noeh_file_edit = QLineEdit()
        self.noeh_file_edit.setPlaceholderText("Select eigenvalues_noeh.dat …")
        self.noeh_file_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.noeh_file_edit.setMinimumWidth(450)

        noeh_browse = QPushButton("Browse…")
        noeh_browse.clicked.connect(
            lambda: self._browse_file(self.noeh_file_edit, main=False)
        )

        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.status_label.setMinimumWidth(180)

        noeh_file_row.addWidget(noeh_file_label)
        noeh_file_row.addWidget(self.noeh_file_edit)
        noeh_file_row.addWidget(noeh_browse)
        noeh_file_row.addStretch(1)
        noeh_file_row.addWidget(self.status_label)

        main_layout.addLayout(noeh_file_row)

        # We'll show/hide this whole row when toggling
        self._noeh_file_row_layout = noeh_file_row
        self._noeh_widgets = [noeh_file_label, self.noeh_file_edit, noeh_browse]
        self._update_noeh_visibility()

        # --- Main area: left options, right plot ---
        center_row = QHBoxLayout()
        center_row.setSpacing(10)

        # Left panel: options + plot button
        left_layout = QVBoxLayout()
        left_layout.setSpacing(10)

        # Broadening options
        broad_row = QHBoxLayout()
        broad_row.setSpacing(6)

        broad_label = QLabel("broadening:")
        self.broad_edit = QLineEdit()
        self.broad_edit.setPlaceholderText("e.g. 0.05")
        self.broad_edit.setMaximumWidth(100)

        broad_row.addWidget(broad_label)
        broad_row.addWidget(self.broad_edit)
        broad_row.addStretch(1)

        left_layout.addLayout(broad_row)

        # Broadening function
        func_row = QHBoxLayout()
        func_row.setSpacing(6)

        func_label = QLabel("broadening function:")
        self.func_combo = QComboBox()
        self.func_combo.addItems(["gaussian", "lorentzian", "voigt"])
        self.func_combo.setMaximumWidth(140)

        func_row.addWidget(func_label)
        func_row.addWidget(self.func_combo)
        func_row.addStretch(1)

        left_layout.addLayout(func_row)

        # x-axis limits
        xlim_row = QHBoxLayout()
        xlim_row.setSpacing(6)

        xmin_label = QLabel("xmin:")
        self.xmin_edit = QLineEdit()
        self.xmin_edit.setPlaceholderText("min E")
        self.xmin_edit.setMaximumWidth(100)

        xmax_label = QLabel("xmax:")
        self.xmax_edit = QLineEdit()
        self.xmax_edit.setPlaceholderText("max E")
        self.xmax_edit.setMaximumWidth(100)

        xlim_row.addWidget(xmin_label)
        xlim_row.addWidget(self.xmin_edit)
        xlim_row.addSpacing(10)
        xlim_row.addWidget(xmax_label)
        xlim_row.addWidget(self.xmax_edit)
        xlim_row.addStretch(1)

        left_layout.addLayout(xlim_row)

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

    def _update_noeh_visibility(self) -> None:
        """Show or hide the second file selector based on the toggle."""
        visible = self.noeh_toggle.isChecked()
        for w in self._noeh_widgets:
            w.setVisible(visible)

    def _browse_file(self, line_edit: QLineEdit, main: bool) -> None:
        """Open a file dialog to select an eigenvalues file."""
        dlg = QFileDialog(self, "Select eigenvalues file")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setNameFilter("All files (*.*)")
        if dlg.exec():
            files = dlg.selectedFiles()
            if files:
                line_edit.setText(files[0])
                # Force reload on next plot
                if main:
                    self._eigs = self._dipole = None
                    self._current_file = None
                else:
                    self._eigs_noeh = self._dipole_noeh = None
                    self._current_file_noeh = None

    # -------------------------
    # Data loading
    # -------------------------

    def _load_main_eigenvalues(self, force: bool = False) -> bool:
        file_path = self.main_file_edit.text().strip()
        if not file_path:
            self._set_status("Main eigenvalues file is required")
            return False

        if (
            not force
            and self._current_file == file_path
            and self._eigs is not None
            and self._dipole is not None
        ):
            self._set_status("Main data loaded")
            return True

        try:
            self._set_status("reading main eigenvalues...")
            eigs, dipole = read_eigenvalues(file_path)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Error: {exc}")
            print("[AbsorptionSpectraPage] Error reading main eigenvalues:", exc)
            self._eigs = self._dipole = None
            self._current_file = None
            return False

        self._eigs, self._dipole = eigs, dipole
        self._current_file = file_path
        self._set_status("Main data loaded")
        return True

    def _load_noeh_eigenvalues(self, force: bool = False) -> bool:
        """Load optional no-e–h eigenvalues if the toggle is on."""
        if not self.noeh_toggle.isChecked():
            # Not requested
            return False

        file_path = self.noeh_file_edit.text().strip()
        if not file_path:
            # Optional – if user toggled but didn't select, we just skip
            self._set_status("Noeh file not set; plotting main spectrum only")
            return False

        if (
            not force
            and self._current_file_noeh == file_path
            and self._eigs_noeh is not None
            and self._dipole_noeh is not None
        ):
            # Already loaded
            return True

        try:
            self._set_status("reading noeh eigenvalues...")
            eigs_noeh, dipole_noeh = read_eigenvalues_noeh(file_path)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Noeh error: {exc}")
            print("[AbsorptionSpectraPage] Error reading noeh eigenvalues:", exc)
            self._eigs_noeh = self._dipole_noeh = None
            self._current_file_noeh = None
            return False

        self._eigs_noeh, self._dipole_noeh = eigs_noeh, dipole_noeh
        self._current_file_noeh = file_path
        return True

    # -------------------------
    # Plotting
    # -------------------------

    def plot(self) -> None:
        """Main entry point for plotting absorption spectra."""
        self._set_status("Computing the absorption spectrum...")
        # Load main data
        if not self._load_main_eigenvalues(force=False):
            return

        if self._eigs is None or self._dipole is None:
            self._set_status("No main data")
            return

        # Parse broadening
        broad_str = self.broad_edit.text().strip()
        if not broad_str:
            self._set_status("Broadening is required")
            return
        try:
            broadening = float(broad_str)
        except ValueError:
            self._set_status("Broadening must be a float")
            return

        # Broadening function
        broad_func = self.func_combo.currentText()

        # x-axis limits (optional: we can allow empty to mean 'auto')
        xmin_str = self.xmin_edit.text().strip()
        xmax_str = self.xmax_edit.text().strip()
        try:
            xmin = float(xmin_str) if xmin_str else None
            xmax = float(xmax_str) if xmax_str else None
        except ValueError:
            self._set_status("xmin/xmax must be floats")
            return

        # Compute main spectrum
        try:
            # self._set_status("computing main spectrum...")
            energy, spectra = get_absorption_spectra(
                self._eigs, self._dipole, broadening, broad_func
            )
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Error in main spectrum: {exc}")
            print("[AbsorptionSpectraPage] Error computing main spectrum:", exc)
            return

        # Optional noeh spectrum
        noeh_energy = None
        noeh_spectra = None
        if self._load_noeh_eigenvalues(force=False):
            if self._eigs_noeh is not None and self._dipole_noeh is not None:
                try:
                    # self._set_status("computing noeh spectrum...")
                    energy_noeh, spectra_noeh = get_absorption_spectra(
                        self._eigs_noeh, self._dipole_noeh, broadening, broad_func
                    )
                    # Assume same energy grid; if not, user should handle resampling upstream
                    noeh_energy = energy_noeh
                    noeh_spectra = spectra_noeh
                except Exception as exc:  # noqa: BLE001
                    self._set_status(f"Noeh spectrum error: {exc}")
                    print("[AbsorptionSpectraPage] Error computing noeh spectrum:", exc)
                    noeh_energy = None
                    noeh_spectra = None

        # Plot
        ax = self.plot_widget.canvas.axes
        ax.cla()

        # For xmin/xmax: if None, let the plotting helper decide or apply later
        try:
            self._set_status("plotting...")
            plot_absorption_spectra(
                ax,
                energy,
                spectra,
                xmin,
                xmax,
                noeh_energy=noeh_energy,
                noeh_spectra=noeh_spectra,
            )
            self.plot_widget.canvas.draw()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[AbsorptionSpectraPage] Error while plotting absorption:", exc)
            return

        self._set_status("Done")
