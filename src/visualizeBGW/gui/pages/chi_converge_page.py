"""
ChiConvergencePage: GUI page for plotting GW chi convergence.

Layout:
- Top-left: "Back to menu" button (calls on_back() callback if provided).
- Top-center: title label "chi converge".
- Next line: data file selector (QLineEdit + "Browse..." button) and,
  on the right, a status label (e.g. "reading...", "plotting...", "Done", etc.).
- Main area: left panel with options, right panel with Matplotlib plot.

Options on the left:
- q-point: QComboBox populated from chi_converge_data keys (qpts).
- G = G' = : ComboBox with ["0", "G_max"].
# - y-axis scale: ComboBox with ["linear", "log"].

Plotting:
- read_chi_converge(chi_converge_file) → chi_converge_data dict-like, keys are q-points.
- plot_chi_converge(ax, chi_converge_data[qpt], G_selected).
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QFileDialog,
    QComboBox,
    QFormLayout,
    QHBoxLayout,
    QLabel,
    QLineEdit,
    QPushButton,
    QSizePolicy,
    QVBoxLayout,
    QWidget,
)

from ..matplotlib_widget import MatplotlibWidget
from ...io.berkeleygw_readers import read_chi_converge
from ...plotting.berkeleygw_plots import plot_chi_converge


class ChiConvergencePage(QWidget):
    """
    Page for chi-convergence plots.

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

        # Data cache
        self._chi_converge_data: Optional[Dict] = None
        self._current_file: Optional[str] = None
        self._qpt_keys: List = (
            []
        )  # ordered list of dict keys corresponding to combo entries

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

        self.back_button = QPushButton("Back to menu")
        self.back_button.setFixedWidth(120)
        self.back_button.clicked.connect(self._handle_back_clicked)

        title_label = QLabel("chi converge")
        title_label.setAlignment(Qt.AlignCenter)
        title_label.setStyleSheet("font-size: 18px; font-weight: bold;")

        top_row.addWidget(self.back_button, alignment=Qt.AlignLeft)
        top_row.addStretch(1)
        top_row.addWidget(title_label, stretch=0, alignment=Qt.AlignCenter)
        top_row.addStretch(1)

        main_layout.addLayout(top_row)

        # --- Second row: file selector + status label ---
        file_row = QHBoxLayout()
        file_row.setSpacing(8)

        file_label = QLabel("chi_converge.dat file:")

        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("Select chi-convergence data file...")
        self.file_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        # Make it noticeably wider
        self.file_edit.setMinimumWidth(450)

        browse_button = QPushButton("Browse…")
        browse_button.clicked.connect(self._browse_file)

        # Status label on the right
        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.status_label.setMinimumWidth(140)

        file_row.addWidget(file_label)
        file_row.addWidget(self.file_edit)
        file_row.addWidget(browse_button)
        file_row.addStretch(1)
        file_row.addWidget(self.status_label)

        main_layout.addLayout(file_row)

        # --- Main area: left options, right plot ---
        center_row = QHBoxLayout()
        center_row.setSpacing(10)

        # Left: options
        options_layout = QVBoxLayout()
        options_layout.setSpacing(12)

        form = QFormLayout()
        form.setLabelAlignment(Qt.AlignRight)
        form.setFormAlignment(Qt.AlignTop)

        # q-point ComboBox (filled after reading data)
        self.qpt_combo = QComboBox()
        self.qpt_combo.setEnabled(False)  # until data is loaded

        # G = G' = ComboBox
        self.G_combo = QComboBox()
        # Display text vs internal value (user can adjust internal values as needed)
        self.G_combo.addItem("0", userData="0")
        self.G_combo.addItem("G_max", userData="Gmax")

        # y-axis scale ComboBox
        # self.scale_combo = QComboBox()
        # self.scale_combo.addItems(["linear", "log"])

        form.addRow("q-point:", self.qpt_combo)
        form.addRow("G = G' =", self.G_combo)
        # form.addRow("y-axis scale:", self.scale_combo)

        options_layout.addLayout(form)

        # Plot button
        self.plot_button = QPushButton("Plot")
        self.plot_button.setEnabled(False)  # until data is loaded
        self.plot_button.clicked.connect(self.plot)

        options_layout.addWidget(self.plot_button)
        options_layout.addStretch(1)

        # Right: Matplotlib plot
        self.plot_widget = MatplotlibWidget(self)
        self.plot_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        center_row.addLayout(options_layout, stretch=0)
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

    def _browse_file(self) -> None:
        """Open a file dialog to select the chi-convergence data file."""
        dlg = QFileDialog(self, "Select chi-convergence data file")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setNameFilter("All files (*.*)")
        if dlg.exec():
            files = dlg.selectedFiles()
            if files:
                self.file_edit.setText(files[0])
                # After selecting a new file, try reading data
                self.load_data(force=True)

    # -------------------------
    # Data loading & plotting
    # -------------------------

    def load_data(self, force: bool = False) -> None:
        """
        Load chi_converge_data from the file in file_edit, if needed.

        Parameters
        ----------
        force : bool
            If True, always reload the file, even if the path hasn't changed.
        """
        file_path = self.file_edit.text().strip()
        if not file_path:
            self._set_status("No file selected")
            return

        if (
            not force
            and self._chi_converge_data is not None
            and self._current_file == file_path
        ):
            # Already loaded from this file
            self._set_status("Data loaded")
            return

        try:
            self._set_status("reading...")
            # Call user-provided reader
            chi_data = read_chi_converge(file_path)
        except Exception as exc:  # noqa: BLE001
            self._chi_converge_data = None
            self._current_file = None
            self._qpt_keys = []
            self.qpt_combo.clear()
            self.qpt_combo.setEnabled(False)
            self.plot_button.setEnabled(False)
            self._set_status(f"Error: {exc}")
            print("[ChiConvergencePage] Error reading file:", exc)
            return

        # Store data and update UI
        self._chi_converge_data = chi_data
        self._current_file = file_path

        # Populate q-point combo box using an explicit list of keys.
        self.qpt_combo.clear()
        self._qpt_keys = list(chi_data.keys())

        try:
            # Try to sort keys for a nicer order if possible
            self._qpt_keys.sort()
        except Exception:
            # If keys are not comparable, just keep original order
            pass

        for key in self._qpt_keys:
            # Use string representation as label; key itself stays in _qpt_keys
            self.qpt_combo.addItem(str(key))

        has_qpts = len(self._qpt_keys) > 0
        self.qpt_combo.setEnabled(has_qpts)
        self.plot_button.setEnabled(has_qpts)

        if has_qpts:
            self._set_status("Data loaded")
        else:
            self._set_status("No q-points found")

    def plot(self) -> None:
        """Main entry point for plotting chi convergence."""
        # Ensure data is loaded
        self.load_data(force=False)
        if self._chi_converge_data is None or not self._qpt_keys:
            return

        idx = self.qpt_combo.currentIndex()
        if idx < 0 or idx >= len(self._qpt_keys):
            self._set_status("No q-point selected")
            return

        # Use the key from our separate list so we never rely on userData
        qpt_key = self._qpt_keys[idx]

        G_value = self.G_combo.currentData()  # internal representation
        # yscale = self.scale_combo.currentText()

        # Extract chi_converge_data for this q-point
        try:
            chi_q_data = self._chi_converge_data[qpt_key]
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Invalid q-point: {exc}")
            print("[ChiConvergencePage] Error accessing q-point data:", exc)
            return

        # Clear axes before plotting
        ax = self.plot_widget.canvas.axes
        ax.cla()

        # Plot
        try:
            self._set_status("plotting...")
            plot_chi_converge(ax, chi_q_data, G_value)
            # ax.set_yscale(yscale)
            self.plot_widget.canvas.draw()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[ChiConvergencePage] Error while plotting:", exc)
            return

        self._set_status("Done")
