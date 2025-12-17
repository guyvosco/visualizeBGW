"""
CHConvergencePage: GUI page for plotting GW CH convergence.

Layout:
- Top-left: "Back to menu" button (calls on_back() callback if provided).
- Top-center: title label "CH converge".
- Next line: data file selector (QLineEdit + "Browse..." button) and,
  on the right, a status label (e.g. "reading...", "plotting...", "Done", etc.).
- Main area: left panel with options, right panel with Matplotlib plot.

Options on the left:
- k-point: QComboBox populated from ch_converge_data keys (kpts).
- quantity: QComboBox with ["vbm", "cbm", "diff"].

Plotting:
- read_ch_converge(ch_converge_file) → ch_converge_data dict-like, keys are k-points.
- plot_ch_convergence(ax, ch_converge_data[kpt], quantity).
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
from ...io.berkeleygw_readers import read_ch_converge
from ...plotting.berkeleygw_plots import plot_ch_convergence


class CHConvergencePage(QWidget):
    """
    Page for CH-convergence plots.

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
        self._ch_converge_data: Optional[Dict] = None
        self._current_file: Optional[str] = None
        self._kpt_keys: List = []  # ordered list of dict keys corresponding to combo entries

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

        title_label = QLabel("CH converge")
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

        file_label = QLabel("ch_converge.dat file:")

        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("Select CH-convergence data file...")
        self.file_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        # Make it comfortably wide
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

        # k-point ComboBox (filled after reading data)
        self.kpt_combo = QComboBox()
        self.kpt_combo.setEnabled(False)  # until data is loaded

        # quantity ComboBox
        self.quantity_combo = QComboBox()
        self.quantity_combo.addItems(["vbm", "cbm", "diff"])

        form.addRow("k-point:", self.kpt_combo)
        form.addRow("quantity:", self.quantity_combo)

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
        """Open a file dialog to select the CH-convergence data file."""
        dlg = QFileDialog(self, "Select CH-convergence data file")
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
        Load ch_converge_data from the file in file_edit, if needed.

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
            and self._ch_converge_data is not None
            and self._current_file == file_path
        ):
            # Already loaded from this file
            self._set_status("Data loaded")
            return

        try:
            self._set_status("reading...")
            # Call user-provided reader
            ch_data = read_ch_converge(file_path)
        except Exception as exc:  # noqa: BLE001
            self._ch_converge_data = None
            self._current_file = None
            self._kpt_keys = []
            self.kpt_combo.clear()
            self.kpt_combo.setEnabled(False)
            self.plot_button.setEnabled(False)
            self._set_status(f"Error: {exc}")
            print("[CHConvergencePage] Error reading file:", exc)
            return

        # Store data and update UI
        self._ch_converge_data = ch_data
        self._current_file = file_path

        # Populate k-point combo box using an explicit list of keys.
        self.kpt_combo.clear()
        self._kpt_keys = list(ch_data.keys())

        try:
            # Try to sort keys for a nicer order if possible
            self._kpt_keys.sort()
        except Exception:
            # If keys are not comparable, just keep original order
            pass

        for key in self._kpt_keys:
            # Use string representation as label; key itself stays in _kpt_keys
            self.kpt_combo.addItem(str(key))

        has_kpts = len(self._kpt_keys) > 0
        self.kpt_combo.setEnabled(has_kpts)
        self.plot_button.setEnabled(has_kpts)

        if has_kpts:
            self._set_status("Data loaded")
        else:
            self._set_status("No k-points found")

    def plot(self) -> None:
        """Main entry point for plotting CH convergence."""
        # Ensure data is loaded
        self.load_data(force=False)
        if self._ch_converge_data is None or not self._kpt_keys:
            return

        idx = self.kpt_combo.currentIndex()
        if idx < 0 or idx >= len(self._kpt_keys):
            self._set_status("No k-point selected")
            return

        # Use the key from our separate list so we never rely on userData
        kpt_key = self._kpt_keys[idx]

        quantity = self.quantity_combo.currentText()

        # Extract ch_converge_data for this k-point
        try:
            ch_kpt_data = self._ch_converge_data[kpt_key]
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Invalid k-point: {exc}")
            print("[CHConvergencePage] Error accessing k-point data:", exc)
            return

        # Clear axes before plotting
        ax = self.plot_widget.canvas.axes
        ax.cla()

        # Plot
        try:
            self._set_status("plotting...")
            plot_ch_convergence(ax, ch_kpt_data, quantity)
            self.plot_widget.canvas.draw()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[CHConvergencePage] Error while plotting:", exc)
            return

        self._set_status("Done")
