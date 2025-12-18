"""
BandstructurePage: GUI page for plotting bandstructures from BerkeleyGW eqp files.

Layout:
- Top-left: "Back to menu" button (calls on_back() callback if provided).
- Top-center: title label "Bandstructure".
- Next line: data file selector (QLineEdit + "Browse..." button) and,
  on the right, a status label (e.g. "reading...", "plotting...", "Done", etc.).

Options:
- First line: checkbox "mean-field energies".
- Second line: checkbox "quasiparticle corrections".
- Third line: open-text box for nv.

High-symmetry path:
- Title label "High-symmetry path:".
- Below: multiple rows, each with
      Label [str] ; k = ([kx], [ky], [kz])
  and from the second row onward, a "-" button on the left to delete the row.
- A "+" button below all rows to add another high-symmetry point row.

Data & plotting:
- read_eqp(eqp_file) -> kpts, emf, eqp
- plot_bandstructure(ax, kpts, emf, eqp, hsp, nv, plot_mf, plot_qp)
  where:
    - hsp is {label: (kx, ky, kz)}
    - nv is an integer (parsed from the nv text box)
    - plot_mf, plot_qp are bools from the checkboxes
"""

from __future__ import annotations

from typing import Callable, Dict, List, Optional, Tuple

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QCheckBox,
    QFileDialog,
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
from ...io.berkeleygw_readers import read_eqp
from ...plotting.berkeleygw_plots import plot_bandstructure


class BandstructurePage(QWidget):
    """
    Page for bandstructure plots.

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

        # Cached data from read_eqp
        self._kpts = None
        self._emf = None
        self._eqp = None
        self._current_file: Optional[str] = None

        # High-symmetry path rows
        self._hsp_rows: List[Dict[str, QWidget]] = []
        self._hsp_rows_layout: Optional[QVBoxLayout] = None

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

        title_label = QLabel("Bandstructure")
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

        file_label = QLabel("eqp.dat file:")

        self.file_edit = QLineEdit()
        self.file_edit.setPlaceholderText("Select Eqp data file...")
        self.file_edit.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
        self.file_edit.setMinimumWidth(450)

        browse_button = QPushButton("Browse…")
        browse_button.clicked.connect(self._browse_file)

        self.status_label = QLabel("")
        self.status_label.setAlignment(Qt.AlignRight | Qt.AlignVCenter)
        self.status_label.setMinimumWidth(140)

        file_row.addWidget(file_label)
        file_row.addWidget(self.file_edit)
        file_row.addWidget(browse_button)
        file_row.addStretch(1)
        file_row.addWidget(self.status_label)

        main_layout.addLayout(file_row)

        # --- Main area: left options + HSP, right plot ---
        center_row = QHBoxLayout()
        center_row.setSpacing(10)

        left_panel = QVBoxLayout()
        left_panel.setSpacing(10)

        # Options (mean-field, qp, nv)
        options_form = QFormLayout()
        options_form.setLabelAlignment(Qt.AlignRight)
        options_form.setFormAlignment(Qt.AlignTop)

        self.meanfield_checkbox = QCheckBox("mean-field energies")
        self.meanfield_checkbox.setChecked(True)

        self.qp_checkbox = QCheckBox("quasiparticle corrections")
        self.qp_checkbox.setChecked(True)

        self.nv_edit = QLineEdit()
        self.nv_edit.setPlaceholderText("nv")

        options_form.addRow("", self.meanfield_checkbox)
        options_form.addRow("", self.qp_checkbox)
        options_form.addRow("nv:", self.nv_edit)

        left_panel.addLayout(options_form)

        # High-symmetry path section
        hsp_title = QLabel("High-symmetry path:")
        hsp_title.setStyleSheet("font-weight: bold;")
        left_panel.addWidget(hsp_title)

        # Container for HSP rows
        self._hsp_rows_layout = QVBoxLayout()
        self._hsp_rows_layout.setSpacing(4)

        # Initialize with 2 rows by default
        for _ in range(2):
            self._add_hsp_row()

        left_panel.addLayout(self._hsp_rows_layout)

        # "+" button to add more HSP points
        self.add_hsp_button = QPushButton("+")
        self.add_hsp_button.setFixedWidth(24)
        self.add_hsp_button.clicked.connect(self._add_hsp_row)
        left_panel.addWidget(self.add_hsp_button, alignment=Qt.AlignLeft)

        # Plot button
        self.plot_button = QPushButton("Plot")
        self.plot_button.clicked.connect(self.plot)
        left_panel.addWidget(self.plot_button)

        left_panel.addStretch(1)

        # Right: Matplotlib plot
        self.plot_widget = MatplotlibWidget(self)
        self.plot_widget.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Expanding)

        center_row.addLayout(left_panel, stretch=0)
        center_row.addWidget(self.plot_widget, stretch=1)

        main_layout.addLayout(center_row)

        self.setLayout(main_layout)

        # Ensure minus-button visibility is consistent
        self._update_hsp_minus_buttons()

    # -------------------------
    # High-symmetry path rows
    # -------------------------

    def _add_hsp_row(self) -> None:
        """Add a new HSP row: Label [str] ; k = ([kx], [ky], [kz])."""
        if self._hsp_rows_layout is None:
            return

        row_layout = QHBoxLayout()
        row_layout.setSpacing(4)

        label_edit = QLineEdit()
        label_edit.setPlaceholderText("Label")
        label_edit.setFixedWidth(75)

        lbl_k_equals = QLabel(" ; k = (")
        kx_edit = QLineEdit()
        kx_edit.setPlaceholderText("kx")
        kx_edit.setFixedWidth(100)
        ky_edit = QLineEdit()
        ky_edit.setPlaceholderText("ky")
        ky_edit.setFixedWidth(100)
        kz_edit = QLineEdit()
        kz_edit.setPlaceholderText("kz")
        kz_edit.setFixedWidth(100)
        lbl_close = QLabel(")")

        minus_button = QPushButton("-")
        minus_button.setFixedWidth(24)
        minus_button.clicked.connect(
            lambda checked=False, btn=minus_button: self._remove_hsp_row(btn)
        )

        row_layout.addWidget(label_edit)
        row_layout.addWidget(lbl_k_equals)
        row_layout.addWidget(kx_edit)
        row_layout.addWidget(ky_edit)
        row_layout.addWidget(kz_edit)
        row_layout.addWidget(lbl_close)
        row_layout.addWidget(minus_button)
        row_layout.addStretch(1)

        # Wrap layout in a QWidget so it's easier to remove
        row_widget = QWidget()
        row_widget.setLayout(row_layout)

        self._hsp_rows_layout.addWidget(row_widget)

        self._hsp_rows.append(
            {
                "container": row_widget,
                "label": label_edit,
                "kx": kx_edit,
                "ky": ky_edit,
                "kz": kz_edit,
                "minus": minus_button,
            }
        )

        self._update_hsp_minus_buttons()

    def _remove_hsp_row(self, minus_button: QPushButton) -> None:
        """Remove the row associated with the given minus button."""
        for i, row in enumerate(self._hsp_rows):
            if row["minus"] is minus_button:
                container: QWidget = row["container"]
                # Remove from layout and delete
                if self._hsp_rows_layout is not None:
                    self._hsp_rows_layout.removeWidget(container)
                container.deleteLater()
                self._hsp_rows.pop(i)
                break

        self._update_hsp_minus_buttons()

    def _update_hsp_minus_buttons(self) -> None:
        """Show '-' only from second row onwards."""
        for i, row in enumerate(self._hsp_rows):
            minus: QPushButton = row["minus"]
            minus.setVisible(i > 0)

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
        """Open a file dialog to select the eqp data file."""
        dlg = QFileDialog(self, "Select Eqp data file")
        dlg.setFileMode(QFileDialog.ExistingFile)
        dlg.setNameFilter("All files (*.*)")
        if dlg.exec():
            files = dlg.selectedFiles()
            if files:
                self.file_edit.setText(files[0])
                # Force reload next time we plot
                self._kpts = self._emf = self._eqp = None
                self._current_file = None

    # -------------------------
    # Data loading & plotting
    # -------------------------

    def _load_eqp(self, force: bool = False) -> bool:
        """Load kpts, emf, eqp from the eqp file. Returns True on success."""
        file_path = self.file_edit.text().strip()
        if not file_path:
            self._set_status("No file selected")
            return False

        if (
            not force
            and self._current_file == file_path
            and self._kpts is not None
            and self._emf is not None
            and self._eqp is not None
        ):
            self._set_status("Data loaded")
            return True

        try:
            self._set_status("reading...")
            kpts, emf, eqp = read_eqp(file_path)
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Error: {exc}")
            print("[BandstructurePage] Error reading eqp file:", exc)
            self._kpts = self._emf = self._eqp = None
            self._current_file = None
            return False

        self._kpts, self._emf, self._eqp = kpts, emf, eqp
        self._current_file = file_path
        self._set_status("Data loaded")
        return True

    def _collect_hsp(self) -> Optional[List[str, Tuple[float, float, float]]]:
        """
        Build the high-symmetry path List from the UI.

        Returns
        -------
        list or None
            [(label, (kx, ky, kz))] on success, or None if invalid input encountered.
        """
        hsp: List[Tuple[str, Tuple[float, float, float]]] = []

        for row in self._hsp_rows:
            label_edit: QLineEdit = row["label"]
            kx_edit: QLineEdit = row["kx"]
            ky_edit: QLineEdit = row["ky"]
            kz_edit: QLineEdit = row["kz"]

            label = label_edit.text().strip()
            kx_str = kx_edit.text().strip()
            ky_str = ky_edit.text().strip()
            kz_str = kz_edit.text().strip()

            # Skip completely empty rows
            if not label and not (kx_str or ky_str or kz_str):
                continue

            if not label:
                self._set_status("HSP error: missing label")
                return None

            try:
                kx = float(kx_str)
                ky = float(ky_str)
                kz = float(kz_str)
            except ValueError:
                self._set_status("HSP error: invalid k components")
                return None

            hsp.append((label, (kx, ky, kz)))

        if not hsp:
            self._set_status("No HSP points defined")
            return None

        return hsp

    def plot(self) -> None:
        """Main entry point for plotting bandstructure."""
        # Load data
        if not self._load_eqp(force=False):
            return

        if self._kpts is None or self._emf is None or self._eqp is None:
            self._set_status("No data")
            return

        # Get nv
        nv_str = self.nv_edit.text().strip()
        if not nv_str:
            self._set_status("nv is required")
            return
        try:
            nv = int(nv_str)
        except ValueError:
            self._set_status("nv must be an integer")
            return

        # Collect HSP dictionary
        hsp = self._collect_hsp()
        if hsp is None:
            return

        plot_mf = self.meanfield_checkbox.isChecked()
        plot_qp = self.qp_checkbox.isChecked()

        ax = self.plot_widget.canvas.axes
        ax.cla()

        try:
            self._set_status("plotting...")
            plot_bandstructure(
                ax,
                self._kpts,
                self._emf,
                self._eqp,
                hsp,
                nv,
                plot_mf,
                plot_qp,
            )
            self.plot_widget.canvas.draw()
        except Exception as exc:  # noqa: BLE001
            self._set_status(f"Plot error: {exc}")
            print("[BandstructurePage] Error while plotting bandstructure:", exc)
            return

        self._set_status("Done")
