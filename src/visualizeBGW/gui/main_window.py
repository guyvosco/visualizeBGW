"""
visualizeBGW.gui.main_window
PySide6 main window + menu screen for the visualizeBGW GUI.
"""

from __future__ import annotations

from typing import Callable, Dict, Optional

from PySide6.QtCore import Qt
from PySide6.QtWidgets import (
    QApplication,
    QHBoxLayout,
    QLabel,
    QMainWindow,
    QPushButton,
    QSpacerItem,
    QSizePolicy,
    QStackedWidget,
    QVBoxLayout,
    QWidget,
)
from .pages.crystal_structure_page import CrystalStructurePage
from .pages.wavefunction_projection_page import WavefunctionProjectionPage
from .pages.chi_converge_page import ChiConvergencePage
from .pages.ch_converge_page import CHConvergencePage
from .pages.bandstructure_page import BandstructurePage
from .pages.dielectric_page import DielectricFunctionPage
from .pages.absorption_spectra_page import AbsorptionSpectraPage
from .pages.exciton_components_page import ExcitonComponentsPage

# -------------------------
# Menu page
# -------------------------

class MenuPage(QWidget):
    """
    Main menu page with title/description/author and three rows:
        Structure: [Crystal] [Wavefunction projection] [Exciton projection]
        GW:        [chi converge] [CH converge] [Bandstructure] [Dielectric function]
        BSE:       [Absorption spectra] [Exciton components] [Bands contribution]

    When a button is clicked, it calls the `on_navigate(page_name)` callback.
    """

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        on_navigate: Optional[Callable[[str], None]] = None,
    ) -> None:
        super().__init__(parent)
        self._on_navigate = on_navigate
        self._build_ui()

    # ---- UI helpers ----

    def _build_ui(self) -> None:
        main_layout = QVBoxLayout()
        main_layout.setSpacing(16)
        main_layout.setContentsMargins(24, 24, 24, 24)

        # Top: Title / description / author
        header_layout = QVBoxLayout()
        header_layout.setSpacing(4)

        title_label = QLabel("visualizeBGW")
        title_label.setAlignment(Qt.AlignHCenter)
        title_label.setStyleSheet("font-size: 22px; font-weight: bold;")

        description_label = QLabel(
            "A Python toolkit and GUI for exploring BerkeleyGW outputs"
        )
        description_label.setAlignment(Qt.AlignHCenter)
        description_label.setWordWrap(True)

        author_label = QLabel("Author: Guy Vosco")
        author_label.setAlignment(Qt.AlignHCenter)

        header_layout.addWidget(title_label)
        header_layout.addWidget(description_label)
        header_layout.addWidget(author_label)

        main_layout.addLayout(header_layout)

        # Small vertical spacer between header and button rows
        main_layout.addSpacerItem(QSpacerItem(0, 20, QSizePolicy.Minimum, QSizePolicy.Minimum))

        # Rows: Structure, GW, BSE
        main_layout.addLayout(
            self._make_section_row(
                "Structure",
                [
                    ("Crystal structure", "structure_crystal"),
                    ("Wavefunction projection", "structure_wavefunction"),
                ],
            )
        )

        main_layout.addLayout(
            self._make_section_row(
                "GW",
                [
                    ("chi converge", "gw_chi"),
                    ("CH converge", "gw_ch"),
                    ("Bandstructure", "gw_bandstructure"),
                    ("Dielectric function", "gw_dielectric"),
                ],
            )
        )

        main_layout.addLayout(
            self._make_section_row(
                "BSE",
                [
                    ("Absorption spectra", "bse_absorption"),
                    ("Exciton components", "bse_components"),
                ],
            )
        )

        # main_layout.addLayout(
        #     self._make_section_row(
        #         "More features",
        #         [
        #             ("Cool feature", "cool_feature"),
        #         ],
        #     )
        # )

        # Add stretch at the bottom so content stays near the top
        main_layout.addStretch(1)

        self.setLayout(main_layout)

    def _make_section_row(
        self,
        label_text: str,
        buttons: list[tuple[str, str]],
    ) -> QHBoxLayout:
        """
        Create a horizontal row with:
            [Label]   [Button1] [Button2] [Button3] ...
        `buttons` is a list of (button_text, page_name) pairs.
        """
        row_layout = QHBoxLayout()
        row_layout.setSpacing(12)

        label = QLabel(label_text)
        label.setStyleSheet("font-weight: bold;")
        label.setFixedWidth(90)  # keeps labels nicely aligned
        row_layout.addWidget(label)

        # Buttons in the same row
        for text, page_name in buttons:
            btn = QPushButton(text)
            btn.setSizePolicy(QSizePolicy.Expanding, QSizePolicy.Fixed)
            btn.clicked.connect(lambda checked=False, name=page_name: self._navigate(name))
            row_layout.addWidget(btn)

        # Spacer so buttons don't stretch to full width if not needed
        row_layout.addStretch(1)

        return row_layout

    # ---- Navigation ----

    def _navigate(self, page_name: str) -> None:
        """Call the navigation callback if set."""
        if self._on_navigate is not None:
            self._on_navigate(page_name)


# -------------------------
# Main window
# -------------------------

class MainWindow(QMainWindow):
    """
    Main application window.

    - Uses a QStackedWidget as central widget.
    - Index 0 is the MenuPage.
    - Other indices are placeholders; you can later replace them with real pages.
    """

    def __init__(self, parent: Optional[QWidget] = None) -> None:
        super().__init__(parent)
        self.setWindowTitle("visualizeBGW")
        self.resize(900, 500)

        self._stack = QStackedWidget(self)
        self.setCentralWidget(self._stack)

        # Keep mapping from page name to index in the stacked widget
        self._page_indices: Dict[str, int] = {}

        self._create_pages()

    def _create_pages(self) -> None:
        # Page 0: menu
        menu_page = MenuPage(on_navigate=self.show_page)
        menu_index = self._stack.addWidget(menu_page)
        self._page_indices["menu"] = menu_index

        crystal_page = CrystalStructurePage(on_back=lambda: self.show_page("menu"))
        idx = self._stack.addWidget(crystal_page)
        self._page_indices["structure_crystal"] = idx
        wf_page = WavefunctionProjectionPage(on_back=lambda: self.show_page("menu"))
        idx = self._stack.addWidget(wf_page)
        self._page_indices["structure_wavefunction"] = idx
        
        chi_page = ChiConvergencePage(on_back=lambda: self.show_page("menu"))
        chi_index = self._stack.addWidget(chi_page)
        self._page_indices["gw_chi"] = chi_index
        ch_page = CHConvergencePage(on_back=lambda: self.show_page("menu"))
        ch_index = self._stack.addWidget(ch_page)
        self._page_indices["gw_ch"] = ch_index
        band_page = BandstructurePage(on_back=lambda: self.show_page("menu"))
        idx = self._stack.addWidget(band_page)
        self._page_indices["gw_bandstructure"] = idx
        dielectric_page = DielectricFunctionPage(on_back=lambda: self.show_page("menu"))
        dielectric_index = self._stack.addWidget(dielectric_page)
        self._page_indices["gw_dielectric"] = dielectric_index

        abs_page = AbsorptionSpectraPage(on_back=lambda: self.show_page("menu"))
        idx = self._stack.addWidget(abs_page)
        self._page_indices["bse_absorption"] = idx
        excit_comp_page = ExcitonComponentsPage(on_back=lambda: self.show_page("menu"))
        idx = self._stack.addWidget(excit_comp_page)
        self._page_indices["bse_components"] = idx

        # self._add_placeholder_page("cool_feature", "Cool feature (under development)")

        # Start with the menu
        self.show_page("menu")

    def _add_placeholder_page(self, page_name: str, text: str) -> None:
        """
        Creates a very simple placeholder page with a centered label.
        This is where you will later plug in your real plotting pages.
        """
        page = QWidget()
        layout = QVBoxLayout()
        label = QLabel(text)
        label.setAlignment(Qt.AlignCenter)
        layout.addWidget(label)

        # Simple "Back to menu" button
        back_btn = QPushButton("Back to menu")
        back_btn.setFixedWidth(140)
        back_btn.clicked.connect(lambda checked=False: self.show_page("menu"))
        layout.addWidget(back_btn, alignment=Qt.AlignHCenter | Qt.AlignBottom)

        page.setLayout(layout)

        index = self._stack.addWidget(page)
        self._page_indices[page_name] = index

    # ---- Public API ----

    def show_page(self, page_name: str) -> None:
        """
        Switch to a given page by name.
        The menu page is called "menu".
        """
        index = self._page_indices.get(page_name)
        if index is None:
            # If unknown, fall back to menu
            print(f"[MainWindow] Unknown page '{page_name}', falling back to menu.")
            index = self._page_indices["menu"]
        self._stack.setCurrentIndex(index)
