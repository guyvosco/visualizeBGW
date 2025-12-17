"""
Application entry point for the visualizeBGW GUI.

This module exposes a `main()` function that is used by the console
script entry point defined in pyproject.toml:

    visualizeBGW-gui = "visualizeBGW.gui.main:main"
"""

from __future__ import annotations

import sys

from PySide6.QtWidgets import QApplication

from .main_window import MainWindow


def main() -> None:
    """
    Start the visualizeBGW GUI.

    This function is called by the console script entry point.
    """
    app = QApplication(sys.argv)

    window = MainWindow()
    window.show()

    # Start the Qt event loop
    sys.exit(app.exec())


if __name__ == "__main__":
    # Allow running `python -m visualizeBGW.gui.main` directly
    main()
