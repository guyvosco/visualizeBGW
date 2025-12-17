"""
Matplotlib widget for the visualizeBGW GUI.

Provides:
- MatplotlibCanvas: a FigureCanvasQTAgg with a Figure and single Axes.
- MatplotlibWidget: a QWidget that contains a canvas + navigation toolbar.

Typical usage in a page:

    from .matplotlib_widget import MatplotlibWidget

    class MyPage(QWidget):
        def __init__(self, parent=None):
            super().__init__(parent)
            layout = QVBoxLayout(self)
            self.plot_widget = MatplotlibWidget(self)
            layout.addWidget(self.plot_widget)

            ax = self.plot_widget.canvas.axes
            ax.plot([0, 1, 2], [0, 1, 0])
            self.plot_widget.canvas.draw()

You can reuse MatplotlibWidget in all your plotting pages.
"""

from __future__ import annotations

from typing import Optional, Tuple

from matplotlib.backends.backend_qtagg import FigureCanvasQTAgg
from matplotlib.backends.backend_qtagg import NavigationToolbar2QT
from matplotlib.figure import Figure

from PySide6.QtWidgets import QVBoxLayout, QWidget


class MatplotlibCanvas(FigureCanvasQTAgg):
    """
    A basic Matplotlib canvas with a single Axes.

    Attributes
    ----------
    figure : matplotlib.figure.Figure
        The Matplotlib Figure instance used for plotting.
    axes : matplotlib.axes.Axes
        The single Axes contained in the Figure.
    """

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        width: float = 6.0,
        height: float = 4.0,
        dpi: int = 100,
    ) -> None:
        fig = Figure(figsize=(width, height), dpi=dpi)
        self.axes = fig.add_subplot(111)

        super().__init__(fig)
        self.setParent(parent)

        # Optional: enable tight layout by default for nicer spacing
        fig.tight_layout()

    def clear(self) -> None:
        """Clear the current axes."""
        self.axes.cla()
        self.draw()

    def plot_dummy(self) -> None:
        """
        Draw a simple dummy plot.
        Useful as a placeholder while wiring up the real logic.
        """
        import numpy as np

        self.axes.cla()
        x = np.linspace(0, 2 * np.pi, 300)
        y = np.sin(x)
        self.axes.plot(x, y, label="sin(x)")
        self.axes.set_xlabel("x")
        self.axes.set_ylabel("sin(x)")
        self.axes.legend()
        self.draw()


class MatplotlibWidget(QWidget):
    """
    A QWidget that wraps a MatplotlibCanvas and its NavigationToolbar.

    This is what you typically add to your page layouts.
    """

    def __init__(
        self,
        parent: Optional[QWidget] = None,
        canvas_size: Tuple[float, float] = (6.0, 4.0),
        dpi: int = 100,
    ) -> None:
        super().__init__(parent)

        width, height = canvas_size
        self.canvas = MatplotlibCanvas(self, width=width, height=height, dpi=dpi)
        self.toolbar = NavigationToolbar2QT(self.canvas, self)

        layout = QVBoxLayout()
        layout.setContentsMargins(0, 0, 0, 0)
        layout.setSpacing(2)

        layout.addWidget(self.toolbar)
        layout.addWidget(self.canvas)

        self.setLayout(layout)

    def clear(self) -> None:
        """Clear the canvas axes."""
        self.canvas.clear()

    def plot_dummy(self) -> None:
        """Forward to the canvas dummy plot."""
        self.canvas.plot_dummy()
