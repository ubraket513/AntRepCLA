"""Matplotlib backend selection.

Imported before pyplot anywhere else in the package. A cluster batch job has no
display, and calling ``plt.show()`` there blocks until the job is killed -- a
failure that previously cost more wall time than the entire analysis.
"""

from __future__ import annotations

import os
import sys

import matplotlib


def has_display() -> bool:
    """Whether an interactive plot window could actually be shown.

    macOS is treated as interactive because it does not use DISPLAY.
    """
    if not sys.stdout.isatty():
        return False
    if sys.platform == "darwin":
        return True
    return bool(os.environ.get("DISPLAY") or os.environ.get("WAYLAND_DISPLAY"))


if not has_display():
    matplotlib.use("Agg")
