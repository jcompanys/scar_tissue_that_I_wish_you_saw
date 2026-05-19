"""Backward-compatible wrapper for scar characterization and analysis helpers.

Prefer importing ``scar_characterization`` for descriptor creation and
``scar_analysis`` for statistical tests/associations in new code.
"""

from src.scar_analysis import *  # noqa: F401,F403
from src.scar_characterization import *  # noqa: F401,F403
