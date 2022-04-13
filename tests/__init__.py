"""Module containing some utilities for testing."""

from pathlib import Path
import pkg_resources as pkg

__all__ = ["PATH_QMCBLIP", "PATH_TEST"]

# Environment data
PATH_QMCBLIP = Path(pkg.resource_filename('qmcblip', ''))
ROOT = PATH_QMCBLIP.parent

PATH_TEST = ROOT / "tests"
