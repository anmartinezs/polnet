"""Membranes sub-package: lipid bilayer generators.

Re-exports the membrane data holder, abstract generator, factory,
configuration file loader, and all registered concrete generators
(sphere, ellipsoid, toroid, curvatubes).

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

from .mb import (
    Mb,
    MbError,
    MbGen,
    MbSetResult,
)
from .mb_factory import MbFactory
from .mb_file import MbFile

# Needed for the dynamic registration of the generators in the factory
from .mb_ellipsoid import EllipGen
from .mb_sphere import SphGen
from .mb_toroid import TorGen

try:
    from .mb_curvatubes import CurvatubesGen  # requires external/cvtub
except ImportError:
    CurvatubesGen = None  # type: ignore[assignment,misc]

__all__ = ["Mb", "MbGen", "MbSetResult", "MbFactory", "MbFile"]
