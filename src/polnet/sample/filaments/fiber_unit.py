"""Abstract fiber unit and concrete implementations.

Defines :class:`FiberUnit` (abstract), :class:`FiberUnitSDimer`
(actin-like S-dimer), and :class:`MTUnit` (microtubule-like unit).
Each subclass generates the VTK surface and density sub-volume
for a single polymer repeat unit.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import math
from abc import (
    ABC,
    abstractmethod,
)

import numpy as np

from ...utils.affine import (
    poly_scale,
    poly_translate,
)
from ...utils.utils import (
    lin_map,
    iso_surface,
)


class FiberUnit(ABC):
    """Abstract base class for fiber structural units.

    Concrete subclasses generate the VTK surface and density
    sub-volume for a single polymer repeat unit.
    """

    @property
    @abstractmethod
    def vtp(self):
        """VTK polygon surface of the structural unit."""
        raise NotImplementedError

    @property
    @abstractmethod
    def tomo(self):
        """3-D density sub-volume of the structural unit."""
        raise NotImplementedError


class FiberUnitSDimer(FiberUnit):
    """Actin-like fiber unit modelled as a sphere dimer.

    The two overlapping spheres are defined via a logistic
    function on a 3-D grid, producing a smooth density envelope
    suitable for cryo-ET simulation.
    """

    def __init__(self, sph_rad, v_size=1):
        """Initialise the sphere-dimer fiber unit.

        Args:
            sph_rad (float): Sphere radius in Angstroms; must be
                > 0.
            v_size (float): Voxel size in Angstroms (default 1);
                must be > 0.

        Raises:
            ValueError: If sph_rad <= 0 or v_size <= 0.
        """
        if sph_rad <= 0 or v_size <= 0:
            raise ValueError("sph_rad and v_size must be > 0.")
        self.__sph_rad, self.__v_size = float(sph_rad), float(v_size)
        self.__size = int(math.ceil(6.0 * (sph_rad / v_size)))
        if self.__size % 2 != 0:
            self.__size += 1
        self.__tomo, self.__surf = None, None
        self.__gen_sdimer()

    @property
    def vtp(self):
        """VTK polygon surface of the sphere-dimer unit."""
        return self.__surf

    @property
    def tomo(self):
        """3-D density sub-volume of the sphere-dimer unit."""
        return self.__tomo

    def __gen_sdimer(self):
        """Generate a sphere-dimer density envelope using logistic functions.

        Two overlapping spheres of radius ``sph_rad`` are placed
        along the Y axis and blended via logistic fall-off.
        The resulting density is normalised and an isosurface
        extracted.
        """

        sph_rad_v = self.__sph_rad / self.__v_size
        sph_rad_v2 = sph_rad_v * sph_rad_v * 0.5625  # (0.75*rad)^2

        self.__tomo = np.zeros(
            shape=(self.__size, self.__size, self.__size), dtype=np.float32
        )
        dx, dy, dz = (
            float(self.__tomo.shape[0]),
            float(self.__tomo.shape[1]),
            float(self.__tomo.shape[2]),
        )
        dx2, dy2, dz2 = (
            math.floor(0.5 * dx),
            math.floor(0.5 * dy),
            math.floor(0.5 * dz),
        )
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(
            np.arange(x_l, x_h),
            np.arange(y_l, y_h),
            np.arange(z_l, z_h),
            indexing="xy",
        )
        X += 0.5
        Y += 0.5
        Z += 0.5

        Yh = Y + sph_rad_v
        R = X * X + Yh * Yh + Z * Z - sph_rad_v2

        self.__tomo += 1.0 / (1.0 + np.exp(-R))

        Yh = Y - sph_rad_v
        R = X * X + Yh * Yh + Z * Z - sph_rad_v2
        self.__tomo += 1.0 / (1.0 + np.exp(-R))

        self.__tomo = lin_map(self.__tomo, lb=1, ub=0)
        self.__surf = iso_surface(self.__tomo, 0.25)
        self.__surf = poly_scale(self.__surf, self.__v_size)
        self.__surf = poly_translate(
            self.__surf,
            -0.5 * self.__v_size * (np.asarray(self.__tomo.shape) - 0.5),
        )


class MTUnit(FiberUnit):
    """Microtubule-like fiber unit composed of n_units protofilaments.

    Each protofilament is represented as a sphere placed at the
    appropriate angular position on a circle of radius mt_rad.
    The logistic-function density envelope approximates the
    hollow cylindrical cross-section of a microtubule.
    """

    def __init__(self, sph_rad=40, mt_rad=100.5, n_units=13, v_size=1):
        """Initialise the microtubule fiber unit.

        Args:
            sph_rad (float): Tubulin monomer sphere radius in
                Angstroms (default 40).
            mt_rad (float): Microtubule protofilament ring radius
                in Angstroms (default 100.5).
            n_units (int): Number of protofilaments (default 13).
            v_size (float): Voxel size in Angstroms (default 1).

        Raises:
            ValueError: If any parameter is <= 0.
        """
        if sph_rad <= 0 or mt_rad <= 0 or n_units <= 0 or v_size <= 0:
            raise ValueError(
                "sph_rad, mt_rad, n_units and " "v_size must be > 0."
            )
        self.__sph_rad, self.__mt_rad, self.__n_units, self.__v_size = (
            float(sph_rad),
            float(mt_rad),
            int(n_units),
            float(v_size),
        )
        self.__size = int(math.ceil(6.0 * (sph_rad / v_size))) - 6
        if self.__size % 2 != 0:
            self.__size += 1
        self.__tomo, self.__surf = None, None
        self.__gen_sdimer()

    @property
    def vtp(self):
        """VTK polygon surface of the microtubule unit."""
        return self.__surf

    @property
    def tomo(self):
        """3-D density sub-volume of the microtubule unit."""
        return self.__tomo

    def __gen_sdimer(self):
        """Generate a protofilament-ring density using logistic functions.

        Places ``n_units`` spheres of radius ``sph_rad`` equally
        spaced around a circle of radius ``mt_rad``, blending
        each with a logistic fall-off.  The resulting density is
        normalised and an isosurface extracted.
        """

        sph_rad_v, mt_rad_v = (
            self.__sph_rad / self.__v_size,
            self.__mt_rad / self.__v_size,
        )
        sph_rad_v2 = sph_rad_v * sph_rad_v * 0.25  # (0.9*rad)^2

        self.__tomo = np.zeros(
            shape=(self.__size, self.__size, self.__size), dtype=np.float32
        )
        dx, dy, dz = (
            float(self.__tomo.shape[0]),
            float(self.__tomo.shape[1]),
            float(self.__tomo.shape[2]),
        )
        dx2, dy2, dz2 = (
            math.floor(0.5 * dx),
            math.floor(0.5 * dy),
            math.floor(0.5 * dz),
        )
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(
            np.arange(x_l, x_h),
            np.arange(y_l, y_h),
            np.arange(z_l, z_h),
            indexing="xy",  # Indexing xy inside the local grid
        )
        X += 0.5
        Y += 0.5
        Z += 0.5

        Z2 = Z * Z
        ang_step = 2.0 * np.pi / self.__n_units
        ang = 0
        while ang <= 2.0 * np.pi:
            x, y = mt_rad_v * math.cos(ang), mt_rad_v * math.sin(ang)
            Xh, Yh = X + x, Y + y
            R = Xh * Xh + Yh * Yh + Z2 - sph_rad_v2
            F = 1.0 / (1.0 + np.exp(-R))
            self.__tomo += -F + 1
            ang += ang_step

        # Generating the surfaces
        self.__tomo = lin_map(self.__tomo, lb=0, ub=1)
        self.__surf = iso_surface(self.__tomo, 0.25)
        self.__surf = poly_scale(self.__surf, self.__v_size)
        self.__surf = poly_translate(
            self.__surf,
            -0.5 * self.__v_size * (np.asarray(self.__tomo.shape) - 0.5),
        )
