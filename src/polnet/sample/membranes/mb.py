"""Module for membrane modeling classes.

Defines:
    Mb: Concrete membrane data holder with insertion methods.
    MbGen: Abstract generator with _build() (per-shape) and generate_set() (shared loop).
    MbSetResult: Lightweight result container for generate_set().
    MbError: Custom exception for membrane-related errors.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

from abc import (
    ABC,
    abstractmethod,
)
import random
from typing import NamedTuple

import numpy as np
import vtk
from scipy.ndimage import binary_dilation

from ...logging_conf import _LOGGER as logger
from ...utils.affine import (
    poly_scale,
)
from ...utils.poly import (
    poly_mask,
)
from ...utils.utils import (
    insert_svol_tomo,
)


class Mb:
    """A single membrane: density volume, binary mask, and VTK surface.

    Constructed by MbGen._build(), not by user code directly.
    """

    def __init__(
        self,
        voi_shape: tuple[int, int, int],
        v_size: float,
        thick: float,
        layer_s: float,
        density: np.ndarray,
        mask: np.ndarray,
        surf: vtk.vtkPolyData,
    ) -> None:
        """Constructor.

        Args:
            voi_shape: Reference volume shape (X, Y, Z) in voxels.
            v_size: Voxel size in angstroms.
            thick: Membrane bilayer thickness in angstroms.
            layer_s: Gaussian sigma for each layer in angstroms.
            density: 3D float32 array with the membrane density profile.
            mask: 3D bool array — binary occupancy mask.
            surf: VTK polydata surface at the membrane midplane.

        Raises:
            TypeError: If voi_shape does not have length 3.
            ValueError: If v_size is not positive or array shapes don't match voi_shape.
        """
        if not hasattr(voi_shape, "__len__") or len(voi_shape) != 3:
            raise TypeError("voi_shape must have length 3")
        if v_size <= 0:
            raise ValueError("v_size must be positive")
        if density.shape != tuple(voi_shape) or mask.shape != tuple(voi_shape):
            raise ValueError("density and mask shapes must match voi_shape")

        self._voi_shape = tuple(voi_shape)
        self._v_size = float(v_size)
        self._thick = float(thick)
        self._layer_s = float(layer_s)
        self._density = density
        self._mask = mask
        self._surf = surf

    @property
    def thick(self) -> float:
        """Membrane bilayer thickness in angstroms."""
        return self._thick

    @property
    def layer_s(self) -> float:
        """Gaussian sigma for each layer in angstroms."""
        return self._layer_s

    @property
    def vol(self) -> float:
        """Membrane volume in cubic angstroms."""
        return float(self._mask.sum()) * self._v_size**3

    @property
    def density(self) -> np.ndarray:
        """3D density array (copy)."""
        return self._density.copy()

    @property
    def mask(self) -> np.ndarray:
        """3D binary mask (copy)."""
        return self._mask.copy()

    @property
    def vtp(self) -> vtk.vtkPolyData:
        """VTK polydata surface."""
        return self._surf

    def masking(self, ext_mask: np.ndarray) -> None:
        """Zero out membrane voxels where ext_mask is False.

        Args:
            ext_mask: Bool array with the same shape as voi_shape.

        Raises:
            TypeError: If ext_mask is not a boolean ndarray.
            ValueError: If ext_mask shape doesn't match voi_shape.
        """
        if not isinstance(ext_mask, np.ndarray) or ext_mask.dtype != bool:
            raise TypeError("ext_mask must be a boolean ndarray")
        if ext_mask.shape != self._voi_shape:
            raise ValueError("ext_mask shape must match voi_shape")
        self._density[~ext_mask] = 0
        self._mask[~ext_mask] = False
        self._surf = poly_mask(self._surf, ext_mask)

    def insert_density_svol(
        self,
        tomo: np.ndarray,
        merge: str = "max",
        mode: str = "tomo",
        grow: int = 0,
    ) -> None:
        """Stamp this membrane into a target volume (in place).

        Args:
            tomo: Target 3D array, modified in place.
            merge: Merge strategy — 'min', 'max', 'sum', or 'insert'.
            mode: Data source — 'tomo' (density), 'mask', or 'voi' (inverted mask).
            grow: Dilation iterations for 'voi' mode.

        Raises:
            TypeError: If tomo is not a 3D ndarray.
            ValueError: If shape mismatch or invalid merge/mode.
        """
        if not isinstance(tomo, np.ndarray) or tomo.ndim != 3:
            raise TypeError("tomo must be a 3D ndarray")
        if tomo.shape != self._voi_shape:
            raise ValueError("tomo shape must match voi_shape")
        if merge not in ("min", "max", "sum", "insert"):
            raise ValueError(f"Invalid merge mode: {merge}")
        if mode not in ("tomo", "mask", "voi"):
            raise ValueError(f"Invalid mode: {mode}")

        if mode == "tomo":
            hold = self._density
        elif mode == "mask":
            hold = self._mask
        else:
            hold = (
                ~binary_dilation(self._mask, iterations=grow)
                if grow >= 1
                else ~self._mask
            )
        insert_svol_tomo(
            hold, tomo, 0.5 * np.asarray(self._voi_shape), merge=merge
        )


class MbSetResult(NamedTuple):
    """Output of MbGen.generate_set()."""

    voi: np.ndarray
    density: np.ndarray
    mask: np.ndarray
    vtp: vtk.vtkPolyData
    num_mbs: int
    mb_occupancy: float


class MbGen(ABC):
    """Abstract membrane generator.

    Subclasses implement ``_build()`` for a single membrane.
    The concrete ``generate_set()`` handles the occupancy loop, overlap
    checking, and density/mask/VOI accumulation.
    """

    def __init__(
        self,
        thick_rg: tuple[float, float],
        layer_s_rg: tuple[float, float],
        occ_rg: tuple[float, float],
        over_tol: float,
        mb_den_cf_rg: tuple[float, float],
    ) -> None:
        """Constructor.

        Args:
            thick_rg: (min, max) bilayer thickness in angstroms.
            layer_s_rg: (min, max) Gaussian layer sigma in angstroms.
            occ_rg: (min, max) target membrane occupancy fraction.
            over_tol: Overlap tolerance as a fraction of membrane voxels.
            mb_den_cf_rg: (min, max) density contrast factor.

        Raises:
            ValueError: If any range is invalid.
        """
        if thick_rg[0] <= 0 or thick_rg[1] <= 0:
            raise ValueError("thick_rg values must be positive")
        if thick_rg[0] > thick_rg[1]:
            raise ValueError("thick_rg must be (min, max)")
        if layer_s_rg[0] < 0 or layer_s_rg[1] < 0:
            raise ValueError("layer_s_rg values must be non-negative")
        if layer_s_rg[0] > layer_s_rg[1]:
            raise ValueError("layer_s_rg must be (min, max)")
        if occ_rg[0] < 0 or occ_rg[1] < 0:
            raise ValueError("occ_rg values must be non-negative")
        if occ_rg[0] > occ_rg[1]:
            raise ValueError("occ_rg must be (min, max)")
        if mb_den_cf_rg[0] < 0 or mb_den_cf_rg[1] < 0:
            raise ValueError("mb_den_cf_rg values must be non-negative")
        if mb_den_cf_rg[0] > mb_den_cf_rg[1]:
            raise ValueError("mb_den_cf_rg must be (min, max)")

        self._thick_rg = thick_rg
        self._layer_s_rg = layer_s_rg
        self._occ_rg = occ_rg
        self._over_tol = over_tol
        self._mb_den_cf_rg = mb_den_cf_rg

    def rnd_occ(self) -> float:
        """Random occupancy within the configured range."""
        return random.uniform(self._occ_rg[0], self._occ_rg[1])

    def rnd_cf(self) -> float:
        """Random density contrast factor within the configured range."""
        return random.uniform(self._mb_den_cf_rg[0], self._mb_den_cf_rg[1])

    @property
    def over_tolerance(self) -> float:
        """Overlap tolerance fraction."""
        return self._over_tol

    @classmethod
    @abstractmethod
    def from_params(cls, params: dict) -> "MbGen":
        """Create a generator from a parameter dictionary."""
        raise NotImplementedError

    @abstractmethod
    def _build(self, voi_shape: tuple[int, int, int], v_size: float) -> Mb:
        """Generate a single membrane with random parameters.

        Args:
            voi_shape: Volume shape in voxels.
            v_size: Voxel size in angstroms.

        Returns:
            A fully constructed Mb instance.
        """
        raise NotImplementedError

    def generate_set(
        self,
        *,
        voi: np.ndarray,
        v_size: float,
        bg_voi: np.ndarray | None = None,
        max_mbtries: int = 10,
        grow: int = 0,
    ) -> MbSetResult:
        """Build a set of membranes until the occupancy target is reached.

        Uses _build() to generate individual membranes, checks overlap,
        and accumulates density/mask/VOI/VTP. All state is local — the
        generator instance is not mutated.

        Args:
            voi: 3D bool array — volume of interest (copied internally).
            v_size: Voxel size in angstroms.
            bg_voi: Optional background VOI mask for additional masking.
            max_mbtries: Maximum consecutive failures before stopping.
            grow: VOI dilation iterations per inserted membrane.

        Returns:
            MbSetResult with the final voi, density, mask, vtp, count, and occupancy.
        """
        voi = voi.copy()
        shape = voi.shape
        density = np.zeros(shape, dtype=np.float32)
        mask = np.zeros(shape, dtype=bool)
        app_vtp = vtk.vtkAppendPolyData()
        count_mbs = 0

        cf = self.rnd_cf()
        max_occ = self.rnd_occ()
        over_tol = self._over_tol

        def _occupancy():
            voi_sum = voi.sum()
            return mask.sum() / voi_sum if voi_sum > 0 else 0.0

        def _overlaps(mb):
            mb_mask = mb.mask
            available = np.logical_and(voi, ~mask)
            overlap = np.logical_and(mb_mask, ~available)
            mb_voxels = mb_mask.sum()
            return (
                (overlap.sum() / mb_voxels) > over_tol
                if mb_voxels > 0
                else False
            )

        count_fails = 0
        while _occupancy() < max_occ:
            try:
                mb = self._build(shape, v_size)

                if bg_voi is not None:
                    mb.masking(bg_voi)
                if mb.vol <= 0:
                    raise MbError("Generated membrane has zero volume.")
                if _overlaps(mb):
                    raise MbError("Membrane overlaps with existing set.")

                mb.insert_density_svol(density, merge="max", mode="tomo")
                mb.insert_density_svol(mask, merge="max", mode="mask")
                mb.insert_density_svol(voi, merge="min", mode="voi", grow=grow)
                app_vtp.AddInputData(mb.vtp)
                count_mbs += 1
                count_fails = 0

            except MbError as e:
                logger.debug(
                    "Membrane generation try %d failed: %s",
                    count_fails + 1,
                    str(e),
                )
                count_fails += 1
                if count_fails >= max_mbtries:
                    logger.warning(
                        "Failed to insert membrane after %d consecutive "
                        "attempts.",
                        max_mbtries,
                    )
                    break
                continue

        density *= cf

        if count_mbs > 0:
            app_vtp.Update()
            surfs = poly_scale(app_vtp.GetOutput(), v_size)
        else:
            surfs = vtk.vtkPolyData()

        return MbSetResult(
            voi=voi,
            density=density,
            mask=mask,
            vtp=surfs,
            num_mbs=count_mbs,
            mb_occupancy=_occupancy(),
        )


class MbError(Exception):
    """Custom exception for membrane-related errors."""

    def __init__(self, message: str) -> None:
        super().__init__(message)
