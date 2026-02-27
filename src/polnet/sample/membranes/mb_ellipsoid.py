"""Ellipsoidal membrane generator.

Defines :class:`EllipGen`, registered under ``"ellipsoid"`` in
:class:`~polnet.sample.membranes.mb_factory.MbFactory`.
Generates randomly oriented ellipsoidal lipid bilayers with
eccentricity drawn from a bounded exponential distribution.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import math
import random

import numpy as np
from scipy.ndimage import gaussian_filter

from .mb import (
    Mb,
    MbError,
    MbGen,
)
from .mb_factory import MbFactory
from ...utils.affine import (
    gen_rand_unit_quaternion,
    tomo_rotate,
)
from ...utils.poly import (
    add_sfield_to_poly,
)
from ...utils.utils import (
    density_norm,
    gen_bounded_exp,
    iso_surface,
    lin_map,
    poly_threshold,
)


@MbFactory.register("ellipsoid")
class EllipGen(MbGen):
    """Ellipsoidal membrane generator."""

    _ECC_MAX_TRIES = int(1e6)

    def __init__(
        self,
        thick_rg: tuple[float, float],
        layer_s_rg: tuple[float, float],
        occ_rg: tuple[float, float],
        over_tol: float,
        mb_den_cf_rg: tuple[float, float],
        max_ecc: float,
        min_axis: float,
        max_axis: float = None,
    ) -> None:
        """Constructor.

        Args:
            thick_rg: (min, max) bilayer thickness in angstroms.
            layer_s_rg: (min, max) Gaussian layer sigma in angstroms.
            occ_rg: (min, max) target occupancy fraction.
            over_tol: Overlap tolerance fraction.
            mb_den_cf_rg: (min, max) density contrast factor.
            max_ecc: Maximum eccentricity constraint.
            min_axis: Minimum semi-axis length in angstroms.
            max_axis: Maximum semi-axis length in angstroms (None = auto from VOI diagonal).
        """
        super().__init__(
            thick_rg=thick_rg,
            layer_s_rg=layer_s_rg,
            occ_rg=occ_rg,
            over_tol=over_tol,
            mb_den_cf_rg=mb_den_cf_rg,
        )
        self._max_ecc = max_ecc
        self._min_axis = min_axis
        self._max_axis = max_axis

    @classmethod
    def from_params(cls, params: dict) -> "EllipGen":
        """Create an EllipGen from a parameter dictionary.

        Args:
            params: Dictionary with membrane parameters.

        Returns:
            EllipGen instance.
        """
        return cls(
            thick_rg=params.get("MB_THICK_RG", (25.0, 35.0)),
            layer_s_rg=params.get("MB_LAYER_S_RG", (0.5, 2.0)),
            occ_rg=params.get("MB_OCC_RG", (0.001, 0.003)),
            over_tol=params.get("MB_OVER_TOL", 0.0),
            mb_den_cf_rg=params.get("MB_DEN_CF_RG", (0.3, 0.5)),
            max_ecc=params.get("MB_MAX_ECC", 0.75),
            min_axis=params.get("MB_MIN_AXIS", 75.0),
            max_axis=params.get("MB_MAX_AXIS", None),
        )

    def _build(self, voi_shape: tuple[int, int, int], v_size: float) -> Mb:
        """Generate a single ellipsoidal membrane with random parameters.

        Args:
            voi_shape: Volume shape (X, Y, Z) in voxels.
            v_size: Voxel size in angstroms.

        Returns:
            Mb: Constructed membrane.

        Raises:
            MbError: If eccentricity constraint cannot be satisfied or membrane has zero volume.
        """
        max_axis = (
            self._max_axis
            if self._max_axis is not None
            else math.sqrt(3) * max(voi_shape) * v_size
        )

        center = tuple(
            (
                voi_shape[0] * v_size * random.random(),
                voi_shape[1] * v_size * random.random(),
                voi_shape[2] * v_size * random.random(),
            )
        )
        thick = random.uniform(self._thick_rg[0], self._thick_rg[1])
        layer_s = random.uniform(self._layer_s_rg[0], self._layer_s_rg[1])
        rot_q = gen_rand_unit_quaternion()

        # Find semi-axes satisfying eccentricity constraint
        for _ in range(self._ECC_MAX_TRIES):
            axes = np.sort(
                np.array(
                    [
                        gen_bounded_exp(
                            8.0 * self._min_axis, self._min_axis, max_axis
                        ),
                        gen_bounded_exp(
                            8.0 * self._min_axis, self._min_axis, max_axis
                        ),
                        gen_bounded_exp(
                            8.0 * self._min_axis, self._min_axis, max_axis
                        ),
                    ]
                )
            )[::-1]
            ecc1 = math.sqrt(1 - (axes[1] / axes[0]) ** 2)
            ecc2 = math.sqrt(1 - (axes[2] / axes[0]) ** 2)
            if ecc1 <= self._max_ecc and ecc2 <= self._max_ecc:
                break
        else:
            raise MbError(
                "Could not generate ellipsoid with the specified maximum eccentricity"
            )

        a, b, c = float(axes[0]), float(axes[1]), float(axes[2])
        center_arr = np.array([float(x) for x in center])

        # --- geometry (unchanged) ---

        t_v = 0.5 * thick / v_size
        s_v = layer_s / v_size
        a_v, b_v, c_v = a / v_size, b / v_size, c / v_size
        ao_v, bo_v, co_v = a_v + t_v, b_v + t_v, c_v + t_v
        ai_v, bi_v, ci_v = a_v - t_v, b_v - t_v, c_v - t_v
        ao_v_p1, bo_v_p1, co_v_p1 = ao_v + 1, bo_v + 1, co_v + 1
        ao_v_m1, bo_v_m1, co_v_m1 = ao_v - 1, bo_v - 1, co_v - 1
        ai_v_p1, bi_v_p1, ci_v_p1 = ai_v + 1, bi_v + 1, ci_v + 1
        ai_v_m1, bi_v_m1, ci_v_m1 = ai_v - 1, bi_v - 1, ci_v - 1
        p0_v = center_arr / v_size

        dx, dy, dz = (
            float(voi_shape[0]),
            float(voi_shape[1]),
            float(voi_shape[2]),
        )
        dx2, dy2, dz2 = (
            math.floor(0.5 * dx),
            math.floor(0.5 * dy),
            math.floor(0.5 * dz),
        )
        p0_v[0] -= dx2
        p0_v[1] -= dy2
        p0_v[2] -= dz2
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(
            np.arange(x_l, x_h),
            np.arange(y_l, y_h),
            np.arange(z_l, z_h),
            indexing="ij",
        )

        R_o = (
            ((X - p0_v[0]) / ao_v) ** 2
            + ((Y - p0_v[1]) / bo_v) ** 2
            + ((Z - p0_v[2]) / co_v) ** 2
        )
        R_i = (
            ((X - p0_v[0]) / ai_v) ** 2
            + ((Y - p0_v[1]) / bi_v) ** 2
            + ((Z - p0_v[2]) / ci_v) ** 2
        )
        mask = tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), rot_q, order=0)
        if mask.sum() == 0:
            raise MbError(
                "Generated membrane has zero volume. Try changing its parameters."
            )

        R_i = (
            ((X - p0_v[0]) / a_v) ** 2
            + ((Y - p0_v[1]) / b_v) ** 2
            + ((Z - p0_v[2]) / c_v) ** 2
        )
        R_i = tomo_rotate(R_i, rot_q, mode="reflect")
        surf = iso_surface(R_i, 1)
        add_sfield_to_poly(
            surf, mask, "mb_mask", dtype="int", interp="NN", mode="points"
        )
        surf = poly_threshold(surf, "mb_mask", mode="points", low_th=0.5)

        # Outer layer
        R_o = (
            ((X - p0_v[0]) / ao_v_p1) ** 2
            + ((Y - p0_v[1]) / bo_v_p1) ** 2
            + ((Z - p0_v[2]) / co_v_p1) ** 2
        )
        R_i = (
            ((X - p0_v[0]) / ao_v_m1) ** 2
            + ((Y - p0_v[1]) / bo_v_m1) ** 2
            + ((Z - p0_v[2]) / co_v_m1) ** 2
        )
        G = tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), rot_q, order=0)

        # Inner layer
        R_o = (
            ((X - p0_v[0]) / ai_v_p1) ** 2
            + ((Y - p0_v[1]) / bi_v_p1) ** 2
            + ((Z - p0_v[2]) / ci_v_p1) ** 2
        )
        R_i = (
            ((X - p0_v[0]) / ai_v_m1) ** 2
            + ((Y - p0_v[1]) / bi_v_m1) ** 2
            + ((Z - p0_v[2]) / ci_v_m1) ** 2
        )
        G += tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), rot_q, order=0)

        if G.sum() == 0:
            raise MbError(
                "Generated membrane bilayer has no voxels within the VOI."
            )

        density = lin_map(
            density_norm(gaussian_filter(G.astype(float), s_v), inv=True),
            ub=0,
            lb=1,
        )

        return Mb(
            voi_shape=voi_shape,
            v_size=v_size,
            thick=thick,
            layer_s=layer_s,
            density=density,
            mask=mask,
            surf=surf,
        )
