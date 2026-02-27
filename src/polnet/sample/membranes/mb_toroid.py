"""Toroidal membrane generator.

Defines :class:`TorGen`, registered under ``"toroid"`` in
:class:`~polnet.sample.membranes.mb_factory.MbFactory`.
Generates randomly oriented toroidal (donut-shaped) lipid
bilayers for modelling organelle membranes with hole topology.

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
    iso_surface,
    lin_map,
    poly_threshold,
)


@MbFactory.register("toroid")
class TorGen(MbGen):
    """Toroidal membrane generator."""

    def __init__(
        self,
        thick_rg: tuple[float, float],
        layer_s_rg: tuple[float, float],
        occ_rg: tuple[float, float],
        over_tol: float,
        mb_den_cf_rg: tuple[float, float],
        min_rad: float,
        max_rad: float = None,
    ) -> None:
        """Constructor.

        Args:
            thick_rg: (min, max) bilayer thickness in angstroms.
            layer_s_rg: (min, max) Gaussian layer sigma in angstroms.
            occ_rg: (min, max) target occupancy fraction.
            over_tol: Overlap tolerance fraction.
            mb_den_cf_rg: (min, max) density contrast factor.
            min_rad: Minimum radius (both major and minor) in angstroms.
            max_rad: Maximum radius in angstroms (None = auto from VOI diagonal).
        """
        super().__init__(
            thick_rg=thick_rg,
            layer_s_rg=layer_s_rg,
            occ_rg=occ_rg,
            over_tol=over_tol,
            mb_den_cf_rg=mb_den_cf_rg,
        )
        self._min_radius = min_rad
        self._max_radius = max_rad

    @classmethod
    def from_params(cls, params: dict) -> "TorGen":
        """Create a TorGen from a parameter dictionary.

        Args:
            params: Dictionary with membrane parameters.

        Returns:
            TorGen instance.
        """
        return cls(
            thick_rg=params.get("MB_THICK_RG", (25.0, 35.0)),
            layer_s_rg=params.get("MB_LAYER_S_RG", (0.5, 2.0)),
            occ_rg=params.get("MB_OCC_RG", (0.001, 0.003)),
            over_tol=params.get("MB_OVER_TOL", 0.0),
            mb_den_cf_rg=params.get("MB_DEN_CF_RG", (0.3, 0.5)),
            min_rad=params.get("MB_MIN_RAD", 75.0),
            max_rad=params.get("MB_MAX_RAD", None),
        )

    def _build(self, voi_shape: tuple[int, int, int], v_size: float) -> Mb:
        """Generate a single toroidal membrane with random parameters.

        Args:
            voi_shape: Volume shape (X, Y, Z) in voxels.
            v_size: Voxel size in angstroms.

        Returns:
            Mb: Constructed membrane.

        Raises:
            MbError: If generated membrane has zero volume.
        """
        max_radius = (
            self._max_radius
            if self._max_radius is not None
            else math.sqrt(3) * max(voi_shape) * v_size
        )

        center = np.array(
            [
                voi_shape[0] * v_size * random.random(),
                voi_shape[1] * v_size * random.random(),
                voi_shape[2] * v_size * random.random(),
            ]
        )
        thick = random.uniform(self._thick_rg[0], self._thick_rg[1])
        layer_s = random.uniform(self._layer_s_rg[0], self._layer_s_rg[1])
        rot_q = gen_rand_unit_quaternion()

        major_radius = random.uniform(self._min_radius, max_radius)
        minor_radius = random.uniform(self._min_radius, max_radius)
        if minor_radius >= major_radius:
            major_radius, minor_radius = minor_radius, major_radius

        # --- geometry (unchanged) ---

        t_v = 0.5 * thick / v_size
        s_v = layer_s / v_size
        rad_a_v = major_radius / v_size
        rad_b_v = minor_radius / v_size
        bo_v, bi_v = rad_b_v + t_v, rad_b_v - t_v

        if bi_v <= 0:
            raise MbError(
                f"Inner tube radius bi_v={bi_v:.2f} vx is non-positive "
                f"(minor_radius={minor_radius:.1f} A, thick={thick:.1f} A, "
                f"v_size={v_size:.1f} A). Minor radius must exceed half "
                "the bilayer thickness."
            )
        bo_v_p1, bo_v_m1 = bo_v + 1, bo_v - 1
        bi_v_p1, bi_v_m1 = bi_v + 1, bi_v - 1
        p0_v = center / v_size

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
            (rad_a_v - np.sqrt((X - p0_v[0]) ** 2 + (Y - p0_v[1]) ** 2)) ** 2
            + (Z - p0_v[2]) ** 2
            - bo_v * bo_v
        ) <= 1
        R_i = (
            (rad_a_v - np.sqrt((X - p0_v[0]) ** 2 + (Y - p0_v[1]) ** 2)) ** 2
            + (Z - p0_v[2]) ** 2
            - bi_v * bi_v
        ) >= 1
        mask = tomo_rotate(np.logical_and(R_i, R_o), rot_q, order=0)

        if mask.sum() == 0:
            raise MbError("Generated membrane has zero volume.")

        R_i = (
            (rad_a_v - np.sqrt((X - p0_v[0]) ** 2 + (Y - p0_v[1]) ** 2)) ** 2
            + (Z - p0_v[2]) ** 2
            - rad_b_v * rad_b_v
        )
        R_i = tomo_rotate(R_i, rot_q, mode="reflect")
        surf = iso_surface(R_i, 1)
        add_sfield_to_poly(
            surf, mask, "mb_mask", dtype="int", interp="NN", mode="points"
        )
        surf = poly_threshold(surf, "mb_mask", mode="points", low_th=0.5)

        # Outer layer
        R_o = (
            (rad_a_v - np.sqrt((X - p0_v[0]) ** 2 + (Y - p0_v[1]) ** 2)) ** 2
            + (Z - p0_v[2]) ** 2
            - bo_v_p1 * bo_v_p1
        ) <= 1
        R_i = (
            (rad_a_v - np.sqrt((X - p0_v[0]) ** 2 + (Y - p0_v[1]) ** 2)) ** 2
            + (Z - p0_v[2]) ** 2
            - bo_v_m1 * bo_v_m1
        ) >= 1
        G = tomo_rotate(np.logical_and(R_i, R_o), rot_q, order=0)

        # Inner layer
        R_o = (
            (rad_a_v - np.sqrt((X - p0_v[0]) ** 2 + (Y - p0_v[1]) ** 2)) ** 2
            + (Z - p0_v[2]) ** 2
            - bi_v_p1 * bi_v_p1
        ) <= 1
        R_i = (
            (rad_a_v - np.sqrt((X - p0_v[0]) ** 2 + (Y - p0_v[1]) ** 2)) ** 2
            + (Z - p0_v[2]) ** 2
            - bi_v_m1 * bi_v_m1
        ) >= 1
        G += tomo_rotate(np.logical_and(R_i, R_o), rot_q, order=0)

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
