"""Helical flexible fiber polymer model.

Defines :class:`HelixFiber`, a
:class:`~polnet.sample.polymers.Polymer` subclass that models a
random worm-like helical chain.  Each monomer is placed following
helix geometry with persistence-length-controlled bending.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import math

import numpy as np

from ..polymers import (
    Monomer,
    Polymer,
)
from ...utils.affine import (
    angle_axis_to_quat,
    quat_mult,
    rot_to_quat,
    rot_vect_quat,
    vector_module,
    vect_to_zmat,
    wrap_angle,
)
from ...utils.utils import (
    gen_uni_s2_sample,
    points_distance,
)


class HelixFiber(Polymer):
    """Random flexible helical fiber polymer.

    Models a worm-like chain (WLC) with helical geometry and
    persistence-length-controlled bending.  The helix axis and
    azimuthal step size are derived from the persistence length,
    helix period, and monomer z-length.
    """

    def __init__(
        self,
        l_length,
        m_surf,
        p_length,
        hp_length,
        mz_length,
        z_length_f=0,
        p0=(0, 0, 0),
        vz=(0, 0, 1),
        rot_rand=True,
    ):
        """Initialise a helical flexible fiber.

        Args:
            l_length (float): Monomer link length in Angstroms;
                must be > 0.
            m_surf (vtk.vtkPolyData): Monomer reference surface.
            p_length (float): Persistence length (bending
                stiffness) in Angstroms; must be > 0.
            hp_length (float): Helix period length: distance
                along the chain axis for one full 360° azimuthal
                rotation; must be > 0.
            mz_length (float): Axial rise per monomer; must be
                > 0.
            z_length_f (float): Z-axis elevation factor
                (slope); must be >= 0 (default 0).
            p0 (array-like): Seed 3-D coordinate
                (default origin).
            vz (array-like): Reference z-axis direction
                (default (0, 0, 1)).
            rot_rand (bool): If True (default), the seed
                orientation is drawn randomly; otherwise aligned
                to vz.

        Raises:
            ValueError: If any length parameter is out of range,
                or vz does not have exactly 3 elements.
        """
        super(HelixFiber, self).__init__(m_surf)
        if (
            l_length <= 0
            or p_length <= 0
            or z_length_f < 0
            or hp_length <= 0
            or mz_length <= 0
        ):
            raise ValueError(
                "l_length, p_length, hp_length, "
                "mz_length must be > 0 and "
                "z_length_f must be >= 0."
            )
        self.__l, self.__lp, self.__lz = l_length, p_length, z_length_f
        self.__hp, self.__mz_length = hp_length, mz_length
        self.__hp_astep = (360.0 * self.__mz_length) / self.__hp
        self.__compute_helical_parameters()
        if not hasattr(vz, "__len__") or len(vz) != 3:
            raise ValueError("vz must have exactly 3 elements.")
        self.__vz = np.asarray(vz, dtype=float)
        # Curve state member variables
        self.__ct, self.__za, self.__rq = (
            0.0,
            0.0,
            np.asarray((1.0, 0.0, 0.0, 0.0)),
        )  # z-aligned curve time (considering speed 1)
        self.set_reference(np.asarray(p0), self.__vz, rot_rand=rot_rand)

    def set_reference(
        self, p0=(0.0, 0.0, 0.0), vz=(0.0, 0.0, 1.0), rot_rand=True
    ):
        """Reset the fiber to a new seed position.

        Discards any previously placed monomers and reinitialises
        the internal helix state (curve time, azimuthal angle,
        reference quaternion).

        Args:
            p0 (array-like): Seed 3-D coordinate
                (default origin).
            vz (array-like): Reference z-axis direction
                (default (0, 0, 1)).
            rot_rand (bool): If True (default), seed orientation
                is random; otherwise aligned to vz.

        Raises:
            ValueError: If p0 does not have exactly 3 elements.
        """
        if not hasattr(p0, "__len__") or len(p0) != 3:
            raise ValueError("p0 must have exactly 3 elements.")
        self._p = np.asarray(p0)
        if rot_rand:
            t = gen_uni_s2_sample(np.asarray((0.0, 0.0, 0.0)), 1.0)
            M = vect_to_zmat(t, mode="passive")
            self.__rq = rot_to_quat(M)
        else:
            self.__rq = np.asarray((1.0, 0.0, 0.0, 0.0))
        t = self.__compute_tangent(self.__ct)
        t = t * (self.__mz_length / vector_module(t))
        self.__ct += self.__l
        q1 = angle_axis_to_quat(self.__za, t[0], t[1], t[2])
        M = vect_to_zmat(t, mode="passive")
        q = rot_to_quat(M)
        hold_q = quat_mult(q, q1)
        hold_monomer = Monomer(self._m_surf, self._m_diam)
        hold_monomer.rotate_q(hold_q)
        hold_monomer.translate(p0)
        self.add_monomer(p0, t, hold_q, hold_monomer)

    def gen_new_monomer(
        self,
        over_tolerance=0,
        voi=None,
        v_size=1,
        net=None,
        branch=None,
        max_dist=None,
    ):
        """Advance the helix by one monomer step.

        Computes the next position and orientation from the
        Frenet-Serret helix parametrisation.  The attempt fails
        if it would overlap an existing monomer, leave the VOI,
        or collide with other polymers in net.

        Args:
            over_tolerance (float): Allowed overlap fraction
                (default 0).
            voi (numpy.ndarray, optional): Boolean VOI mask.
            v_size (float): VOI voxel size in Angstroms
                (default 1).
            net (Network, optional): Polymer network to avoid;
                None disables inter-polymer collision checking.
            branch (object, optional): Branch reference used to
                skip network collision near branching points.
            max_dist (float, optional): Maximum inter-centre
                distance for network collision; defaults to 1.2x
                monomer diameter.

        Returns:
            tuple | None: A 4-tuple (r, t, q, monomer) on success,
                or None if placement failed.
        """

        hold_m = Monomer(self._m_surf, self._m_diam)

        t = self.__compute_tangent(self.__ct)
        t = t * (self.__mz_length / vector_module(t))
        self.__za = wrap_angle(self.__za + self.__hp_astep)
        q1 = angle_axis_to_quat(self.__za, t[0], t[1], t[2])
        M = vect_to_zmat(t, mode="passive")
        q = rot_to_quat(M)
        hold_m.rotate_q(quat_mult(q, q1))

        hold_r = self._r[-1]
        self.__ct += self.__l
        r = hold_r + t
        hold_m.translate(r)

        if voi is not None:
            if hold_m.overlap_voi(voi, v_size, over_tolerance=over_tolerance):
                return None
        if branch is None:
            if self.overlap_polymer(hold_m, over_tolerance=over_tolerance):
                return None
            if net is not None:
                if hold_m.overlap_net(
                    net, over_tolerance=over_tolerance, max_dist=max_dist
                ):
                    return None
        else:
            branch_dst = points_distance(branch.point, hold_m.center_mass)
            if branch_dst > hold_m.diameter:
                if self.overlap_polymer(hold_m, over_tolerance=over_tolerance):
                    return None
                if net is not None:
                    if hold_m.overlap_net(net, over_tolerance=over_tolerance):
                        return None

        return r, t, q, hold_m

    def __compute_helical_parameters(self):
        """Compute helix curvature and elevation parameters.

        Derives the circular radius (a) from the persistence
        length via k = arccos(exp(-l/lp))/l, and the Z-axis
        elevation coefficient (b) from z_length_f.  Updates
        self.__a and self.__b in place.
        """

        # Compute curvature from persistence
        k = math.acos(math.exp(-self.__l / self.__lp)) / self.__l

        # Compute circular parameter, a, from curvature
        self.__a = 1 / k

        # Compute Z-axis elevation from the circular parameter
        self.__b = self.__lz * self.__a

    def __compute_tangent(self, t):
        """Compute the z-aligned helix tangent vector at curve time t.

        Uses the Frenet-Serret parametrisation of a helix with
        radius a and elevation b, then rotates the result by the
        reference quaternion __rq.

        Args:
            t (float): Curve time parameter (assuming unit speed).

        Returns:
            numpy.ndarray: Normalised tangent vector (3 elements).
        """
        sq = math.sqrt(self.__a * self.__a + self.__b * self.__b)
        s = t
        s_sq = s / sq
        t = (1.0 / sq) * np.asarray(
            (-self.__a * math.sin(s_sq), self.__a * math.cos(s_sq), self.__b)
        )
        return rot_vect_quat(t, self.__rq)
