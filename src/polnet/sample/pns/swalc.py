"""Self-Avoiding Worm-Like Chain (SAWLC) polymer models.

Defines :class:`SAWLC` (chain in free 3-D space) and
:class:`SAWLCPoly` (chain constrained to a VTK PolyData surface),
the core stochastic polymer models for placing cytosolic and
membrane-bound protein complexes.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import numpy as np
import vtk

from ..polymers import (
    Monomer,
    Polymer,
)
from ...utils.affine import gen_rand_unit_quaternion
from ...utils.poly import (
    find_point_on_poly,
    gen_rand_quaternion_on_vector,
    gen_uni_s2_sample_on_poly,
)
from ...utils.utils import gen_uni_s2_sample


class SAWLC(Polymer):
    """Self-Avoiding Worm-Like Chain polymer in free 3-D space.

    Each monomer is placed at a uniformly random direction
    (on the unit sphere) at the fixed link length from the
    previous monomer, subject to self-avoidance and optional
    VOI constraints.
    """

    def __init__(
        self, l_length, m_surf, p0=(0, 0, 0), id0=0, code0="", rot=None
    ):
        """Initialise a SAWLC chain.

        Args:
            l_length (float): Link length in Angstroms; must be > 0.
            m_surf (vtk.vtkPolyData): Monomer reference surface.
            p0 (array-like): Seed 3-D coordinate (default origin).
            id0 (int): Type ID of the seed monomer (default 0).
            code0 (str): Code string for the seed monomer
                (default '').
            rot (array-like, optional): External rotation
                quaternion [w, x, y, z] for the seed monomer;
                None for random (default).

        Raises:
            ValueError: If l_length <= 0.
        """
        super(SAWLC, self).__init__(m_surf)
        if l_length <= 0:
            raise ValueError("l_length must be > 0.")
        self.__l = l_length
        self.set_reference(p0, id0=id0, code0=code0, rot=rot)

    def set_reference(self, p0=(0.0, 0.0, 0), id0=0, code0="", rot=None):
        """Reset the chain to a new seed point.

        Discards any previously added monomers and initialises the
        chain with a single monomer at p0.

        Args:
            p0 (array-like): Seed 3-D coordinate (default origin).
            id0 (int): Type ID for the seed monomer (default 0).
            code0 (str): Code string for the seed monomer
                (default '').
            rot (array-like, optional): Rotation quaternion
                [w, x, y, z]; None for random (default).

        Raises:
            ValueError: If p0 does not have exactly 3 elements, or
                rot does not have exactly 4 elements.
        """
        if not hasattr(p0, "__len__") or len(p0) != 3:
            raise ValueError("p0 must have exactly 3 elements.")
        self._p = np.asarray(p0)
        hold_monomer = Monomer(self._m_surf, self._m_diam)
        if rot is None:
            hold_q = gen_rand_unit_quaternion()
        else:
            if not hasattr(rot, "__len__") or len(rot) != 4:
                raise ValueError("rot must have exactly 4 elements.")
            hold_q = rot
        hold_monomer.rotate_q(hold_q)
        hold_monomer.translate(p0)
        self.add_monomer(
            p0,
            np.asarray((0.0, 0.0, 0.0)),
            hold_q,
            hold_monomer,
            mmer_id=id0,
            code=code0,
        )

    def gen_new_monomer(
        self,
        over_tolerance=0,
        voi=None,
        v_size=1,
        fix_dst=None,
        ext_surf=None,
        rot=None,
    ):
        """Try to place the next monomer in free 3-D space.

        Samples a random point within a hollow shell around the
        current polymer tail and applies a random rotation.

        Args:
            over_tolerance (float): Allowed overlap fraction
                (default 0).
            voi (numpy.ndarray, optional): Boolean VOI mask;
                None disables the VOI check (default).
            v_size (float): VOI voxel size in Angstroms
                (default 1).
            fix_dst (float, optional): Override the link length.
            ext_surf (vtk.vtkPolyData, optional): Alternative
                monomer surface.
            rot (array-like, optional): External rotation
                quaternion [w, x, y, z].

        Returns:
            tuple | None: A 4-tuple (r, t, q, monomer) on success,
                or None if placement failed.
        """

        if fix_dst is None:
            hold_l = self.__l
        else:
            hold_l = fix_dst
        t = gen_uni_s2_sample(np.asarray((0.0, 0.0, 0.0)), hold_l)
        r = self._r[-1] + t

        if rot is None:
            q = gen_rand_unit_quaternion()
        else:
            q = rot

        if ext_surf is None:
            hold_m = Monomer(self._m_surf, self._m_diam)
        else:
            hold_m = Monomer(ext_surf, self._m_diam)
        hold_m.rotate_q(q)
        hold_m.translate(r)

        if self.overlap_polymer(hold_m, over_tolerance=over_tolerance):
            return None
        elif voi is not None:
            if hold_m.overlap_voi(voi, v_size, over_tolerance=over_tolerance):
                return None

        return r, t, q, hold_m


class SAWLCPoly(Polymer):
    """SAWLC polymer constrained to a vtkPolyData surface.

    Each monomer is placed on the membrane surface by sampling
    within a hollow sphere around the previous monomer's position,
    then snapping to the nearest surface point.
    """

    def __init__(
        self, poly, l_length, m_surf, p0=(0, 0, 0), id0=0, code="", rot=None
    ):
        """Initialise a surface-constrained SAWLC chain.

        Args:
            poly (vtk.vtkPolyData): Reference membrane surface.
            l_length (float): Link length in Angstroms; must be > 0.
            m_surf (vtk.vtkPolyData): Monomer reference surface.
            p0 (array-like): Seed 3-D coordinate (default origin).
            id0 (int): Type ID for the seed monomer (default 0).
            code (str): Code string for the seed monomer
                (default '').
            rot (array-like, optional): External rotation
                quaternion [w, x, y, z]; None for random
                (default).

        Raises:
            TypeError: If poly is not a vtkPolyData.
            ValueError: If l_length <= 0.
        """
        super(SAWLCPoly, self).__init__(m_surf)
        if not isinstance(poly, vtk.vtkPolyData):
            raise TypeError("poly must be a vtkPolyData.")
        if l_length <= 0:
            raise ValueError("l_length must be > 0.")
        self.__l = l_length
        self.__poly = poly
        self.set_reference(p0, id0=id0, code0=code, rot=rot)

    def set_reference(self, p0=(0.0, 0.0, 0), id0=0, code0="", rot=None):
        """Reset the surface chain to a new seed point.

        The seed is snapped to the nearest point on the membrane
        surface; the monomer orientation is aligned with the local
        surface normal.

        Args:
            p0 (array-like): Requested seed coordinate; snapped to
                the surface.
            id0 (int): Type ID for the seed monomer (default 0).
            code0 (str): Code string for the seed monomer
                (default '').
            rot (array-like, optional): External rotation
                quaternion [w, x, y, z]; None uses the surface
                normal (default).

        Raises:
            ValueError: If p0 does not have exactly 3 elements, or
                rot does not have exactly 4 elements.
        """
        if not hasattr(p0, "__len__") or len(p0) != 3:
            raise ValueError("p0 must have exactly 3 elements.")
        self._p, hold_n = find_point_on_poly(np.asarray(p0), self.__poly)
        hold_monomer = Monomer(self._m_surf, self._m_diam)
        if rot is None:
            hold_q = gen_rand_quaternion_on_vector(hold_n)
        else:
            if not hasattr(rot, "__len__") or len(rot) != 4:
                raise ValueError("rot must have exactly 4 elements.")
            hold_q = rot
        hold_monomer.rotate_q(hold_q)
        hold_monomer.translate(self._p)
        self.add_monomer(
            self._p,
            np.asarray((0.0, 0.0, 0.0)),
            hold_q,
            hold_monomer,
            mmer_id=id0,
            code=code0,
        )

    def gen_new_monomer(
        self,
        over_tolerance=0,
        voi=None,
        v_size=1,
        fix_dst=None,
        ext_surf=None,
        rot=None,
    ):
        """Try to place the next surface-constrained monomer.

        Samples a point on the membrane within a hollow shell
        around the current tail and orients the monomer along the
        local normal.

        Args:
            over_tolerance (float): Allowed overlap fraction
                (default 0).
            voi (numpy.ndarray, optional): Boolean VOI mask;
                None disables the VOI check (default).
            v_size (float): VOI voxel size in Angstroms
                (default 1).
            fix_dst (float, optional): Override the link length.
            ext_surf (vtk.vtkPolyData, optional): Alternative
                monomer surface.
            rot (array-like, optional): External rotation
                quaternion [w, x, y, z].

        Returns:
            tuple | None: A 4-tuple (r, t, q, monomer) on success,
                or None if placement failed.
        """

        if fix_dst is None:
            hold_l = self.__l
        else:
            hold_l = fix_dst
        r = gen_uni_s2_sample_on_poly(self._r[-1], hold_l, 2, self.__poly)
        if r is None:
            return None
        r = np.asarray(r)
        t = r - self._r[-1]

        if rot is None:
            hold_n = find_point_on_poly(r, self.__poly)[1]
            q = gen_rand_quaternion_on_vector(hold_n)
        else:
            q = rot

        if ext_surf is None:
            hold_m = Monomer(self._m_surf, self._m_diam)
        else:
            hold_m = Monomer(ext_surf, self._m_diam)
        hold_m.rotate_q(q)
        hold_m.translate(r)

        if self.overlap_polymer(hold_m, over_tolerance=over_tolerance):
            return None
        elif voi is not None:
            if hold_m.overlap_voi(voi, v_size, over_tolerance=over_tolerance):
                return None

        return r, t, q, hold_m
