"""Single polymer repeat unit with VTK surface geometry.

Defines :class:`Monomer`, which wraps a vtkPolyData surface and
provides placement (rotation + translation) of the unit inside
a 3-D tomographic volume.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import numpy as np
import vtk

from ...utils.affine import (
    poly_rotate_wxyz,
    poly_translate,
    quat_to_angle_axis,
    tomo_rotate,
    vect_rotate,
)
from ...utils.poly import (
    poly_center_mass,
    poly_diam,
    poly_volume,
)
from ...utils.utils import (
    insert_svol_tomo,
    points_distance,
    VTK_RAY_TOLERANCE,
)

MB_DOMAIN_FIELD_STR = "mb_domain"


class Monomer:
    """Single polymer repeat unit with VTK surface geometry.

    Wraps a vtkPolyData surface and tracks the ordered sequence of
    rigid transformations (rotations and translations) applied to it
    so that the matching sub-volume density can be reproduced.

    Attributes:
        __m_surf (vtk.vtkPolyData): Copy of the reference surface,
            updated after each transformation.
        __diam (float): Diameter of the monomer in Angstroms.
        __rot_angs (numpy.ndarray): Euler rotation angles in degrees
            (kept for legacy/debug use).
        __bounds (numpy.ndarray): Axis-aligned bounding box in
            format [xmin, xmax, ymin, ymax, zmin, zmax].
        __trans (list[tuple]): Ordered transformation queue.  Each
            element is a 2-tuple (type, value) where type is 'r'
            (rotation quaternion) or 't' (translation vector).
    """

    def __init__(self, m_surf, diam):
        """Initialise a Monomer from a VTK surface and diameter.

        Args:
            m_surf (vtk.vtkPolyData): Monomer reference surface.
            diam (float): Monomer diameter in Angstroms; must be
                > 0.

        Raises:
            TypeError: If m_surf is not a vtkPolyData instance.
            ValueError: If diam is not positive.
        """
        if not isinstance(m_surf, vtk.vtkPolyData):
            raise TypeError("m_surf must be a vtkPolyData instance.")
        if diam <= 0:
            raise ValueError("diam must be greater than 0.")
        self.__m_surf = vtk.vtkPolyData()
        self.__m_surf.DeepCopy(m_surf)
        self.__diam = diam
        self.__rot_angs = np.asarray((0.0, 0.0, 0.0))
        self.__bounds = np.zeros(shape=6)
        self.compute_bounds()
        # Ordered transformation queue, each entry is a 2-tuple
        # (str in ['t', 'r'], transform value in [vector, quaternion])
        self.__trans = list()

    @property
    def vtp(self):
        """Return the current vtkPolyData surface of the monomer."""
        return self.__m_surf

    @property
    def center_mass(self):
        """Compute and return the monomer centre of mass.

        Returns:
            numpy.ndarray: Centre of mass as a 3-element array.
        """
        return np.asarray(poly_center_mass(self.__m_surf))

    @property
    def diameter(self):
        """Return the monomer diameter in Angstroms."""
        return self.__diam

    @property
    def trans_list(self):
        """Return the ordered transformation queue.

        Returns:
            list[tuple]: Each element is a 2-tuple (type, value)
                where type is 'r' (rotation quaternion [w,x,y,z])
                or 't' (translation vector [x,y,z]).
        """
        return self.__trans

    def compute_bounds(self):
        """Recompute the axis-aligned bounding box from current geometry."""
        arr = self.__m_surf.GetPoints().GetData()
        self.__bounds[0], self.__bounds[1] = arr.GetRange(0)
        self.__bounds[2], self.__bounds[3] = arr.GetRange(1)
        self.__bounds[4], self.__bounds[5] = arr.GetRange(2)

    def rotate_q(self, q):
        """Apply a rigid rotation about the origin using a quaternion.

        The surface geometry and bounding box are updated in place,
        and the transform is appended to the transformation queue.

        Args:
            q (array-like): Unit quaternion as [w, x, y, z].
        """
        w, v_axis = quat_to_angle_axis(q[0], q[1], q[2], q[3])
        self.__m_surf = poly_rotate_wxyz(
            self.__m_surf, w, v_axis[0], v_axis[1], v_axis[2]
        )
        self.compute_bounds()
        self.__trans.append(("r", q))

    def translate(self, t_v):
        """Apply a rigid translation to the monomer surface.

        The surface geometry and bounding box are updated in place,
        and the transform is appended to the transformation queue.

        Args:
            t_v (array-like): Translation vector (x, y, z) in
                Angstroms.
        """
        self.__m_surf = poly_translate(self.__m_surf, t_v)
        self.compute_bounds()
        self.__trans.append(("t", t_v))

    def point_in_bounds(self, point):
        """Test whether a 3-D point lies inside the monomer's AABB.

        Args:
            point (array-like): 3-element coordinate to test.

        Returns:
            bool: True if the point is inside the bounding box.
        """
        x_over, y_over, z_over = True, True, True
        if (self.__bounds[0] > point[0]) or (self.__bounds[1] < point[0]):
            x_over = False
        if (self.__bounds[2] > point[1]) or (self.__bounds[3] < point[1]):
            y_over = False
        if (self.__bounds[4] > point[2]) or (self.__bounds[5] < point[2]):
            z_over = False
        return x_over and y_over and z_over

    def bound_in_bounds(self, bounds):
        """Test whether this monomer's AABB overlaps another AABB.

        Args:
            bounds (array-like): Reference bounding box in the form
                [xmin, xmax, ymin, ymax, zmin, zmax].

        Returns:
            bool: True if the two bounding boxes overlap on all
                three axes.
        """
        x_over, y_over, z_over = True, True, True
        if (self.__bounds[0] > bounds[1]) or (self.__bounds[1] < bounds[0]):
            x_over = False
        if (self.__bounds[2] > bounds[3]) or (self.__bounds[3] < bounds[2]):
            y_over = False
        if (self.__bounds[4] > bounds[5]) or (self.__bounds[5] < bounds[4]):
            z_over = False
        return x_over and y_over and z_over

    def overlap_voi(self, voi, v_size=1, over_tolerance=0):
        """Test whether this monomer exceeds the VOI boundary.

        Returns True when the fraction of surface vertices outside
        the VOI exceeds over_tolerance.  Vertices belonging to
        membrane domains (label > 0 in MB_DOMAIN_FIELD_STR) are
        excluded from the check.

        Args:
            voi (numpy.ndarray): Boolean 3-D VOI mask (True inside).
            v_size (float): Voxel size in Angstroms (default 1).
            over_tolerance (float): Maximum outside fraction before
                overlap is declared (default 0).

        Returns:
            bool: True if the monomer is outside the VOI beyond the
                tolerance, False otherwise.

        Raises:
            ValueError: If v_size <= 0.
            TypeError: If voi is not a boolean numpy array.
        """

        if v_size <= 0:
            raise ValueError("v_size must be greater than 0.")
        if not isinstance(voi, np.ndarray) or voi.dtype != "bool":
            raise TypeError("voi must be a boolean numpy array.")
        nx, ny, nz = voi.shape
        v_size_i = 1.0 / v_size
        mbd_prop = self.__m_surf.GetPointData().GetArray(MB_DOMAIN_FIELD_STR)

        count, n_points = 0.0, self.__m_surf.GetNumberOfPoints()
        n_points_if = 1.0 / float(n_points)
        if mbd_prop is None:
            for i in range(self.__m_surf.GetNumberOfPoints()):
                pt = np.asarray(self.__m_surf.GetPoint(i)) * v_size_i
                x, y, z = np.round(pt).astype(int)
                if (
                    (x < nx)
                    and (y < ny)
                    and (z < nz)
                    and (x >= 0)
                    and (y >= 0)
                    and (z >= 0)
                ):
                    if not voi[x, y, z]:
                        count += 1
                        over = count * n_points_if
                        if over > over_tolerance:
                            return True
                else:
                    count += 1
                    over = count * n_points_if
                    if over > over_tolerance:
                        return True
        else:
            for i in range(self.__m_surf.GetNumberOfPoints()):
                if mbd_prop.GetValue(i) == 0:
                    pt = np.asarray(self.__m_surf.GetPoint(i)) * v_size_i
                    x, y, z = np.round(pt).astype(int)
                    if (
                        (x < nx)
                        and (y < ny)
                        and (z < nz)
                        and (x >= 0)
                        and (y >= 0)
                        and (z >= 0)
                    ):
                        if not voi[x, y, z]:
                            count += 1
                            over = count * n_points_if
                            if over > over_tolerance:
                                return True
                    else:
                        count += 1
                        over = count * n_points_if
                        if over > over_tolerance:
                            return True

        return False

    @property
    def vol(self):
        """Return the enclosed volume of the monomer surface.

        Returns:
            float: Volume in voxel units cubed.
        """
        return poly_volume(self.__m_surf)

    @property
    def area(self):
        """Estimate the projected cross-sectional area of the monomer.

        The monomer is approximated as a sphere with diameter equal
        to its poly_diam, yielding area = pi*(d/2)^2.

        Returns:
            float: Projected area in voxel units squared.
        """
        diam = poly_diam(self.__m_surf)
        return np.pi * diam * diam * 0.25

    def get_copy(self):
        """Return a deep copy of the current Monomer.

        Returns:
            Monomer: New Monomer instance with an independent copy
                of the surface geometry.
        """
        return Monomer(self.__m_surf, self.__diam)

    def insert_density_svol(
        self, m_svol, tomo, v_size=1, merge="max", off_svol=None
    ):
        """Stamp a monomer sub-volume into a tomogram.

        Replays the transformation queue on a copy of the reference
        sub-volume (rotations applied via tomo_rotate, translations
        accumulated) and then blends it into tomo.

        Args:
            m_svol (numpy.ndarray): Reference monomer density
                sub-volume.
            tomo (numpy.ndarray): Target tomogram modified in place.
            v_size (float): Voxel size in Angstroms (default 1).
            merge (str): Blending strategy: 'max' (default), 'min',
                'sum', or 'insert'.
            off_svol (array-like, optional): Additional (x, y, z)
                offset applied to the sub-volume centre after
                replaying rotations.
        """
        v_size_i = 1.0 / v_size
        tot_v = np.asarray((0.0, 0.0, 0.0))
        hold_svol = m_svol
        for trans in self.__trans:
            if trans[0] == "t":
                tot_v += trans[1] * v_size_i
            elif trans[0] == "r":
                if merge == "min":
                    if hold_svol.dtype == bool:
                        hold_svol = tomo_rotate(
                            hold_svol,
                            trans[1],
                            order=0,
                            mode="constant",
                            cval=hold_svol.max(),
                        )
                    else:
                        hold_svol = tomo_rotate(
                            hold_svol,
                            trans[1],
                            mode="constant",
                            cval=hold_svol.max(),
                        )
                else:
                    if hold_svol.dtype == bool:
                        hold_svol = tomo_rotate(
                            hold_svol,
                            trans[1],
                            order=0,
                            mode="constant",
                            cval=hold_svol.min(),
                        )
                    else:
                        hold_svol = tomo_rotate(
                            hold_svol,
                            trans[1],
                            mode="constant",
                            cval=hold_svol.min(),
                        )
                if off_svol is not None:
                    off_svol = vect_rotate(off_svol, trans[1])
        if off_svol is not None:
            tot_v += off_svol
        insert_svol_tomo(hold_svol, tomo, tot_v, merge=merge)

    def overlap_mmer(self, mmer, over_tolerance=0):
        """Test whether this monomer overlaps another monomer.

        Uses vtkSelectEnclosedPoints to count how many vertices of
        mmer lie inside this monomer's closed surface.

        Args:
            mmer (Monomer): The monomer to test against self.
            over_tolerance (float): Maximum allowed overlap
                fraction (default 0).

        Returns:
            bool: True if the overlap fraction exceeds
                over_tolerance.
        """
        selector = vtk.vtkSelectEnclosedPoints()
        selector.SetTolerance(VTK_RAY_TOLERANCE)
        selector.Initialize(self.vtp)

        dist = points_distance(self.center_mass, mmer.center_mass)
        if dist <= self.diameter:
            poly_b = mmer.vtp
            count, n_points = 0.0, poly_b.GetNumberOfPoints()
            n_points_if = 1.0 / float(n_points)
            for i in range(n_points):
                if selector.IsInsideSurface(poly_b.GetPoint(i)) > 0:
                    count += 1
                    over = count * n_points_if
                    if over > over_tolerance:
                        return True

        return False

    def overlap_net(self, net, over_tolerance=0, max_dist=None):
        """Test whether this monomer overlaps any monomer in a network.

        Iterates over all polymers and monomers in net and performs
        the same point-in-surface test as :meth:`overlap_mmer`.

        Args:
            net (Network): Polymer network to check against.
            over_tolerance (float): Maximum allowed overlap
                fraction (default 0).
            max_dist (float, optional): Search radius in Angstroms.
                Defaults to 1.2 × monomer diameter.

        Returns:
            bool: True if any overlap beyond over_tolerance is
                found.
        """
        selector = vtk.vtkSelectEnclosedPoints()
        selector.SetTolerance(VTK_RAY_TOLERANCE)
        selector.Initialize(self.vtp)

        for pmer in net.pmers_list:
            for mmer in pmer.mmers_list:
                dist = points_distance(self.center_mass, mmer.center_mass)
                if max_dist is None:
                    max_dist_h = self.diameter * 1.2
                else:
                    max_dist_h = max_dist
                if dist <= max_dist_h:
                    poly_b = mmer.vtp
                    count, n_points = 0.0, poly_b.GetNumberOfPoints()
                    n_points_if = 1.0 / float(n_points)
                    for i in range(n_points):
                        if selector.IsInsideSurface(poly_b.GetPoint(i)) > 0:
                            count += 1
                            over = count * n_points_if
                            if over > over_tolerance:
                                return True

        return False
