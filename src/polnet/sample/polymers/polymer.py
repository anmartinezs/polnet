"""Abstract polymer chain composed of ordered monomers.

Defines :class:`Polymer`, the abstract base class for all linear
chain models (helical fibers, SAWLC chains).  Subclasses implement
:meth:`set_reference` and :meth:`gen_new_monomer` to drive chain
growth.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

from abc import (
    ABC,
    abstractmethod,
)

import numpy as np
import vtk

from .monomer import Monomer
from ...utils.poly import (
    GTRUTH_VTP_LBLS,
    add_label_to_poly,
    merge_polys,
    points_to_poly_spheres,
    poly_diam,
)
from ...utils.utils import (
    VTK_RAY_TOLERANCE,
    points_distance,
)


class Polymer(ABC):
    """Abstract class for modeling a Polymer (a sequence of monomers).

    Subclasses must implement :meth:`set_reference` and :meth:`gen_new_monomer`.

    Protected attributes (available to subclasses via ``self._name``):
        _m_surf:   monomer reference surface (vtkPolyData).
        _m_diam:   monomer diameter derived from the reference surface.
        _m:        list of placed Monomer instances.
        _p:        starting point set by ``set_reference`` (ndarray or None).
        _t:        list of tangent vectors, one per monomer.
        _r:        list of centre-of-mass position vectors, one per monomer.
        _q:        list of rotation quaternions (w,x,y,z), one per monomer.
        _t_length: total length of the polymer (float).
        _ids:      list of integer IDs, one per monomer.
        _codes:    list of code strings, one per monomer.
    """

    def __init__(self, m_surf, id0=0, code0=""):
        """Initialise a Polymer from a monomer surface template.

        Args:
            m_surf (vtk.vtkPolyData): Monomer reference surface;
                its diameter is used to initialise _m_diam.
            id0 (int): ID for the seed monomer (default 0).
            code0 (str): Code string for the seed monomer
                (default '').

        Raises:
            ValueError: If the computed monomer diameter is <= 0.
        """
        h_diam = poly_diam(m_surf)
        if h_diam <= 0:
            raise ValueError("Monomer diameter must be > 0.")
        self._m_surf = m_surf
        self._m_diam = h_diam
        self._m = list()
        self._p = None
        self._t, self._r, self._q = list(), list(), list()
        self._t_length = 0
        self._ids = list()
        self._codes = list()

    @property
    def num_monomers(self):
        """Return the number of monomers in the polymer."""
        return len(self._m)

    @property
    def num_mmers(self):
        """Alias for :attr:`num_monomers`."""
        return self.num_monomers

    @property
    def mmer_ids(self):
        """Return the list of monomer IDs for this polymer.

        Returns:
            list[int]: Monomer ID values in insertion order.
        """
        return self._ids

    def get_mmer_id(self, m_id):
        """Return the integer ID of the monomer at position m_id.

        Args:
            m_id (int): Zero-based position in the polymer chain.

        Returns:
            int: The monomer ID at that position.
        """
        return self._ids[m_id]

    def get_mmer_code(self, m_id):
        """Return the code string of the monomer at position m_id.

        Args:
            m_id (int): Zero-based position in the polymer chain.

        Returns:
            str: The monomer code string.
        """
        return self._codes[m_id]

    @property
    def vol(self):
        """Compute the total enclosed volume of all monomers.

        Returns:
            float: Sum of individual monomer volumes; 0 if the
                polymer is empty.
        """
        vol = 0
        if len(self._m) == 0:
            return vol
        else:
            for m in self._m:
                vol += m.vol
            return vol

    def get_area(self, mode="sphere"):
        """Compute the total projected area of all monomers.

        Args:
            mode (str): Approximation mode.  Currently only
                'sphere' is supported (each monomer is treated as a
                sphere). Defaults to 'sphere'.

        Returns:
            float: Sum of monomer projected areas; 0 if empty.
        """
        area = 0
        if len(self._m) == 0:
            return area
        else:
            for m in self._m:
                area += m.area
            return area

    @property
    def total_len(self):
        """Return the total accumulated chain length in Angstroms."""
        return self._t_length

    def get_monomer(self, m_id):
        """Return the Monomer at position m_id.

        Args:
            m_id (int): Zero-based index in [0, get_num_monomers()-1].

        Returns:
            Monomer: The monomer instance at that position.
        """
        return self._m[m_id]

    @property
    def mmers_list(self):
        """Return all monomers in the polymer chain.

        Returns:
            list[Monomer]: Monomers in insertion order.
        """
        return self._m

    def get_mmer_center(self, m_id):
        """Return the centre-of-mass coordinates of monomer m_id.

        Args:
            m_id (int): Zero-based monomer index.

        Returns:
            numpy.ndarray: 3-D coordinate of the monomer centre.
        """
        return self._r[m_id]

    def get_mmer_rotation(self, m_id):
        """Return the rotation quaternion of monomer m_id.

        Args:
            m_id (int): Zero-based monomer index.

        Returns:
            numpy.ndarray: Unit quaternion [w, x, y, z].
        """
        return self._q[m_id]

    @property
    def tail_point(self):
        """Return the centre coordinate of the last monomer added.

        Returns:
            numpy.ndarray: 3-D tail coordinate.
        """
        return self._r[-1]

    @property
    def vtp(self):
        """Build a labelled vtkPolyData merging all monomer surfaces.

        Each monomer surface is annotated with its integer ID via
        :func:`add_label_to_poly` before merging.

        Returns:
            vtk.vtkPolyData: Combined polygon dataset.
        """

        app_flt = vtk.vtkAppendPolyData()

        for m_id in range(len(self._m)):
            m_poly = self._m[m_id].vtp
            add_label_to_poly(m_poly, self._ids[m_id], GTRUTH_VTP_LBLS)
            app_flt.AddInputData(m_poly)
        app_flt.Update()

        return app_flt.GetOutput()

    def get_skel(self, add_verts=True, add_lines=True, verts_rad=0):
        """Build a skeleton vtkPolyData for the polymer chain.

        Represents each monomer as a vertex (or sphere) and
        connects consecutive centres with line segments.

        Args:
            add_verts (bool): Include vertex glyphs (default True).
            add_lines (bool): Include connecting lines
                (default True).
            verts_rad (float): Sphere radius for vertex glyphs;
                use <= 0 for point vertices (default 0).

        Returns:
            vtk.vtkPolyData: Skeleton polygon dataset.
        """

        poly, points = vtk.vtkPolyData(), vtk.vtkPoints()
        verts, lines = vtk.vtkCellArray(), vtk.vtkCellArray()
        sph_points = list()

        if len(self._r) == 1:
            sph_points.append(self._r[0])
            id_p0 = points.InsertNextPoint(self._r[0])
            if add_verts and (verts_rad <= 0):
                verts.InsertNextCell(1)
                verts.InsertCellPoint(id_p0)
        else:
            # Insert first point for sphere glyphs and lines
            sph_points.append(self._r[0])
            prev_id = points.InsertNextPoint(self._r[0])
            if add_verts and (verts_rad <= 0):
                verts.InsertNextCell(1)
                verts.InsertCellPoint(prev_id)
            # Connect subsequent points
            for i in range(1, len(self._r)):
                sph_points.append(self._r[i])
                curr_id = points.InsertNextPoint(self._r[i])
                if add_verts and (verts_rad <= 0):
                    verts.InsertNextCell(1)
                    verts.InsertCellPoint(curr_id)
                if add_lines:
                    lines.InsertNextCell(2)
                    lines.InsertCellPoint(prev_id)
                    lines.InsertCellPoint(curr_id)
                prev_id = curr_id

        # Construct poly
        poly.SetPoints(points)
        if add_verts and (verts_rad <= 0):
            poly.SetVerts(verts)
        if add_lines:
            poly.SetLines(lines)

        # Spheres case
        if add_verts and (verts_rad > 0):
            sph_vtp = points_to_poly_spheres(sph_points, verts_rad)
            poly = merge_polys(sph_vtp, poly)

        return poly

    def add_monomer(self, r, t, q, m, mmer_id=0, code=""):
        """Append a pre-transformed monomer to the polymer chain.

        Updates the internal lists of positions, tangents,
        quaternions, and IDs, and recalculates the total chain
        length.

        Args:
            r (numpy.ndarray): Centre-of-mass position (x, y, z).
            t (numpy.ndarray): Tangent direction (x, y, z).
            q (numpy.ndarray): Rotation quaternion [w, x, y, z].
            m (Monomer): The placed Monomer instance.
            mmer_id (int): Monomer type ID (default 0).
            code (str): Monomer code string (default '').

        Raises:
            TypeError: If r, t, q are not numpy arrays or m is not
                a Monomer.
            ValueError: If r, t are not length-3 or q is not
                length-4.
        """
        if not (
            isinstance(r, np.ndarray)
            and isinstance(t, np.ndarray)
            and isinstance(q, np.ndarray)
            and isinstance(m, Monomer)
        ):
            raise TypeError(
                "r, t, q must be numpy arrays and " "m must be a Monomer."
            )
        if len(r) != 3 or len(t) != 3 or len(q) != 4:
            raise ValueError("r and t must have 3 elements, " "q must have 4.")
        self._r.append(r)
        self._t.append(t)
        self._q.append(q)
        self._m.append(m)
        self._ids.append(mmer_id)
        self._codes.append(code)
        if self.num_monomers <= 1:
            self._t_length = 0
        else:
            self._t_length += points_distance(self._r[-1], self._r[-2])

    def insert_density_svol(
        self, m_svol, tomo, v_size=1, merge="max", off_svol=None
    ):
        """Stamp all monomers of this polymer into a tomogram.

        Args:
            m_svol (numpy.ndarray | list[numpy.ndarray]): Reference
                monomer sub-volume(s).  A single array is wrapped in
                a list; list index corresponds to monomer ID.
            tomo (numpy.ndarray): Target tomogram modified in place.
            v_size (float): Voxel size in Angstroms (default 1).
            merge (str): Blending strategy: 'max' (default), 'min',
                'sum', or 'insert'.
            off_svol (array-like, optional): Offset for sub-volume
                centre coordinates.
        """
        if isinstance(m_svol, np.ndarray):
            m_svol = [
                m_svol,
            ]
        for mmer, mmer_id in zip(self._m, self._ids):
            mmer.insert_density_svol(
                m_svol[mmer_id], tomo, v_size, merge=merge, off_svol=off_svol
            )

    @abstractmethod
    def set_reference(self):
        """Set the reference monomer for this polymer."""
        raise NotImplementedError

    @abstractmethod
    def gen_new_monomer(self):
        """Generate a new monomer and return its placement data."""
        raise NotImplementedError

    def overlap_polymer(self, monomer, over_tolerance=0):
        """Test whether a monomer overlaps any monomer in this polymer.

        Iterates backwards over the chain (most recent monomers first)
        and short-circuits once the inter-centre distance exceeds
        the test monomer's diameter.

        Args:
            monomer (Monomer): The candidate monomer to test.
            over_tolerance (float): Maximum allowed overlap fraction
                (default 0).

        Returns:
            bool: True if an overlap exceeding over_tolerance is
                found.
        """

        selector = vtk.vtkSelectEnclosedPoints()
        selector.SetTolerance(VTK_RAY_TOLERANCE)
        selector.Initialize(monomer.vtp)

        # Polymer loop, no need to process monomer beyond diameter distance
        diam, center = monomer.diameter, monomer.center_mass
        for i in range(len(self._m) - 1, -1, -1):
            hold_monomer = self._m[i]
            dist = points_distance(center, hold_monomer.center_mass)
            if dist <= diam:
                poly_b = hold_monomer.vtp
                count, n_points = 0.0, poly_b.GetNumberOfPoints()
                n_points_if = 1.0 / float(n_points)
                for j in range(n_points):
                    if selector.IsInsideSurface(poly_b.GetPoint(j)) > 0:
                        count += 1
                        over = count * n_points_if
                        if over > over_tolerance:
                            return True

        return False
