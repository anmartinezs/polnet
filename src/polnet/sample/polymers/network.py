"""Abstract network of polymer chains inside a 3-D VOI.

Defines :class:`Network`, the abstract base class for all polymer
network types.  Subclasses implement :meth:`build_network` to
distribute chains within a volumetric occupancy limit.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

from abc import (
    ABC,
    abstractmethod,
)

import numpy as np
import vtk
from scipy.ndimage import binary_dilation as sp_binary_dilation

from ...utils.affine import tomo_rotate
from ...utils.utils import insert_svol_tomo

NET_TYPE_STR = "net_type"


class Network(ABC):
    """Abstract base class for a network of polymers.

    Subclasses must implement :meth:`build_network`.

    Protected attributes (available to subclasses via ``self._name``):
        _voi:         boolean 3D VOI array.
        _vol:         VOI volume in cubic Angstroms.
        _v_size:      voxel size in Angstroms.
        _pl_occ:      current polymer occupancy percentage.
        _pl:          list of Polymer instances.
        _pl_nmmers:   list of monomer counts (one per polymer).
        _svol:        monomer subvolume reference (ndarray or None).
        _min_nmmer:   minimum monomers per polymer to keep it.
        _poly_area:   total surface area for membrane-bound occupancy.
        _pmer_fails:  count of failed polymer placement attempts.
    """

    def __init__(self, voi, v_size, svol=None):
        """Initialise a polymer network inside a 3-D VOI.

        Args:
            voi (numpy.ndarray): 3-D boolean (or numeric) array
                defining the placement volume; non-zero voxels are
                treated as True.
            v_size (float): Voxel size in Angstroms.
            svol (numpy.ndarray | list | None): Reference monomer
                sub-volume(s) (default None).

        Raises:
            TypeError: If svol is provided but is not a numpy array.
        """
        self.set_voi(voi)
        # float cast avoids overflow warning on Windows for large VOIs
        self._vol = float((self._voi > 0).sum()) * v_size * v_size * v_size
        self._v_size = v_size
        self._pl_occ = 0
        self._pl = list()
        self._pl_nmmers = list()
        self._svol = svol
        if self._svol is not None:
            if not hasattr(svol, "__len__"):
                if not isinstance(self._svol, np.ndarray):
                    raise TypeError("svol must be a numpy array.")
        self._min_nmmer = 1
        self._poly_area = 0
        self._pmer_fails = 0

    def set_min_nmmer(self, min_nmmer):
        """Set the minimum monomer count for a polymer to be kept.

        Polymers shorter than this threshold are discarded during
        :meth:`build_network`.

        Args:
            min_nmmer (int): Minimum number of monomers per chain.
        """
        self._min_nmmer = int(min_nmmer)

    def set_voi(self, voi):
        """Replace the network's Volume Of Interest.

        Non-boolean inputs are binarised as voi > 0.

        Args:
            voi (numpy.ndarray): New VOI array (bool or numeric).
        """
        if voi.dtype != bool:
            self._voi = voi > 0
        else:
            self._voi = voi

    @property
    def pmer_fails(self):
        """Return the number of failed polymer placement attempts."""
        return self._pmer_fails

    @property
    def pmers_list(self):
        """Return the list of all polymer chains in the network."""
        return self._pl

    @property
    def num_pmers(self):
        """Return the number of polymer chains in the network.

        Returns:
            int: Polymer count.
        """
        return len(self._pl)

    @property
    def num_mmers(self):
        """Return the total number of monomers across all chains.

        Returns:
            int: Total monomer count.
        """
        count_mmers = 0
        for pl in self._pl:
            count_mmers += pl.num_mmers
        return count_mmers

    @property
    def polymer_occupancy(self):
        """Return the current polymer occupancy as a percentage."""
        return self._pl_occ

    @property
    def voi(self):
        """Return the current Volume Of Interest array.

        Returns:
            numpy.ndarray: Boolean 3-D VOI mask.
        """
        return self._voi

    def add_polymer(self, polymer, occ_mode="volume"):
        """Register a polymer chain and update occupancy tracking.

        Args:
            polymer (Polymer): The polymer to register.
            occ_mode (str): Occupancy accumulation basis: 'volume'
                (default, vol/total_vol × 100) or 'area' for
                membrane-bound polymers (area/total_area × 100).

        Raises:
            ValueError: If occ_mode is neither 'volume' nor 'area'.
        """
        if occ_mode not in ("volume", "area"):
            raise ValueError("occ_mode must be 'volume' or 'area'.")
        self._pl.append(polymer)
        self._pl_nmmers.append(polymer.num_mmers)
        if occ_mode == "volume":
            self._pl_occ += 100.0 * (polymer.vol / self._vol)
        # NOTE: Polymer.get_area(mode) is a method, not a property
        else:
            self._pl_occ += 100.0 * (polymer.get_area() / self._poly_area)

    @abstractmethod
    def build_network(self):
        """Populate the network with polymers up to the target occupancy.

        Subclasses must implement this method to drive chain growth
        using the specific geometric model (e.g. SAWLC fibers,
        membrane-attached chains).

        Raises:
            NotImplementedError: Always, when called on the base
                class.
        """
        raise NotImplementedError(
            "Subclasses must implement build_network method."
        )

    @property
    def vtp(self):
        """Build a labelled vtkPolyData merging all polymer surfaces.

        Returns:
            vtk.vtkPolyData: Combined polygon dataset, or an empty
                poly if the network has no polymers.
        """
        if len(self._pl) == 0:
            return vtk.vtkPolyData()
        app_flt = vtk.vtkAppendPolyData()
        for pol in self._pl:
            app_flt.AddInputData(pol.vtp)
        app_flt.Update()
        return app_flt.GetOutput()

    def get_skel(self, add_verts=True, add_lines=True, verts_rad=0):
        """Build a skeleton vtkPolyData for all polymers in the network.

        Args:
            add_verts (bool): Include vertex glyphs (default True).
            add_lines (bool): Include connecting lines (default True).
            verts_rad (float): Sphere radius for vertex glyphs (default 0).

        Returns:
            vtk.vtkPolyData: Combined skeleton dataset.
        """
        if len(self._pl) == 0:
            return vtk.vtkPolyData()
        app_flt = vtk.vtkAppendPolyData()
        for pol in self._pl:
            app_flt.AddInputData(pol.get_skel(add_verts, add_lines, verts_rad))
        app_flt.Update()
        return app_flt.GetOutput()

    def get_gtruth(self, thick=1):
        """Generate a binary ground-truth tomogram for this network.

        Project the skeleton voxel positions into a boolean volume
        and optionally dilate to produce a thicker ground truth.

        Args:
            thick (int): Dilation iterations for thickness in
                voxels (default 1, meaning 1-voxel skeleton).

        Returns:
            numpy.ndarray: 3-D boolean ground-truth array.
        """
        hold_gtruth = self.gen_vtp_points_tomo()
        if thick >= 1:
            hold_gtruth = sp_binary_dilation(
                hold_gtruth, iterations=int(thick)
            )
        return hold_gtruth

    def gen_vtp_points_tomo(self):
        """Project skeleton vertices into a binary VOI-shaped tomogram.

        Each skeleton point coordinate is rounded to the nearest
        voxel and set to True; out-of-bounds points are silently
        ignored.

        Returns:
            numpy.ndarray: Boolean 3-D array with the same shape as
                the VOI.
        """
        nx, ny, nz = self._voi.shape
        hold_tomo = np.zeros(shape=(nx, ny, nz), dtype=bool)
        hold_vtp_skel = self.get_skel()
        for i in range(hold_vtp_skel.GetNumberOfPoints()):
            x, y, z = hold_vtp_skel.GetPoint(i)
            x, y, z = int(round(x)), int(round(y)), int(round(z))
            if (
                (x >= 0)
                and (y >= 0)
                and (z >= 0)
                and (x < nx)
                and (y < ny)
                and (z < nz)
            ):
                hold_tomo[x, y, z] = True
        return hold_tomo

    def insert_density_svol(
        self, m_svol, tomo, v_size=1, merge="max", off_svol=None
    ):
        """Stamp all polymers in the network into a target tomogram.

        Args:
            m_svol (numpy.ndarray | list[numpy.ndarray]): Reference
                monomer sub-volume(s).
            tomo (numpy.ndarray): Target tomogram modified in place.
            v_size (float): Voxel size in Angstroms (default 1).
            merge (str): Blending strategy: 'max' (default), 'min',
                'sum', or 'insert'.
            off_svol (array-like, optional): Offset for sub-volume
                centre coordinates.

        Raises:
            TypeError: If m_svol, tomo, or off_svol are of wrong
                type.
            ValueError: If merge is invalid or v_size <= 0.
        """
        if not hasattr(m_svol, "__len__"):
            if not isinstance(m_svol, np.ndarray) or m_svol.ndim != 3:
                raise TypeError("m_svol must be a 3-D numpy array.")
        if not isinstance(tomo, np.ndarray) or tomo.ndim != 3:
            raise TypeError("tomo must be a 3-D numpy array.")
        if merge not in ("max", "min", "sum", "insert"):
            raise ValueError(
                "merge must be 'max', 'min', 'sum' " "or 'insert'."
            )
        if v_size <= 0:
            raise ValueError("v_size must be > 0.")
        if off_svol is not None:
            if not isinstance(off_svol, np.ndarray) or len(off_svol) != 3:
                raise TypeError(
                    "off_svol must be a numpy array " "with 3 elements."
                )

        for pl in self._pl:
            pl.insert_density_svol(
                m_svol, tomo, v_size, merge=merge, off_svol=off_svol
            )

    def add_monomer_to_voi(self, mmer, mmer_svol=None):
        """Carve a monomer's footprint out of the VOI mask.

        Replays the monomer's transformation queue on the binary
        sub-volume and stamps it into the VOI with merge='min' to
        mark occupied voxels as False.

        Args:
            mmer (Monomer): Monomer whose transformation queue is
                replayed.
            mmer_svol (numpy.ndarray): Binary boolean sub-volume
                representing the monomer's spatial footprint.

        Raises:
            TypeError: If mmer_svol is not a boolean numpy array.
        """
        if not isinstance(mmer_svol, np.ndarray) or mmer_svol.dtype != bool:
            raise TypeError("mmer_svol must be a boolean numpy " "array.")
        v_size_i = 1.0 / self._v_size
        tot_v = np.asarray((0.0, 0.0, 0.0))
        hold_svol = mmer_svol > 0
        for trans in mmer.trans_list:
            if trans[0] == "t":
                tot_v += trans[1] * v_size_i
            elif trans[0] == "r":
                hold_svol = tomo_rotate(
                    hold_svol,
                    trans[1],
                    order=0,
                    mode="constant",
                    cval=hold_svol.max(),
                    prefilter=False,
                )
        insert_svol_tomo(hold_svol, self._voi, tot_v, merge="min")

    def count_proteins(self):
        """Collect per-ID monomer count statistics for this network.

        Returns:
            dict[int, int]: Mapping from monomer ID to the number
                of monomers with that ID placed in the network.
        """
        counts = dict()
        for pl in self._pl:
            ids = pl.mmer_ids
            for mmer_id in ids:
                try:
                    counts[mmer_id] += 1
                except KeyError:
                    counts[mmer_id] = 1
        return counts
