"""Networks of helical flexible fibers for cryo-ET volumes.

Defines :class:`NetHelixFiber` and :class:`NetHelixFiberB`, which
place random helical polymer chains inside a 3-D volume of
interest (VOI) up to a user-specified occupancy fraction.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import math
import random
from dataclasses import (
    dataclass,
    field,
)

import numpy as np
import vtk

from .flms_fiber import HelixFiber
from .flms_gen import (
    FlmsParamGen,
    HxParamGenBranched,
)
from ..polymers import Network
from ...utils.affine import poly_translate
from ...utils.utils import (
    points_distance,
    point_to_poly,
)

NET_TYPE_STR = "net_type"


class NetHelixFiber(Network):
    """Network of isolated, unconnected helical fiber polymers.

    Randomly distributes :class:`HelixFiber` chains throughout a
    3-D VOI until a target occupancy percentage is reached.
    """

    def __init__(
        self,
        voi,
        v_size,
        l_length,
        m_surf,
        gen_hfib_params,
        occ,
        min_p_len,
        hp_len,
        mz_len,
        mz_len_f,
        over_tolerance=0,
        unit_diam=None,
    ):
        """Initialise the helical fiber network.

        Args:
            voi (numpy.ndarray): 3-D boolean/numeric VOI array.
            v_size (float): Voxel size in Angstroms.
            l_length (float): Monomer link length; must be > 0.
            m_surf (vtk.vtkPolyData): Monomer reference surface.
            gen_hfib_params (FlmsParamGen): Stochastic parameter
                generator for fiber geometry.
            occ (float): Target occupancy in [0, 100] percent.
            min_p_len (float): Minimum persistence length;
                must be > 0.
            hp_len (float): Helix period length; must be > 0.
            mz_len (float): Monomer axial rise; must be > 0.
            mz_len_f (float): Z-axis elevation factor;
                must be >= 0.
            over_tolerance (float): Allowed surface overlap
                fraction in [0, 1) (default 0).
            unit_diam (float, optional): Structural unit diameter
                for collision tests; None uses the monomer
                diameter.

        Raises:
            ValueError: If any length or occupancy parameter is
                out of range.
            TypeError: If m_surf or gen_hfib_params is the wrong
                type.
        """
        super(NetHelixFiber, self).__init__(voi, v_size)
        if l_length <= 0:
            raise ValueError("l_length must be > 0.")
        if not isinstance(m_surf, vtk.vtkPolyData):
            raise TypeError("m_surf must be a vtkPolyData.")
        if not isinstance(gen_hfib_params, FlmsParamGen):
            raise TypeError(
                "gen_hfib_params must be a " "FlmsParamGen instance."
            )
        if occ < 0 or occ > 100:
            raise ValueError("occ must be in [0, 100].")
        if over_tolerance < 0 or over_tolerance > 100:
            raise ValueError("over_tolerance must be in [0, 100].")
        if min_p_len <= 0 or hp_len <= 0 or mz_len <= 0 or mz_len_f < 0:
            raise ValueError(
                "min_p_len, hp_len, mz_len must be "
                "> 0 and mz_len_f must be >= 0."
            )

        self.__l_length, self.__m_surf = l_length, m_surf
        self.__gen_hfib_params = gen_hfib_params
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__min_p_len, self.__hp_len = min_p_len, hp_len
        self.__mz_len, self.__mz_len_f = mz_len, mz_len_f
        self.__unit_diam = unit_diam

    def build_network(self):
        """Grow helix fiber chains until the occupancy target is met.

        Each fiber starts at a random point inside the VOI and
        extends until it reaches the diagonal length of the VOI or
        monomer placement fails.  Short fibers below _min_nmmer
        monomers are discarded.
        """

        MAX_TRIES = 1000
        tries_count = 0

        while self._pl_occ < self.__occ and tries_count < MAX_TRIES:
            tries_count += 1

            p0 = np.asarray(
                (
                    self._voi.shape[0] * self._v_size * random.random(),
                    self._voi.shape[1] * self._v_size * random.random(),
                    self._voi.shape[2] * self._v_size * random.random(),
                )
            )
            max_length = (
                math.sqrt(
                    self._voi.shape[0] ** 2
                    + self._voi.shape[1] ** 2
                    + self._voi.shape[2] ** 2
                )
                * self._v_size
            )
            p_len = self.__gen_hfib_params.gen_persistence_length(
                self.__min_p_len
            )
            z_len_f = self.__gen_hfib_params.gen_zf_length(self.__mz_len_f)
            hold_polymer = HelixFiber(
                self.__l_length,
                self.__m_surf,
                p_len,
                self.__hp_len,
                self.__mz_len,
                z_len_f,
                p0,
            )

            not_finished = True
            while (hold_polymer.total_len < max_length) and not_finished:
                monomer_data = hold_polymer.gen_new_monomer(
                    self.__over_tolerance,
                    self._voi,
                    self._v_size,
                    net=self,
                    max_dist=self.__unit_diam,
                )
                if monomer_data is None:
                    not_finished = False
                else:
                    new_len = points_distance(
                        monomer_data[0], hold_polymer.tail_point
                    )
                    if hold_polymer.total_len + new_len < max_length:
                        hold_polymer.add_monomer(
                            monomer_data[0],
                            monomer_data[1],
                            monomer_data[2],
                            monomer_data[3],
                        )
                    else:
                        not_finished = False

            # Updating polymer
            if hold_polymer.num_mmers >= self._min_nmmer:
                self.add_polymer(hold_polymer)


class NetHelixFiberB(Network):
    """Network of branched helical fiber polymers.

    Like :class:`NetHelixFiber`, but new fibers can optionally
    branch from existing ones based on a Bernoulli branching
    probability.
    """

    def __init__(
        self,
        voi,
        v_size,
        l_length,
        m_surf,
        gen_hfib_params,
        occ,
        min_p_len,
        hp_len,
        mz_len,
        mz_len_f,
        b_prop,
        max_p_branch=0,
        over_tolerance=0,
        unit_diam=None,
    ):
        """Initialise the branched helical fiber network.

        Args:
            voi (numpy.ndarray): 3-D boolean/numeric VOI array.
            v_size (float): Voxel size in Angstroms.
            l_length (float): Monomer link length; must be > 0.
            m_surf (vtk.vtkPolyData): Monomer reference surface.
            gen_hfib_params (HxParamGenBranched): Branching-aware
                stochastic parameter generator.
            occ (float): Target occupancy in [0, 100] percent.
            min_p_len (float): Minimum persistence length;
                must be > 0.
            hp_len (float): Helix period length; must be > 0.
            mz_len (float): Monomer axial rise; must be > 0.
            mz_len_f (float): Z-axis elevation factor;
                must be >= 0.
            b_prop (float): Branching probability in [0, 1);
                checked each time a new polymer is placed.
            max_p_branch (int): Maximum branches per polymer;
                0 disables branching (default 0).
            over_tolerance (float): Allowed surface overlap
                fraction in [0, 1) (default 0).
            unit_diam (float, optional): Structural unit diameter
                for collision tests; None uses the monomer
                diameter.

        Raises:
            ValueError: If any parameter is out of range.
            TypeError: If m_surf or gen_hfib_params is the wrong
                type.
        """

        super(NetHelixFiberB, self).__init__(voi, v_size)

        if l_length <= 0:
            raise ValueError("l_length must be > 0.")
        if not isinstance(m_surf, vtk.vtkPolyData):
            raise TypeError("m_surf must be a vtkPolyData.")
        if not isinstance(gen_hfib_params, HxParamGenBranched):
            raise TypeError(
                "gen_hfib_params must be a " "HxParamGenBranched instance."
            )
        if occ < 0 or occ > 100:
            raise ValueError("occ must be in [0, 100].")
        if over_tolerance < 0 or over_tolerance > 100:
            raise ValueError("over_tolerance must be in [0, 100].")
        if min_p_len <= 0 or hp_len <= 0 or mz_len <= 0 or mz_len_f < 0:
            raise ValueError(
                "min_p_len, hp_len, mz_len must be "
                "> 0 and mz_len_f must be >= 0."
            )
        if max_p_branch < 0 or b_prop < 0:
            raise ValueError("max_p_branch and b_prop must be " ">= 0.")

        self.__l_length, self.__m_surf = l_length, m_surf
        self.__gen_hfib_params = gen_hfib_params
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__min_p_len, self.__hp_len = min_p_len, hp_len
        self.__mz_len, self.__mz_len_f = mz_len, mz_len_f
        self.__max_p_branch, self.__p_branches, self.__b_prop = (
            max_p_branch,
            list(),
            b_prop,
        )
        self.__unit_diam = unit_diam

    def build_network(self):
        """Grow branched helix fibers until the occupancy target is met.

        On each iteration, randomly chooses to branch from an
        existing fiber or plant a fresh fiber at a random VOI
        location.  Continues until _pl_occ >= occ or MAX_TRIES
        iterations are exhausted.
        """

        MAX_TRIES = 1000
        tries_count = 0

        while self._pl_occ < self.__occ and tries_count < MAX_TRIES:
            tries_count += 1

            max_length = (
                math.sqrt(
                    self._voi.shape[0] ** 2
                    + self._voi.shape[1] ** 2
                    + self._voi.shape[2] ** 2
                )
                * self._v_size
            )
            p_len = self.__gen_hfib_params.gen_persistence_length(
                self.__min_p_len
            )
            z_len_f = self.__gen_hfib_params.gen_zf_length(self.__mz_len_f)
            branch = None
            if (self.__max_p_branch > 0) and self.__gen_hfib_params.gen_branch(
                self.__b_prop
            ):
                branch = self.__gen_random_branch()
            if branch is None:
                p0 = np.asarray(
                    (
                        self._voi.shape[0] * self._v_size * random.random(),
                        self._voi.shape[1] * self._v_size * random.random(),
                        self._voi.shape[2] * self._v_size * random.random(),
                    )
                )
            else:
                p0 = branch.point
            hold_polymer = HelixFiber(
                self.__l_length,
                self.__m_surf,
                p_len,
                self.__hp_len,
                self.__mz_len,
                z_len_f,
                p0,
            )

            not_finished = True
            while (hold_polymer.total_len < max_length) and not_finished:
                monomer_data = hold_polymer.gen_new_monomer(
                    self.__over_tolerance,
                    self._voi,
                    self._v_size,
                    net=self,
                    max_dist=self.__unit_diam,
                )
                if monomer_data is None:
                    not_finished = False
                else:
                    new_len = points_distance(
                        monomer_data[0], hold_polymer.tail_point
                    )
                    if hold_polymer.total_len + new_len < max_length:
                        hold_polymer.add_monomer(
                            monomer_data[0],
                            monomer_data[1],
                            monomer_data[2],
                            monomer_data[3],
                        )
                    else:
                        not_finished = False

            # Updating polymer
            if hold_polymer.num_mmers >= self._min_nmmer:
                if branch is not None:
                    self.add_polymer(hold_polymer)
                    self.__p_branches.append(list())
                    self.__add_branch(hold_polymer, branch)
                else:
                    self.add_polymer(hold_polymer)
                    self.__p_branches.append(list())

    @property
    def branch_list(self):
        """Flatten all per-polymer branch lists into a single list.

        Returns:
            list[Branch]: All Branch objects in the network.
        """
        hold_list = list()
        for bl in self.__p_branches:
            for b in bl:
                hold_list.append(b)
        return hold_list

    def get_skel(self, add_verts=True, add_lines=True, verts_rad=0):
        """Build a vtkPolyData skeleton including polymer lines and branch points.

        Polymer lines are annotated with NET_TYPE_STR = 1;
        branch vertices are annotated with NET_TYPE_STR = 2.

        Args:
            add_verts (bool): Include vertex glyphs (default True).
            add_lines (bool): Include connecting lines (default True).
            verts_rad (float): Sphere radius for vertex glyphs (default 0).

        Returns:
            vtk.vtkPolyData: Combined skeleton dataset.
        """
        if len(self._pl) == 0:
            return vtk.vtkPolyData()

        app_flt_l, app_flt_v, app_flt = (
            vtk.vtkAppendPolyData(),
            vtk.vtkAppendPolyData(),
            vtk.vtkAppendPolyData(),
        )

        p_type_l = vtk.vtkIntArray()
        p_type_l.SetName(NET_TYPE_STR)
        p_type_l.SetNumberOfComponents(1)
        for pol in self._pl:
            app_flt_l.AddInputData(
                pol.get_skel(add_verts, add_lines, verts_rad)
            )
        app_flt_l.Update()
        out_vtp_l = app_flt_l.GetOutput()
        for _ in range(out_vtp_l.GetNumberOfCells()):
            p_type_l.InsertNextTuple((1,))
        out_vtp_l.GetCellData().AddArray(p_type_l)

        p_type_v = vtk.vtkIntArray()
        p_type_v.SetName(NET_TYPE_STR)
        p_type_v.SetNumberOfComponents(1)
        for _, branch in enumerate(self.branch_list):
            app_flt_v.AddInputData(branch.get_vtp())
        app_flt_v.Update()
        out_vtp_v = app_flt_v.GetOutput()
        for _ in range(out_vtp_v.GetNumberOfCells()):
            p_type_v.InsertNextTuple((2,))
        out_vtp_v.GetCellData().AddArray(p_type_v)

        app_flt.AddInputData(out_vtp_l)
        app_flt.AddInputData(out_vtp_v)
        app_flt.Update()

        return app_flt.GetOutput()

    def get_branches_vtp(self, shape_vtp=None):
        """Build a vtkPolyData with one glyph per branch.

        Args:
            shape_vtp (vtk.vtkPolyData, optional): Glyph shape to
                place at each branch; None (default) uses a single
                point vertex.

        Returns:
            vtk.vtkPolyData: Dataset with branch glyphs.
        """

        app_flt = vtk.vtkAppendPolyData()

        for branch in self.branch_list:
            app_flt.AddInputData(branch.get_vtp(shape_vtp))
        app_flt.Update()

        return app_flt.GetOutput()

    def __gen_random_branch(self):
        """Draw a random branch site from an existing polymer.

        Selects a polymer weighted by monomer count, picks a
        random monomer, and returns a new :class:`Branch` only if
        the site does not overlap an existing branch.  Returns
        None if no valid site is found.

        Returns:
            Branch | None: New branch object, or None.
        """

        count, branch = 0, None
        while (count < len(self._pl)) and (branch is None):
            hold_pid = random.choices(
                range(0, len(self._pl)),
                weights=self._pl_nmmers,
            )[0]
            if len(self.__p_branches[hold_pid]) < self.__max_p_branch:
                hold_pol = self._pl[hold_pid]
                hold_mid = random.randint(0, hold_pol.num_monomers - 1)
                hold_m = hold_pol.get_monomer(hold_mid)
                found = True
                for branch in self.__p_branches[hold_pid]:
                    if (
                        points_distance(hold_m.center_mass, branch.point)
                        <= 2 * hold_m.diameter
                    ):
                        found = False
                if found:
                    branch = Branch(hold_m.center_mass, hold_pid, hold_mid)
            count += 1

        return branch

    def __add_branch(self, polymer, branch):
        """Register a branch into the per-polymer branch registry.

        Args:
            polymer (Polymer): The newly placed polymer that
                originates from the branch.
            branch (Branch): Branch object recording the source
                polymer and monomer indices.
        """
        branch.set_t_pmer(len(self._pl) - 1)
        self.__p_branches[branch.s_pmer_id].append(branch)


@dataclass
class Branch:
    """Data record for a branching site in a fiber network."""

    point: np.ndarray = field(default_factory=lambda: np.zeros(3))
    s_pmer_id: int = 0
    s_mmer_id: int = 0
    t_pmer_id: int = None

    def __post_init__(self):
        if not hasattr(self.point, "__len__") or len(self.point) != 3:
            raise ValueError("point must have exactly 3 elements.")
        self.point = np.asarray(self.point, dtype=float)
        if self.s_pmer_id < 0 or self.s_mmer_id < 0:
            raise ValueError("s_pmer_id and s_mmer_id must be " ">= 0.")
        if self.t_pmer_id is not None and self.t_pmer_id < 0:
            raise ValueError("t_pmer_id must be >= 0.")

    def set_t_pmer(self, t_pmer_id):
        """Set the target (child) polymer ID.

        Args:
            t_pmer_id (int): Non-negative polymer index.

        Raises:
            ValueError: If t_pmer_id < 0.
        """
        if t_pmer_id < 0:
            raise ValueError("t_pmer_id must be >= 0.")
        self.t_pmer_id = t_pmer_id

    def get_vtp(self, shape_vtp=None):
        """Return a vtkPolyData glyph for this branch.

        Args:
            shape_vtp (vtk.vtkPolyData, optional): Glyph shape
                translated to the branch point.  None (default)
                returns a single-vertex poly at the branch point.

        Returns:
            vtk.vtkPolyData: Glyph polygon dataset.
        """
        if shape_vtp is None:
            return point_to_poly(self.point)
        else:
            return poly_translate(shape_vtp, self.point)
