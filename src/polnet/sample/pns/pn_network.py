"""SAWLC polymer network for cytosolic protein placement.

Defines :class:`NetSAWLC`, which distributes self-avoiding
worm-like chain (SAWLC) polymers inside a 3-D VOI up to a
target occupancy fraction.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import random

import numpy as np
import vtk

from .swalc import (
    SAWLC,
    SAWLCPoly,
)
from ..polymers import Network
from ...utils.poly import poly_surface_area
from ...utils.utils import points_distance


class NetSAWLC(Network):
    """Network of Self-Avoiding Worm-Like Chain (SAWLC) polymers.

    Distributes cytosolic SAWLC polymer chains inside a 3-D VOI
    until a target occupancy percentage is reached.
    """

    def __init__(
        self,
        voi,
        v_size,
        l_length,
        m_surf,
        max_p_length,
        occ,
        over_tolerance=0,
        poly=None,
        svol=None,
        tries_mmer=50,
        tries_pmer=10,
        rots=None,
        rot_id=0,
    ):
        """Initialise a SAWLC polymer network.

        Args:
            voi (numpy.ndarray): 3-D boolean/numeric VOI array.
            v_size (float): Voxel size in Angstroms.
            l_length (float): Link (step) length in Angstroms;
                must be > 0.
            m_surf (vtk.vtkPolyData): Monomer reference surface.
            max_p_length (float): Maximum allowed polymer length
                in Angstroms.
            occ (float): Target occupancy percentage in [0, 100].
            over_tolerance (float): Allowed surface overlap
                fraction in [0, 1) (default 0).
            poly (vtk.vtkPolyData, optional): If provided,
                restricts monomer placement to this surface.
            svol (numpy.ndarray, optional): Monomer sub-volume
                for VOI carving (default None).
            tries_mmer (int): Placement attempts per monomer
                before aborting the polymer (default 50).
            tries_pmer (int): Total polymer placement retries
                (default 10).
            rots (array-like, optional): Pre-defined rotation
                quaternion array; consumed sequentially if
                provided.
            rot_id (int): Starting index into rots (default 0).

        Raises:
            ValueError: If l_length, max_p_length, occ,
                over_tolerance, tries_mmer, or tries_pmer are out
                of range.
            TypeError: If m_surf or poly is not a vtkPolyData.
        """

        super(NetSAWLC, self).__init__(voi, v_size, svol=svol)

        if l_length <= 0:
            raise ValueError("l_length must be > 0.")
        if not isinstance(m_surf, vtk.vtkPolyData):
            raise TypeError("m_surf must be a vtkPolyData.")
        if max_p_length < 0:
            raise ValueError("max_p_length must be >= 0.")
        if occ < 0 or occ > 100:
            raise ValueError("occ must be in [0, 100].")
        if over_tolerance < 0 or over_tolerance > 100:
            raise ValueError("over_tolerance must be in [0, 100].")
        if poly is not None:
            if not isinstance(poly, vtk.vtkPolyData):
                raise TypeError("poly must be a vtkPolyData.")
        if tries_mmer < 1:
            raise ValueError("tries_mmer must be >= 1.")
        if tries_pmer < 1:
            raise ValueError("tries_pmer must be >= 1.")
        self.__rots, self.__rot_id = None, 0
        if rots is not None:
            if len(rots) == 0 or len(rots[0]) != 4:
                raise ValueError(
                    "rots must be a non-empty array "
                    "of 4-element quaternions."
                )
            self.__rots = rots
            self.__rot_id = int(rot_id)

        self.__l_length, self.__m_surf = l_length, m_surf
        self.__max_p_length = max_p_length
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__poly = poly
        self.__tries_mmer = int(tries_mmer)
        self.__tries_pmer = int(tries_pmer)
        if self.__poly is not None:
            self._poly_area = poly_surface_area(self.__poly)

    def build_network(self):
        """Grow SAWLC chains until the occupancy target is met.

        Places polymers one at a time; each polymer grows until
        max_p_length is reached or monomer placement fails after
        tries_mmer attempts.  The whole procedure retries up to
        tries_pmer times.
        """

        c_try = 0
        self._pmer_fails = 0
        if self.__rots is not None:
            rot_id = self.__rot_id

        while (c_try <= self.__tries_pmer) and (self._pl_occ < self.__occ):

            c_try += 1
            if self.__poly:
                p0 = np.asarray(
                    self.__poly.GetPoint(
                        random.randint(0, self.__poly.GetNumberOfPoints() - 1)
                    )
                )
            else:
                p0 = np.asarray(
                    (
                        self._voi.shape[0] * self._v_size * random.random(),
                        self._voi.shape[1] * self._v_size * random.random(),
                        self._voi.shape[2] * self._v_size * random.random(),
                    )
                )
            max_length = random.uniform(0, self.__max_p_length)
            hold_rot = None
            if self.__rots is not None:
                hold_rot = self.__rots[rot_id]
            if self.__poly is None:
                hold_polymer = SAWLC(
                    self.__l_length, self.__m_surf, p0, rot=hold_rot
                )
            else:
                hold_polymer = SAWLCPoly(
                    self.__poly,
                    self.__l_length,
                    self.__m_surf,
                    p0,
                    rot=hold_rot,
                )
            if hold_polymer.get_monomer(-1).overlap_voi(
                self._voi,
                self._v_size,
                over_tolerance=self.__over_tolerance,
            ):
                self._pmer_fails += 1
                continue
            self.add_monomer_to_voi(hold_polymer.get_monomer(-1), self._svol)
            if self.__rots is not None:
                if rot_id >= len(self.__rots) - 1:
                    rot_id = 0
                else:
                    rot_id += 1

            cont_pol = 1
            not_finished = True
            while (hold_polymer.total_len < max_length) and not_finished:
                hold_rot = None
                if self.__rots is not None:
                    hold_rot = self.__rots[rot_id]
                monomer_data = hold_polymer.gen_new_monomer(
                    self.__over_tolerance,
                    self._voi,
                    self._v_size,
                    fix_dst=random.uniform(
                        self.__l_length, 2 * self.__l_length
                    ),
                    rot=hold_rot,
                )

                cont_pol += 1

                if monomer_data is None:
                    if cont_pol >= self.__tries_mmer:
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
                        self.add_monomer_to_voi(
                            hold_polymer.get_monomer(-1), self._svol
                        )
                        hold_occ = self._pl_occ + 100.0 * (
                            hold_polymer.vol / self._vol
                        )
                        if self.__rots is not None:
                            if rot_id >= len(self.__rots) - 1:
                                rot_id = 0
                            else:
                                rot_id += 1
                        if hold_occ >= self.__occ:
                            not_finished = False
                    else:
                        not_finished = False

            # Updating polymer
            if self.__poly is None:
                self.add_polymer(hold_polymer, occ_mode="volume")
                c_try = 0
            else:
                self.add_polymer(hold_polymer, occ_mode="area")
                c_try = 0
