"""
Classes for modeling networks polymers.
A networks is a combination of a polymer in a volume
"""

__author__ = 'Antonio Martinez-Sanchez'

from .utils import *
from .affine import *
from .polymer import SAWLC
from .lrandom import PGen
from abc import ABC, abstractmethod


class Network(ABC):
    """
    General class for a network of polymers
    """

    def __init__(self, voi, v_size):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        """
        assert isinstance(voi, np.ndarray)
        self.__voi = voi
        self.__vol = (self.__voi > 0).sum() * v_size * v_size * v_size
        self.__v_size = v_size
        self.__pl_occ = 0
        self.__pl = list()

    def get_polymer_occupancy(self):
        return self.__pl_occ

    def add_polymer(self, polymer):
        self.__pl.append(polymer)
        self.__pl_occ += 100. * (polymer.get_vol() / self._Network__vol)

    @abstractmethod
    def build_network(self):
        """
        Builds an instance of the network
        :return: None        """
        raise NotImplemented

    def get_voi(self):
        """
        Get the VOI
        :return: an ndarray
        """
        return self.__voi

    def get_vtp(self):
        """
        Get Polymers Network as a vtkPolyData with their surfaces
        :return: a vtkPolyData
        """

        app_flt = vtk.vtkAppendPolyData()

        # Polymers loop
        for pol in self.__pl:
            app_flt.AddInputData(pol.get_vtp())
        app_flt.Update()

        return app_flt.GetOutput()

    def get_skel(self):
        """
        Get Polymers Network as a vtkPolyData as points and lines
        :return: a vtkPolyData
        """

        app_flt = vtk.vtkAppendPolyData()

        # Polymers loop
        for pol in self.__pl:
            app_flt.AddInputData(pol.get_skel())
        app_flt.Update()

        return app_flt.GetOutput()

    def insert_density_svol(self, m_svol, tomo, v_size=1, merge='max'):
        """
        Insert a polymer network as set of subvolumes into a tomogram
        :param m_svol: input monomer sub-volume reference
        :param tomo: tomogram where m_svol is added
        :param v_size: tomogram voxel size (default 1)
        :param merge: merging mode, valid: 'max' (default), 'min', 'sum' and 'insert'
        :return:
        """
        for pl in self.__pl:
            pl.insert_density_svol(m_svol, tomo, v_size, merge=merge)


class NetSAWLC(Network):
    """
    Class for generating a network of SAWLC polymers
    """

    def __init__(self, voi, v_size, l_length, m_surf, gen_pol_lengths, occ, over_tolerance=0):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_length: polymer link length
        :param m_surf: monomer surf
        :param gen_pol_lengths: a instance of a random generation model (random.PGen) to determine the achievable
        lengths for polymers
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        """

        # Initialize abstract varibles
        super(NetSAWLC, self).__init__(voi, v_size)

        # Input parsing
        assert l_length > 0
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert issubclass(gen_pol_lengths.__class__, PGen)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)

        # Variables assignment
        self.__l_length, self.__m_surf = l_length, m_surf
        self.__gen_pol_lengths = gen_pol_lengths
        self.__occ, self.__over_tolerance = occ, over_tolerance

    def build_network(self):
        """
        Add polymers following SAWLC model un an occupancy limit is passed
        :return:
        """

        # Network loop
        while self._Network__pl_occ < self.__occ:

            # Polymer initialization
            p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            max_length = self.__gen_pol_lengths.gen_next_length()
            hold_polymer = SAWLC(self.__l_length, self.__m_surf, p0)

            # Polymer loop
            not_finished = True
            while (hold_polymer.get_total_len() < max_length) and not_finished:
                monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                            self._Network__v_size)
                if monomer_data is None:
                    not_finished = False
                else:
                    new_len = points_distance(monomer_data[0], hold_polymer.get_tail_point())
                    if hold_polymer.get_total_len() + new_len < max_length:
                        hold_polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3])
                    else:
                        not_finished = False

            # Updating polymer
            self.add_polymer(hold_polymer)
            print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
                  ' and length ' + str(hold_polymer.get_total_len()))

