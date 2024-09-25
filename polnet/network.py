"""
Classes for modeling networks polymers.
A networks is a combination of a polymer in a volume
"""

__author__ = 'Antonio Martinez-Sanchez'

import scipy as sp

from polnet.poly import *
from polnet.polymer import SAWLC, SAWLCPoly, HelixFiber
from polnet.lrandom import SGen, SGenUniform, SGenFixed, SGenProp, PGenHelixFiber, PGenHelixFiberB
from abc import ABC, abstractmethod

NET_TYPE_STR = 'net_type'


class Network(ABC):
    """
    General class for a network of polymers
    """

    def __init__(self, voi, v_size, svol=None):
        """
        Construction

        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param svol: monomer subvolume (or list of) as a numpy ndarray (default None)
        :param mb_area: total membrane area within the same VOI as the network (deftault None)
        """
        self.set_voi(voi)
        self.__vol = float((self.__voi > 0).sum()) * v_size * v_size * v_size # withou the float cast is my raise overflow warning in Windows
        self.__v_size = v_size
        self.__pl_occ = 0
        self.__pl = list()
        self.__pl_nmmers = list()
        self.__svol = svol
        if self.__svol is not None:
            if not hasattr(svol, '__len__'):
                assert isinstance(self.__svol, np.ndarray)
        self.__min_nmmer = 1
        self.__poly_area = 0
        self.__pmer_fails = 0

    def set_min_nmmer(self, min_nmmer):
        """
        Set a minimum number of monomers for the generated filaments

        :param min_nmmer: integer with the minimum number of monomenrs per filament
        :return:
        """
        self.__min_nmmer = int(min_nmmer)

    def get_pmer_fails(self):
        return self.__pmer_fails

    def get_pmers_list(self):
        return self.__pl

    def get_num_pmers(self):
        """
        :return: the number of polymers in the network
        """
        return len(self.__pl)

    def get_num_mmers(self):
        """
        :return: the number of monomers in the network
        """
        count_mmers = 0
        for pl in self.__pl:
            count_mmers += pl.get_num_mmers()
        return count_mmers

    def get_polymer_occupancy(self):
        return self.__pl_occ

    def add_polymer(self, polymer, occ_mode='volume'):
        """
        Add a new polymer to the network

        :param polymer: polymer to add
        :param occ_mode: occupancy mode, valid: 'volume' (default), 'area' for membrane-bound polymer
        :return:
        """
        assert (occ_mode == 'volume') or (occ_mode == 'area')
        self.__pl.append(polymer)
        self.__pl_nmmers.append(polymer.get_num_mmers())
        if occ_mode == 'volume':
            self.__pl_occ += 100. * (polymer.get_vol() / self.__vol)
        else:
            self.__pl_occ += 100. * (polymer.get_area() / self._Network__poly_area)
        # print('Occ: ', self.__pl_occ)

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

    def get_gtruth(self, thick=1):
        """
        Get the ground truth tomogram

        :param thick: ground truth tickness in voxels (default 1)
        :return: a binary numpy 3D array
        """
        hold_gtruth = self.gen_vtp_points_tomo()
        if thick >= 1:
            hold_gtruth = sp.ndimage.morphology.binary_dilation(hold_gtruth, iterations=int(thick))
        return hold_gtruth

    def set_voi(self, voi):
        """
        Set the VOI

        :param voi:
        """
        assert isinstance(voi, np.ndarray)
        if voi.dtype is bool:
            self.__voi = voi
        else:
            self.__voi = voi > 0

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

    def get_skel(self, add_verts=True, add_lines=True, verts_rad=0):
        """
        Get the polymer as a skeleton, each monomer is a point or sphere and lines connecting monomers

        :param add_verts: if True (default) the vertices are included in the vtkPolyData
        :param add_lines: if True (default) the lines are included in the vtkPolyData
        :param verts_rad: if verts is True then sets the vertex radius, if <=0 a vertices are just points
        :return: a vtkPolyData
        """
        app_flt = vtk.vtkAppendPolyData()

        # Polymers loop
        for pol in self.__pl:
            app_flt.AddInputData(pol.get_skel(add_verts, add_lines, verts_rad))

        # Update and return
        app_flt.Update()
        return app_flt.GetOutput()

    def gen_vtp_points_tomo(self):
        """
        Generates a binary tomogram where True elements correspond with the polydata closes voxel projection

        :return: a binary VOI shaped numpy array
        """
        nx, ny, nz = self.__voi.shape
        hold_tomo = np.zeros(shape=(nx, ny, nz), dtype=bool)
        hold_vtp_skel = self.get_skel()
        for i in range(hold_vtp_skel.GetNumberOfPoints()):
            x, y, z = hold_vtp_skel.GetPoint(i)
            x, y, z = int(round(x)), int(round(y)), int(round(z))
            if (x >= 0) and (y >= 0) and (z >= 0) and  (x < nx) and (y < ny) and (z < nz):
                hold_tomo[x, y, z] = True
        return hold_tomo

    def insert_density_svol(self, m_svol, tomo, v_size=1, merge='max', off_svol=None):
        """
        Insert a polymer network as set of subvolumes into a tomogram

        :param m_svol: input monomer (or list) sub-volume reference
        :param tomo: tomogram where m_svol is added
        :param v_size: tomogram voxel size (default 1)
        :param merge: merging mode, valid: 'min' (default), 'max', 'sum' and 'insert'
        :param off_svol: offset coordinates for sub-volume monomer center coordinates
        """
        if not hasattr(m_svol, '__len__'):
            assert isinstance(m_svol, np.ndarray) and (len(m_svol.shape) == 3)
        assert isinstance(tomo, np.ndarray) and (len(tomo.shape) == 3)
        assert (merge == 'max') or (merge == 'min') or (merge == 'sum') or (merge == 'insert')
        assert v_size > 0
        if off_svol is not None:
            assert isinstance(off_svol, np.ndarray) and (len(off_svol) == 3)

        for pl in self.__pl:
            pl.insert_density_svol(m_svol, tomo, v_size, merge=merge, off_svol=off_svol)

    def add_monomer_to_voi(self, mmer, mmer_svol=None):
        """
        Adds a monomer to VOI mask

        :param mmer: monomer to define rigid transformations
        :param mmer_voi: subvolume (binary numpy ndarray) with monomer VOI
        """
        assert isinstance(mmer_svol, np.ndarray) and (mmer_svol.dtype == bool)
        v_size_i = 1. / self.__v_size
        tot_v = np.asarray((0., 0., 0.))
        hold_svol = mmer_svol > 0
        for trans in mmer.get_trans_list():
            if trans[0] == 't':
                tot_v += (trans[1] * v_size_i)
            elif trans[0] == 'r':
                hold_svol = tomo_rotate(hold_svol, trans[1], order=0, mode='constant', cval=hold_svol.max(),
                                        prefilter=False)
                # hold_svol = tomo_rotate(hold_svol, trans[1], mode='constant', cval=hold_svol.min())
        insert_svol_tomo(hold_svol, self.__voi, tot_v, merge='min')

    def count_proteins(self):
        """
        Genrrates output statistics for this network

        :return: a dictionary with the number of proteins for protein id
        """
        counts = dict()
        for pl in self.__pl:
            ids = pl.get_mmer_ids()
            for id in ids:
                try:
                    counts[id] += 1
                except KeyError:
                    counts[id] = 0
        return counts


class NetSAWLC(Network):
    """
    Class for generating a network of SAWLC polymers
    """

    def __init__(self, voi, v_size, l_length, m_surf, max_p_length, gen_pol_lengths, occ, over_tolerance=0,
                 poly=None, svol=None, tries_mmer=20, tries_pmer=100, rots=None, rot_id=0):
        """
        Construction

        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_length: polymer link length
        :param m_surf: monomer surf
        :param max_p_length: maximum polymer length
        :param gen_pol_lengths: a instance of a random generation model (random.PGen) to determine the achievable
        lengths for polymers
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        :param poly: it allows to restrict monomer localizations to a polydata (e.g. a membrane)
        :param svol: monomer subvolume as a numpy ndarray (default None)
        :param off_svol: offset coordinates in voxels for shifting sub-volume monomer center coordinates (default None)
        :param tries_mmer: number of tries to place a monomer before starting a new polymer
        :param tries_pmer: number of tries to place a polymer
        :param rots: allow to externally control the rotations of the macromolecules, if not None (default), the
                     rotations (quaternions) are taken sequentially from the input array of rotations.
        :param rot_id: starting index for rots rotations array.
        """

        # Initialize abstract varibles
        super(NetSAWLC, self).__init__(voi, v_size, svol=svol)

        # Input parsing
        assert l_length > 0
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert max_p_length >= 0
        assert isinstance(gen_pol_lengths, PGenHelixFiber)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        if poly is not None:
            assert isinstance(poly, vtk.vtkPolyData)
        assert tries_mmer >= 1
        assert tries_pmer >= 1
        self.__rots, self.__rot_id = None, 0
        if rots is not None:
            assert len(rots) > 0 and len(rots[0]) == 4
            self.__rots = rots
            self.__rot_id = int(rot_id)

        # Variables assignment
        self.__l_length, self.__m_surf = l_length, m_surf
        self.__max_p_length = max_p_length
        self.__gen_pol_lengths = gen_pol_lengths
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__poly = poly
        self.__tries_mmer = int(tries_mmer)
        self.__tries_pmer = int(tries_pmer)
        self.__poly_area = None
        if self.__poly is not None:
            self._Network__poly_area = poly_surface_area(self.__poly)

    def build_network(self):
        """
        Add polymers following SAWLC model until an occupancy limit is passed

        :return:
        """

        c_try = 0
        self._Network__pmer_fails = 0
        if self.__rots is not None:
            rot_id = self.__rot_id

        # Network loop
        while (c_try <= self.__tries_pmer) and (self._Network__pl_occ < self.__occ):

            # Polymer initialization
            c_try += 1
            if self.__poly:
                p0 = np.asarray(self.__poly.GetPoint(random.randint(0, self.__poly.GetNumberOfPoints())))
            else:
                p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                                 self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                                 self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            max_length = self.__gen_pol_lengths.gen_length(0, self.__max_p_length)
            hold_rot = None
            if self.__rots is not None:
                hold_rot = self.__rots[rot_id]
            if self.__poly is None:
                hold_polymer = SAWLC(self.__l_length, self.__m_surf, p0, rot = hold_rot)
            else:
                hold_polymer = SAWLCPoly(self.__poly, self.__l_length, self.__m_surf, p0, rot = hold_rot)
            if hold_polymer.get_monomer(-1).overlap_voi(self._Network__voi, self._Network__v_size,
                                                        over_tolerance=self.__over_tolerance):
                self._Network__pmer_fails += 1
                continue
            self.add_monomer_to_voi(hold_polymer.get_monomer(-1), self._Network__svol)
            if self.__rots is not None:
                if rot_id >= len(self.__rots) - 1:
                    rot_id = 0
                else:
                    rot_id += 1

            # Polymer loop
            cont_pol = 1
            not_finished = True
            while (hold_polymer.get_total_len() < max_length) and not_finished:
                hold_rot = None
                if self.__rots is not None:
                    hold_rot = self.__rots[rot_id]
                monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                            self._Network__v_size,
                                                            fix_dst=self.__gen_pol_lengths.gen_length(self.__l_length, 2*self.__l_length),
                                                            rot = hold_rot)

                cont_pol += 1

                if monomer_data is None:
                    if cont_pol >= self.__tries_mmer:
                        not_finished = False
                    else:
                        c_try += 1
                else:
                    new_len = points_distance(monomer_data[0], hold_polymer.get_tail_point())
                    if hold_polymer.get_total_len() + new_len < max_length:
                        # ) and (monomer_data[3].overlap_voi(self._Network__voi, self._Network__v_size)):
                        hold_polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3])
                        self.add_monomer_to_voi(hold_polymer.get_monomer(-1), self._Network__svol)
                        hold_occ = self._Network__pl_occ + 100. * (hold_polymer.get_vol() / self._Network__vol)
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
                self.add_polymer(hold_polymer, occ_mode='volume')
                c_try = 0
            else:
                self.add_polymer(hold_polymer, occ_mode='area')
                c_try = 0
            # print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
            #       ', length ' + str(hold_polymer.get_total_len()) + ' and occupancy ' +
            #       str(self._Network__pl_occ) + '%')

        # print('Exit with c_try=' + str(c_try) + ' and c_fails=' + str(self._Network__pmer_fails))


class NetSAWLCInter(Network):
    """
    Class for generating a network of SAWLC polymers with intercalated monomers
    """

    def __init__(self, voi, v_size, l_lengths, m_surfs, max_p_length, gen_pol_lengths, gen_seq_mmers,
                 occ, over_tolerance=0, poly=None, svols=None, codes=None, compaq=False, tries_mmer=100, rots_mmer=10):
        """
        Construction

        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_lengths: polymer link lengths
        :param m_surfs: monomers list with the posible monomer surfaces
        :param max_p_length: maximum polymer length
        :param gen_pol_lengths: an instance of a random generation model (lrandom.PGen) to determine the achievable
        lengths for polymers
        :param gen_seq_mmers: an instance of a random generation model (lrandom.SGen) to generate a sequence of monomers
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        :param poly: it allows to restrict monomer localizations to a polydata
        :param svol: monomer subvolumes as list numpy ndarrays (default None)
        :param codes: monomer codes a a list (default None)
        :param off_svol: offset coordinates in voxels for shifting sub-volume monomer center coordinates (default None)
        :param compaq: if a number (default None) encodes the minimum distance between two surfaces
        :param tries_mmer: number of tries to place a monomer before starting a new polymer
        :param rots_mmer: number of rotations for monomer before increasing distance to provious by the 'compaq' value
        """

        # Initialize abstract varibles
        super(NetSAWLCInter, self).__init__(voi, v_size, svol=svols)

        # Input parsing
        assert hasattr(l_lengths, '__len__')
        assert hasattr(m_surfs, '__len__')
        if codes is not None:
            assert hasattr(codes, '__len__') and (len(codes) == len(svols))
        assert max_p_length >= 0
        assert isinstance(gen_pol_lengths, PGenHelixFiber)
        assert issubclass(type(gen_seq_mmers), SGen)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        if poly is not None:
            assert isinstance(poly, vtk.vtkPolyData)
        assert (rots_mmer > 0) and (tries_mmer > rots_mmer)

        # Variables assignment
        self.__l_lengths, self.__m_surfs = l_lengths, m_surfs
        self.__centers = np.zeros(shape=(len(self.__m_surfs), 3))
        for i in range(len(self.__m_surfs)):
            self.__centers[i, :] = poly_center_mass(self.__m_surfs[i])
        self.__max_p_length = max_p_length
        self.__gen_pol_lengths, self.__gen_seq_mmers = gen_pol_lengths, gen_seq_mmers
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__poly = poly
        self.__codes = codes
        self.__compaq = None
        if compaq is not None:
            self.__compaq = compaq
        self.__tries_mmer, self.__rots_mmer = tries_mmer, rots_mmer

    def build_network(self):
        """
        Add polymers following SAWLC model until an occupancy limit is passed

        :return:
        """

        # Computes the maximum number of tries
        mmer_id = None
        c_try = 0
        n_tries = math.ceil(self._Network__vol) # math.ceil(self._Network__vol / poly_volume(self.__m_surfs[prev_id]))

        # Network loop
        while (c_try <= n_tries) and (self._Network__pl_occ < self.__occ):

            # Polymer initialization
            c_try += 1
            p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            max_length = self.__gen_pol_lengths.gen_length(0, self.__max_p_length)
            if mmer_id is None:
                mmer_id = self.__gen_seq_mmers.gen_next_mmer_id(len(self.__m_surfs))
            if self.__poly is None:
                hold_polymer = SAWLC(self.__l_lengths[mmer_id], self.__m_surfs[mmer_id], p0, id0=mmer_id,
                                     code0=self.__codes[mmer_id])
            else:
                hold_polymer = SAWLCPoly(self.__poly, self.__l_lengths[mmer_id], self.__m_surfs[mmer_id], p0,
                                         id0=mmer_id, code0=self.__codes[mmer_id])
            if hold_polymer.get_monomer(-1).overlap_voi(self._Network__voi, self._Network__v_size,
                                                        self.__over_tolerance):
                continue
            self.add_monomer_to_voi(hold_polymer.get_monomer(-1), self._Network__svol[mmer_id])
            prev_id = mmer_id
            mmer_id = None

            # Polymer loop
            cont_pol, hold_compaq = 0, self.__compaq
            not_finished = True
            while (hold_polymer.get_total_len() < max_length) and not_finished:

                if mmer_id is None:
                    mmer_id = self.__gen_seq_mmers.gen_next_mmer_id(len(self.__m_surfs))

                if self.__compaq is None:
                    monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                                self._Network__v_size, id=mmer_id)
                else:
                    hold_off = poly_point_min_dst(self.__m_surfs[mmer_id], self.__centers[mmer_id]) + \
                               poly_point_min_dst(self.__m_surfs[prev_id], self.__centers[prev_id]) + hold_compaq
                    monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                                self._Network__v_size, hold_off,
                                                                ext_surf=self.__m_surfs[mmer_id])
                cont_pol += 1

                if monomer_data is None:
                    if cont_pol >= self.__tries_mmer:
                        not_finished = False
                    else:
                        c_try += 1
                        if cont_pol % self.__rots_mmer == 0:
                            hold_compaq = cont_pol * self.__compaq
                else:
                    new_len = points_distance(monomer_data[0], hold_polymer.get_tail_point())
                    if hold_polymer.get_total_len() + new_len < max_length:
                        # ) and (monomer_data[3].overlap_voi(self._Network__voi, self._Network__v_size)):
                        hold_polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3],
                                                 id=mmer_id, code=self.__codes[mmer_id])
                        self.add_monomer_to_voi(hold_polymer.get_monomer(-1), self._Network__svol[mmer_id])
                        prev_id = mmer_id
                        mmer_id = None
                        cont_pol, hold_compaq = 0, 0
                        hold_occ = self._Network__pl_occ + 100. * (hold_polymer.get_vol() / self._Network__vol)
                        if hold_occ >= self.__occ:
                            not_finished = False
                    else:
                        not_finished = False

            # Updating polymer
            self.add_polymer(hold_polymer)
            # print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
            #       ', length ' + str(hold_polymer.get_total_len()) + ' and occupancy ' + str(self._Network__pl_occ) + '%')


class NetHelixFiber(Network):
    """
    Class for generating a network of isolated helix fibers, unconnected and randomly distributed
    """

    def __init__(self, voi, v_size, l_length, m_surf, gen_hfib_params, occ, min_p_len, hp_len, mz_len, mz_len_f,
                 over_tolerance=0, unit_diam=None):
        """
        Construction

        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_length: polymer link length
        :param m_surf: monomer surf
        :param gen_hfib_params: a instance of a random generation model (random.PGen) to obtain random fiber
        parametrization
        :param min_p_len: minimum persistence length
        :param hp_len: helix period length
        :param mz_len: monomer length
        :param mz_len_f: maximum length factor in z-axis
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        :param unit_diam: structural unit diameter
        """

        # Initialize abstract variables
        super(NetHelixFiber, self).__init__(voi, v_size)

        # Input parsing
        assert l_length > 0
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert isinstance(gen_hfib_params, PGenHelixFiber)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        assert (min_p_len > 0) and (hp_len > 0) and (mz_len > 0) and (mz_len_f >= 0)

        # Variables assignment
        self.__l_length, self.__m_surf = l_length, m_surf
        self.__gen_hfib_params = gen_hfib_params
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__min_p_len, self.__hp_len = min_p_len, hp_len
        self.__mz_len, self.__mz_len_f = mz_len, mz_len_f
        self.__unit_diam = unit_diam

    def build_network(self):
        """
        Add helix fibres until an occupancy limit is passed

        :return:
        """

        # Network loop
        while self._Network__pl_occ < self.__occ:

            # Polymer initialization
            p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                             self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            max_length = math.sqrt(self._Network__voi.shape[0]**2 + self._Network__voi.shape[1]**2
                                  + self._Network__voi.shape[2]**2) * self._Network__v_size
            p_len = self.__gen_hfib_params.gen_persistence_length(self.__min_p_len)
            z_len_f = self.__gen_hfib_params.gen_zf_length(self.__mz_len_f)
            hold_polymer = HelixFiber(self.__l_length, self.__m_surf, p_len, self.__hp_len, self.__mz_len, z_len_f, p0)

            # Polymer loop
            not_finished = True
            while (hold_polymer.get_total_len() < max_length) and not_finished:
                monomer_data = hold_polymer.gen_new_monomer(self.__over_tolerance, self._Network__voi,
                                                            self._Network__v_size, net=self, max_dist=self.__unit_diam)
                if monomer_data is None:
                    not_finished = False
                else:
                    new_len = points_distance(monomer_data[0], hold_polymer.get_tail_point())
                    if hold_polymer.get_total_len() + new_len < max_length:
                        hold_polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3])
                    else:
                        not_finished = False

            # Updating polymer
            if hold_polymer.get_num_mmers() >= self._Network__min_nmmer:
                self.add_polymer(hold_polymer)
                # print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
                #       ' and length ' + str(hold_polymer.get_total_len()) + ': occupancy ' + str(self._Network__pl_occ))


class NetHelixFiberB(Network):
    """
    Class for generating a network of brancked helix fibers randomly distributed
    """

    def __init__(self, voi, v_size, l_length, m_surf, gen_hfib_params, occ, min_p_len, hp_len, mz_len, mz_len_f,
                 b_prop, max_p_branch=0, over_tolerance=0):
        """
        Construction

        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param l_length: polymer link length
        :param m_surf: monomer surf
        :param gen_hfib_params: a instance of a random generation model (random.PGen.NetHelixFiberB) to obtain random fiber
        parametrization for helix with branches
        :param occ: occupancy threshold in percentage [0, 100]%
        :param min_p_len: minimum persistence length
        :param hp_len: helix period length
        :param mz_len: monomer length in z-axis
        :param mz_len_f: maximum length factor in z-axis
        :param b_prob: branching probability, checked every time a new monomer is added
        :param max_p_branch: maximum number of branches per polymer, if 0 (default) then no branches are generated
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        """

        # Initialize abstract variables
        super(NetHelixFiberB, self).__init__(voi, v_size)

        # Input parsing
        assert l_length > 0
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert isinstance(gen_hfib_params, PGenHelixFiberB)
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        assert (min_p_len > 0) and (hp_len > 0) and (mz_len > 0) and (mz_len_f >= 0)
        assert (max_p_branch >= 0) and (b_prop >= 0)

        # Variables assignment
        self.__l_length, self.__m_surf = l_length, m_surf
        self.__gen_hfib_params = gen_hfib_params
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__min_p_len, self.__hp_len = min_p_len, hp_len
        self.__mz_len, self.__mz_len_f = mz_len, mz_len_f
        self.__max_p_branch, self.__p_branches, self.__b_prop = max_p_branch, list(), b_prop

    def build_network(self):
        """
        Add helix fibres until an occupancy limit is passed

        :return:
        """

        # Network loop
        while self._Network__pl_occ < self.__occ:

            # Polymer initialization
            max_length = math.sqrt(self._Network__voi.shape[0]**2 + self._Network__voi.shape[1]**2
                                  + self._Network__voi.shape[2]**2) * self._Network__v_size
            p_len = self.__gen_hfib_params.gen_persistence_length(self.__min_p_len)
            z_len_f = self.__gen_hfib_params.gen_zf_length(self.__mz_len_f)
            branch = None
            if (self.__max_p_branch > 0) and self.__gen_hfib_params.gen_branch(self.__b_prop):
                branch = self.__gen_random_branch()
            if branch is None:
                p0 = np.asarray((self._Network__voi.shape[0] * self._Network__v_size * random.random(),
                                 self._Network__voi.shape[1] * self._Network__v_size * random.random(),
                                 self._Network__voi.shape[2] * self._Network__v_size * random.random()))
            else:
                p0 = branch.get_point()
            hold_polymer = HelixFiber(self.__l_length, self.__m_surf, p_len, self.__hp_len,
                                      self.__mz_len, z_len_f, p0)

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
            if hold_polymer.get_num_mmers() >= self._Network__min_nmmer:
                if branch is not None:
                    self.add_polymer(hold_polymer)
                    self.__p_branches.append(list())
                    self.__add_branch(hold_polymer, branch)
                else:
                    self.add_polymer(hold_polymer)
                    self.__p_branches.append(list())
                # print('build_network: new polymer added with ' + str(hold_polymer.get_num_monomers()) +
                #       ' and length ' + str(hold_polymer.get_total_len()) + ': occupancy ' + str(self._Network__pl_occ))

    def get_branch_list(self):
        """
        Get all branches in a list

        :return: a single list with the branches
        """
        hold_list = list()
        for bl in self.__p_branches:
            for b in bl:
                hold_list.append(b)
        return hold_list

    def get_skel(self):
        """
        Get Polymers Network as a vtkPolyData as points and lines with branches

        :return: a vtkPolyData
        """

        # Initialization
        app_flt_l, app_flt_v, app_flt = vtk.vtkAppendPolyData(), vtk.vtkAppendPolyData(), vtk.vtkAppendPolyData()

        # Polymers loop
        p_type_l = vtk.vtkIntArray()
        p_type_l.SetName(NET_TYPE_STR)
        p_type_l.SetNumberOfComponents(1)
        for pol in self._Network__pl:
            app_flt_l.AddInputData(pol.get_skel())
        app_flt_l.Update()
        out_vtp_l = app_flt_l.GetOutput()
        for i in range(out_vtp_l.GetNumberOfCells()):
            p_type_l.InsertNextTuple((1,))
        out_vtp_l.GetCellData().AddArray(p_type_l)

        # Branches loop
        p_type_v = vtk.vtkIntArray()
        p_type_v.SetName(NET_TYPE_STR)
        p_type_v.SetNumberOfComponents(1)
        for i, branch in enumerate(self.get_branch_list()):
            app_flt_v.AddInputData(branch.get_vtp())
            # print('Point ' + str(i) + ': ' + str(branch.get_point()))
        app_flt_v.Update()
        out_vtp_v = app_flt_v.GetOutput()
        for i in range(out_vtp_v.GetNumberOfCells()):
            p_type_v.InsertNextTuple((2,))
        out_vtp_v.GetCellData().AddArray(p_type_v)

        # Merging branches and polymers
        app_flt.AddInputData(out_vtp_l)
        app_flt.AddInputData(out_vtp_v)
        app_flt.Update()

        return app_flt.GetOutput()

    def get_branches_vtp(self, shape_vtp=None):
        """
        Get Branches as a vtkPolyData with points

        :param shape_vtp: if None (default) the a point is returned, otherwise this shape is used
                          TODO: so far only isotropic shapes are recommended and starting monomer tangent is not considered yet
        :return: a vtkPolyData
        """

        # Initialization
        app_flt_l, app_flt_v, app_flt = vtk.vtkAppendPolyData(), vtk.vtkAppendPolyData(), vtk.vtkAppendPolyData()

        # Branches loop
        for i, branch in enumerate(self.get_branch_list()):
            app_flt_v.AddInputData(branch.get_vtp(shape_vtp))
            # print('Point ' + str(i) + ': ' + str(branch.get_point()))
        app_flt_v.Update()
        out_vtp_v = app_flt_v.GetOutput()

        # Merging branches and polymers
        app_flt.AddInputData(out_vtp_v)
        app_flt.Update()

        return app_flt.GetOutput()

    def __gen_random_branch(self):
        """
        Generates a position point randomly for a branch on the filament network, no more than one branch per polymer

        :return: a branch
        """

        # Loop for polymers
        count, branch = 0, None
        while (count < len(self._Network__pl)) and (branch is None):
            hold_pid = random.choices(range(0, len(self._Network__pl)), weights=self._Network__pl_nmmers)[0]
            if len(self.__p_branches[hold_pid]) < self.__max_p_branch:
                hold_pol = self._Network__pl[hold_pid]
                hold_mid = random.randint(0, len(hold_pol._Polymer__m) - 1)
                hold_m = hold_pol._Polymer__m[hold_mid]
                found = True
                for branch in self.__p_branches[hold_pid]:
                    if points_distance(hold_m.get_center_mass(), branch.get_point()) <= 2 * hold_m.get_diameter():
                        found = False
                if found:
                    branch = Branch(hold_m.get_center_mass(), hold_pid, hold_mid)
            count += 1

        return branch

    def __add_branch(self, polymer, branch):
        """
        Add a new branch to the polymer network

        :param polymer: targeting polymer where the branch is going to be added (starting polumer is obtained from
                        the branch)
        :param branch: branch to be added
        """
        branch.set_t_pmer(len(self._Network__pl) - 1)
        self.__p_branches[branch.get_s_pmer()].append(branch)


class Branch:
    """
    Class to model a brach in a Network
    """

    def __init__(self, point, s_pmer_id, s_mmer_id, t_pmer_id=None):
        """
        Constructor

        :param point: branch point coordinates (3-array)
        :param s_pmer_id: starting polymer ID
        :param s_mmer_id: starting monomer ID
        :param t_pmer_id: targeting polymer ID, (default None) it may be unknown at construction time, the targeting
                          monomer ID is 0.
        """
        assert hasattr(point, '__len__') and (len(point) == 3)
        assert (s_pmer_id >= 0) and (s_mmer_id >= 0)
        if t_pmer_id is not None:
            assert t_pmer_id >= 0
        self.__point = np.asarray(point, dtype=float)
        self.__s_pmer_id, self.__s_mmer_id, self.__t_pmer_id = s_pmer_id, s_mmer_id, t_pmer_id

    def set_t_pmer(self, t_pmer_id):
        """
        Set targeting polymer ID
        """
        assert t_pmer_id >= 0
        self.__t_pmer_id = t_pmer_id

    def get_point(self):
        """
        Get point coordinates
        """
        return self.__point

    def get_s_pmer(self):
        """
        Get starting polymer ID
        """
        return self.__s_pmer_id

    def get_s_mmer(self):
        """
        Get starting monomer ID
        """
        return self.__s_mmer_id

    def get_t_pmer(self):
        """
        Get targeting polymer ID
        """
        return self.__t_pmer_id

    def get_vtp(self, shape_vtp=None):
        """
        Gets a polydata with the branch shape

        :param shape_vtp: if None (default) the a point is returned, otherwise this shape is used
                          TODO: so far only isotropic shapes are recommended and starting monomer tangent is not considered yet
        :return:
        """
        if shape_vtp is None:
            return point_to_poly(self.get_point())
        else:
            return poly_translate(shape_vtp, self.get_point())


