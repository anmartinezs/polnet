"""
Classes for modeling membranes.
A membrane is modelled as two parallel surfaces with Gaussian profile
"""

__author__ = 'Antonio Martinez-Sanchez'

from polnet.utils import *
from polnet.affine import *
from polnet.lrandom import EllipGen
from abc import ABC, abstractmethod


class MbEllipsoid:
    """
    Class for generating a membrane with Ellipsoid shape
    """

    def __init__(self, a, b, c, tomo_shape, v_size=1, center=(0, 0, 0), rot_q=(1, 0, 0, 0), thick=1, layer_s=1):
        """
        Constructor
        :param a: semi axis length in X axis (before rotation)
        :param b: semi axis length in Y axis (before rotation)
        :param c: semi axis length in Z axis (before rotation)
        :param tomo_shape: reference tomogram shape (X, Y and Z dimensions)
        :param v_size: reference tomogram voxel size (default 1)
        :param center: ellipsoid center (VERY IMPORTANT: coordinates are not in voxels)
        :param rot_q: rotation expressed as quaternion with respect ellipsoid center (default [1, 0, 0, 0] no rotation)
        :param thick: membrane thickness (default 1)
        :param layer_s: Gaussian sigma for each layer
        """
        assert hasattr(tomo_shape, '__len__') and (len(tomo_shape) == 3)
        assert v_size > 0
        assert (a > 0) and (b > 0) and (c > 0)
        assert hasattr(center, '__len__') and (len(center) == 3)
        assert hasattr(rot_q, '__len__') and (len(rot_q) == 4)
        assert (thick > 0) and (layer_s > 0)
        self.__tomo_shape, self.__v_size = tomo_shape, v_size
        self.__a, self.__b, self.__c = float(a), float(b), float(c)
        self.__center, self.__rot_q = np.asarray(center, dtype=float), np.asarray(rot_q, dtype=float)
        self.__thick, self.__layer_s = thick, layer_s
        self.__tomos, self.__mask = None, None
        self.__build_tomos()

    def get_vol(self):
        """
        Get the polymer volume
        :param fast: if True (default) the volume monomer is only computed once
        :return: the computed volume
        """
        return self.__mask.sum() * self.__v_size**3

    def get_tomo(self):
        """
        Get the membrane within a tomogram
        :return: a numpy 3D array
        """
        return self.__tomo

    def get_mask(self):
        """
        Get the membrane within a binary tomogram
        :return: a binary numpy 3D array
        """
        return self.__mask

    def get_vtp(self):
        """
        Get the membrane as an VTK surface
        :return: a vtkPolyData object
        """
        return iso_surface(self.__mask.astype(np.float), .5)

    def __build_tomos(self):
        """
        Generates a tomogram that encloses the membrane
        :return: the generated tomogram and its binary mask
        """

        # Input parsing
        t_v, s_v = .5 * (self.__thick / self.__v_size), self.__layer_s / self.__v_size
        a_v, b_v, c_v = self.__a / self.__v_size, self.__b / self.__v_size, self.__c / self.__v_size
        ao_v, bo_v, co_v = a_v + t_v, b_v + t_v, c_v + t_v
        ai_v, bi_v, ci_v = a_v - t_v, b_v - t_v, c_v - t_v
        g_cte = 2 * s_v * s_v

        # Generating the bilayer
        dx, dy, dz = float(self.__tomo_shape[0]), float(self.__tomo_shape[1]), float(self.__tomo_shape[2])
        dx2, dy2, dz2 = math.floor(.5 * dx), math.floor(.5 * dy), math.floor(.5 * dz)
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(np.arange(x_l, x_h), np.arange(y_l, y_h), np.arange(z_l, z_h), indexing='xy')
        # X, Y, Z = X.astype(np.float16), Y.astype(np.float16), X.astype(np.float16)
        X *= self.__v_size
        Y *= self.__v_size
        Z *= self.__v_size
        R_o = np.sqrt(X**2/ao_v**2 + Y**2/bo_v**2 + Z**2/co_v**2)
        G_o = np.exp(-(R_o - 1) ** 2 / g_cte)
        R_i = np.sqrt(X ** 2 / ai_v ** 2 + Y ** 2 / bi_v ** 2 + Z ** 2 / ci_v ** 2)
        G_i = np.exp(-(R_o - 1) ** 2 / g_cte)

        # Tomogram and binary mask
        self.__tomo = np.asarray(G_i + G_o, dtype=float)
        self.__mask = (R_i >= 1) + (R_o <= 1)

        # Tomograms rotation
        self.__tomo = tomo_rotate(self.__tomo, self.__rot_q)
        self.__mask = tomo_rotate(self.__mask, self.__rot_q, order=1)

    def insert_density_svol(self, tomo, merge='min'):
        """
        Insert a the membrane into a tomogram
        :param tomo: tomogram where m_svol is added
        :param merge: merging mode, valid: 'min' (default), 'max', 'sum' and 'insert'
        :return:
        """
        insert_svol_tomo(self.get_tomo(), tomo, self.__center / self.__v_size, merge=merge)


class SetMembranes(ABC):
    """
    Class for modelling a set of membranes within a tomogram
    """

    def __init__(self, voi, v_size):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for membranes
        :param v_size: voxel size
        """
        assert isinstance(voi, np.ndarray)
        self.__voi = voi
        self.__vol = (self.__voi > 0).sum() * v_size * v_size * v_size
        self.__v_size = v_size
        self.__tomo, self.__gtruth = np.zeros(shape=voi.shape, dtype=np.float16), \
                                     np.zeros(shape=voi.shape, dtype=bool)
        self.__surfs, self.__app_vtp = vtk.vtkPolyData(), vtk.vtkAppendPolyData()

    def get_mb_occupancy(self):
        return self.__gtruth.sum() / np.prod(np.asarray(self.__voi.shape, dtype=float))

    @abstractmethod
    def build_set(self):
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

    def get_tomo(self):
        """
        Get the tomogram with the membranes within the VOI
        :return: an ndarray
        """
        return self.__voi * self.__tomo

    def get_gtruth(self):
        """
        Get the ground truth within the VOI
        :return: an ndarray
        """
        return self.__voi * self.__gtruth

    def get_vtp(self):
        """
        Get the set of membranes as a vtkPolyData with their surfaces
        :return: a vtkPolyData
        """
        return self.__surfs

    def check_overlap(self, mb, over_tolerance):
        """
        Determines if the membrane overlaps with any within the membranes set
        :param mb: input Membrane to check for the overlapping
        :param over_tolerance: overlapping tolerance (percentage of membrane voxel overlapping)
        """
        mb_mask = mb.get_mask()
        tomo_mb = np.zeros(shape=mb_mask.shape, dtype=bool)
        mb.insert_density_svol(tomo_mb, merge='max')
        tomo_over = np.logical_and(np.logical_and(tomo_mb, self.__gtruth), self.__voi)
        if 100. * (tomo_over.sum() / self.get_vol()) > over_tolerance:
            return True
        return False

    def check_overlap(self, mb, over_tolerance):
        """
        Determines if the membrane overlaps with any within the membranes set
        :param mb: input Membrane to check for the overlapping
        :param over_tolerance: overlapping tolerance (percentage of membrane voxel overlapping)
        """
        if self.compute_overlap(mb) > over_tolerance:
            return True
        return False

    def compute_overlap(self, mb):
        """
        Computes membrane overlapping with the set
        :param mb: input Membrane to check for the overlapping
        """
        mb_mask = mb.get_mask()
        tomo_mb = np.zeros(shape=mb_mask.shape, dtype=bool)
        mb.insert_density_svol(tomo_mb, merge='max')
        tomo_over = np.logical_and(np.logical_and(tomo_mb, self.__gtruth), self.__voi)
        return 100. * (tomo_over.sum() / self.get_vol())

    def insert_mb(self, mb, merge='min', over_tolerance=None):
        """
        Insert the membrane into the set (tomogram, vtkPolyData and Ground Truth)
        :param mb: input membrane (Mb) object
        :param merge: merging mode for density insertion, valid: 'min' (default), 'max', 'sum' and 'insert'
        :param over_tolerance: overlapping tolerance (percentage of membrane voxel overlapping), if None then disabled
        :return: raises a ValueError if the membrane is not inserted
        """
        if (over_tolerance is None) or (not self.check_overlap(mb, over_tolerance)):
            # Density tomogram insertion
            mb.insert_density_svol(self.__tomo, merge=merge)
            # Ground Truth and VOI insertion
            mb.insert_density_svol(self.__voi, merge='max')
            mb.insert_density_svol(self.__gtruth, merge='max')
            # Surfaces insertion
            self.__app_vtp.AddInputData(mb.get_vtp())
            self.__app_vtp.Update()
            self.__surfs = self.__app_vtp.GetOutput()
        else:
            raise ValueError


class SetEllipMembranes(SetMembranes):
    """
    Class for generating a set of ellipsoid membranes
    """

    def __init__(self, voi, v_size, gen_ellip_mbs, param_rg, thick_rg, layer_rg, occ, over_tolerance=0):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for polymers
        :param v_size: voxel size (default 1)
        :param gen_ellip_mbs: a instance of a random generation model (random.MbGen) to determine the geometrical
        paramters of a membrane
        :param param_rg: ragen for the ellipsoid parameters (the three radii)
        :param thick_rg: membrane thickness range (2-tuple)
        :param layer_s: lipid layer range (2-tuple)
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        """

        # Initialize abstract varibles
        super(SetEllipMembranes, self).__init__(voi, v_size)

        # Input parsing
        assert issubclass(gen_ellip_mbs.__class__, EllipGen)
        assert hasattr(param_rg, '__len__') and (len(param_rg) == 2) and (param_rg[0] <= param_rg[1])
        assert hasattr(thick_rg, '__len__') and (len(thick_rg) == 2) and (thick_rg[0] <= thick_rg[1])
        assert hasattr(layer_rg, '__len__') and (len(layer_rg) == 2) and (layer_rg[0] <= layer_rg[1])
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)

        # Variables assignment
        self.__gen_ellip_mbs = gen_ellip_mbs
        self.__param_rg, self.__thick_rg, self.__layer_rg = param_rg, thick_rg, layer_rg
        self.__occ, self.__over_tolerance = occ, over_tolerance

    def build_set(self):
        """
        Build a set of ellipsoid membranes and insert them in a tomogram and a vtkPolyData object
        :return:
        """

        # Initialization
        cont_mb = 1

        # Network loop
        while self.get_mb_occupancy() < self.__occ:

            # Polymer initialization
            p0 = np.asarray((self._SetMembranes__voi.shape[0] * self._SetMembranes__v_size * random.random(),
                             self._SetMembranes__voi.shape[1] * self._SetMembranes__v_size * random.random(),
                             self._SetMembranes__voi.shape[2] * self._SetMembranes__v_size * random.random()))
            a, b, c = self.__gen_ellip_mbs.gen_parameters(self.__param_rg)
            thick, layer_s = random.uniform(self.__thick_rg[0], self.__thick_rg[1]), \
                             random.uniform(self.__layer_rg[0], self.__layer_rg[1])
            rot_q = gen_rand_unit_quaternion()
            max_a = max((2. * a + 3. * (thick + layer_s), self._SetMembranes__voi.shape[0]))
            max_b = max((2. * b + 3. * (thick + layer_s), self._SetMembranes__voi.shape[1]))
            max_c = max((2. * c + 3. * (thick + layer_s), self._SetMembranes__voi.shape[2]))
            hold_mb = MbEllipsoid(max_a, max_b, max_c, self._SetMembranes__voi.shape, v_size=self._SetMembranes__v_size,
                                  center=p0, rot_q=rot_q, thick=thick, layer_s=layer_s)

            # Insert membrane
            self.insert_mb(hold_mb, merge='min', over_tolerance=self.__over_tolerance)

            print('Membrane ' + str(cont_mb) + ', total occupancy: ' + self.get_mb_occupancy())
            cont_mb += 1
