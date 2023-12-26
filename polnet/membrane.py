"""
Classes for modeling membranes.
A membrane is modelled as two parallel surfaces with Gaussian profile
"""

__author__ = 'Antonio Martinez-Sanchez'

import scipy as sp
from polnet.poly import *
from polnet.affine import *
from polnet.lrandom import *
from abc import ABC, abstractmethod

MAX_TRIES_MB = 10


class MbError(Exception):
    pass


class Mb(ABC):
    """
    Abstract class to model membranes with different geometries
    """

    def __init__(self, tomo_shape, v_size=1, center=(0, 0, 0), rot_q=(1, 0, 0, 0), thick=1, layer_s=1):
        """
        Constructor
        :param tomo_shape: reference tomogram shape (X, Y and Z dimensions)
        :param v_size: reference tomogram voxel size (default 1)
        :param center: ellipsoid center (VERY IMPORTANT: coordinates are not in voxels)
        :param rot_q: rotation expressed as quaternion with respect ellipsoid center (default [1, 0, 0, 0] no rotation)
        :param thick: membrane thickness (default 1)
        :param layer_s: Gaussian sigma for each layer
        """
        assert hasattr(tomo_shape, '__len__') and (len(tomo_shape) == 3)
        assert v_size > 0
        assert (thick > 0) and (layer_s > 0)
        assert hasattr(center, '__len__') and (len(center) == 3)
        assert hasattr(rot_q, '__len__') and (len(rot_q) == 4)
        self.__tomo_shape, self.__v_size = tomo_shape, v_size
        self.__center, self.__rot_q = np.asarray(center, dtype=float), np.asarray(rot_q, dtype=float)
        self.__thick, self.__layer_s = float(thick), float(layer_s)
        self.__tomo, self.__mask, self.__surf = None, None, None

    def get_thick(self):
        """
        Get membrane thickness, bilayer gap
        :return: thickness as a float
        """
        return self.__thick

    def get_layer_s(self):
        """
        Get Gaussian sigma for each layer
        :return: layer sigma as a float
        """
        return self.__layer_s

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
        return self.__surf

    def masking(self, mask):
        """
        Removes membrane voxels in an external mask
        :param mask: the input external mask, binary ndarray with the same shape as the membrane tomogram, tomogram
        voxels at mask 0-valued positions will be set to 0
        :return: None
        """
        assert isinstance(mask, np.ndarray) and mask.dtype == bool
        assert len(mask.shape) == len(self.__tomo.shape) and mask.shape == self.__tomo.shape
        mask_ids = np.invert(mask)
        self.__tomo[mask_ids] = 0
        self.__mask[mask_ids] = False
        self.__surf = poly_mask(self.__surf, mask)


    def insert_density_svol(self, tomo, merge='max', mode='tomo', grow=0):
        """
        Insert a membrane into a tomogram
        :param tomo: tomogram where m_svol is added
        :param merge: merging mode, valid: 'min' (default), 'max', 'sum' and 'insert'
        :param mode: determines which data are inserted, valid: 'tomo' (default), 'mask' and 'voi'
        :param grow: number of voxel to grow the membrane tomogram to insert (default 0), only used in 'voi' mode
        :return:
        """
        assert (mode == 'tomo') or (mode == 'mask') or (mode == 'voi')
        if mode == 'tomo':
            hold = self.__tomo
        elif mode == 'mask':
            hold = self.__mask
        elif mode == 'voi':
            if grow >= 1:
                hold = np.invert(sp.ndimage.morphology.binary_dilation(self.__mask, iterations=grow))
            else:
                hold = np.invert(self.__mask)
        insert_svol_tomo(hold, tomo, .5 * np.asarray(self.__tomo.shape), merge=merge)

    @abstractmethod
    def __build_tomos(self):
        """
        Generates the membrane within a tomogram
        :return: the generated tomogram and its binary mask
        """
        raise NotImplemented


class MbEllipsoid(Mb):
    """
    Class for generating a membrane with Ellipsoid shape
    """

    def __init__(self, tomo_shape, v_size=1, center=(0, 0, 0), rot_q=(1, 0, 0, 0), thick=1, layer_s=1, a=1, b=1, c=1):
        """
        Constructor
        :param tomo_shape: reference tomogram shape (X, Y and Z dimensions)
        :param v_size: reference tomogram voxel size (default 1)
        :param center: ellipsoid center (VERY IMPORTANT: coordinates are not in voxels)
        :param rot_q: rotation expressed as quaternion with respect ellipsoid center (default [1, 0, 0, 0] no rotation)
        :param thick: membrane thickness (default 1)
        :param layer_s: Gaussian sigma for each layer
        :param a: (default 1) semi axis length in X axis (before rotation)
        :param b: (default 1) semi axis length in Y axis (before rotation)
        :param c: (default 1) semi axis length in Z axis (before rotation)
        """
        super(MbEllipsoid, self).__init__(tomo_shape, v_size, center, rot_q, thick, layer_s)
        assert (a > 0) and (b > 0) and (c > 0)
        self.__a, self.__b, self.__c = float(a), float(b), float(c)
        self._Mb__build_tomos()

    def _Mb__build_tomos(self):

        # Input parsing
        t_v, s_v = .5 * self._Mb__thick / self._Mb__v_size, self._Mb__layer_s / self._Mb__v_size
        a_v, b_v, c_v = self.__a / self._Mb__v_size, self.__b / self._Mb__v_size, self.__c / self._Mb__v_size
        ao_v, bo_v, co_v = a_v + t_v, b_v + t_v, c_v + t_v
        ai_v, bi_v, ci_v = a_v - t_v, b_v - t_v, c_v - t_v
        ao_v_p1, bo_v_p1, co_v_p1 = ao_v + 1, bo_v + 1, co_v + 1
        ao_v_m1, bo_v_m1, co_v_m1 = ao_v - 1, bo_v - 1, co_v - 1
        ai_v_p1, bi_v_p1, ci_v_p1 = ai_v + 1, bi_v + 1, ci_v + 1
        ai_v_m1, bi_v_m1, ci_v_m1 = ai_v - 1, bi_v - 1, ci_v - 1
        p0_v = self._Mb__center / self._Mb__v_size

        # Generating the grid
        dx, dy, dz = float(self._Mb__tomo_shape[0]), float(self._Mb__tomo_shape[1]), float(self._Mb__tomo_shape[2])
        dx2, dy2, dz2 = math.floor(.5 * dx), math.floor(.5 * dy), math.floor(.5 * dz)
        p0_v[0] -= dx2
        p0_v[1] -= dy2
        p0_v[2] -= dz2
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(np.arange(x_l, x_h), np.arange(y_l, y_h), np.arange(z_l, z_h), indexing='xy')

        # Mask generation
        R_o = ((X - p0_v[0]) / ao_v) ** 2 + ((Y - p0_v[1]) / bo_v) ** 2 + ((Z - p0_v[2]) / co_v) ** 2
        R_i = ((X - p0_v[0]) / ai_v) ** 2 + ((Y - p0_v[1]) / bi_v) ** 2 + ((Z - p0_v[2]) / ci_v) ** 2
        self._Mb__mask = tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), self._Mb__rot_q, order=0)
        if self._Mb__mask.sum() == 0:
            raise MbError

        # Surface generation
        R_i = ((X - p0_v[0]) / a_v) ** 2 + ((Y - p0_v[1]) / b_v) ** 2 + ((Z - p0_v[2]) / c_v) ** 2
        R_i = tomo_rotate(R_i, self._Mb__rot_q, mode='reflect')
        self._Mb__surf = iso_surface(R_i, 1)
        add_sfield_to_poly(self._Mb__surf, self._Mb__mask, 'mb_mask', dtype='int', interp='NN', mode='points')
        # lio.save_vtp(self._Mb__surf, './out/hold.vtp')
        self._Mb__surf = poly_threshold(self._Mb__surf, 'mb_mask', mode='points', low_th=.5)
        # lio.save_vtp(self._Mb__surf, './out/hold2.vtp')

        # Outer layer
        R_o = ((X - p0_v[0]) / ao_v_p1) ** 2 + ((Y - p0_v[1]) / bo_v_p1) ** 2 + ((Z - p0_v[2]) / co_v_p1) ** 2
        R_i = ((X - p0_v[0]) / ao_v_m1) ** 2 + ((Y - p0_v[1]) / bo_v_m1) ** 2 + ((Z - p0_v[2]) / co_v_m1) ** 2
        G = tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), self._Mb__rot_q, order=0)

        # Inner layer
        R_o = ((X - p0_v[0]) / ai_v_p1) ** 2 + ((Y - p0_v[1]) / bi_v_p1) ** 2 + ((Z - p0_v[2]) / ci_v_p1) ** 2
        R_i = ((X - p0_v[0]) / ai_v_m1) ** 2 + ((Y - p0_v[1]) / bi_v_m1) ** 2 + ((Z - p0_v[2]) / ci_v_m1) ** 2
        G += tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), self._Mb__rot_q, order=0)

        # Smoothing
        self._Mb__tomo = lin_map(density_norm(sp.ndimage.gaussian_filter(G.astype(float), s_v), inv=True), ub=0, lb=1)


class MbSphere(Mb):
    """
    Class for generating a membrane with Spherical shape
    """

    def __init__(self, tomo_shape, v_size=1, center=(0, 0, 0), rot_q=(1, 0, 0, 0), thick=1, layer_s=1, rad=1):
        """
        Constructor
        :param tomo_shape: reference tomogram shape (X, Y and Z dimensions)
        :param v_size: reference tomogram voxel size (default 1)
        :param center: ellipsoid center (VERY IMPORTANT: coordinates are not in voxels)
        :param rot_q: rotation expressed as quaternion with respect ellipsoid center (default [1, 0, 0, 0] no rotation)
        :param thick: membrane thickness (default 1)
        :param layer_s: Gaussian sigma for each layer
        :param rad: (default 1) sphere radius
        """
        super(MbSphere, self).__init__(tomo_shape, v_size, center, rot_q, thick, layer_s)
        assert (rad > 0)
        self.__rad = float(rad)
        self._Mb__build_tomos()

    def _Mb__build_tomos(self):

        # Input parsing
        t_v, s_v = .5 * (self._Mb__thick / self._Mb__v_size), self._Mb__layer_s / self._Mb__v_size
        rad_v = self.__rad / self._Mb__v_size
        ao_v = rad_v + t_v
        ai_v = rad_v - t_v
        ao_v_p1 = ao_v + 1
        ao_v_m1 = ao_v - 1
        ai_v_p1 = ai_v + 1
        ai_v_m1 = ai_v - 1
        p0_v = self._Mb__center / self._Mb__v_size

        # Generating the bilayer
        dx, dy, dz = float(self._Mb__tomo_shape[0]), float(self._Mb__tomo_shape[1]), float(self._Mb__tomo_shape[2])
        dx2, dy2, dz2 = math.floor(.5 * dx), math.floor(.5 * dy), math.floor(.5 * dz)
        p0_v[0] -= dx2
        p0_v[1] -= dy2
        p0_v[2] -= dz2
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(np.arange(x_l, x_h), np.arange(y_l, y_h), np.arange(z_l, z_h), indexing='xy')

        # Mask generation
        R_o = ((X - p0_v[0]) / ao_v) ** 2 + ((Y - p0_v[1]) / ao_v) ** 2 + ((Z - p0_v[2]) / ao_v) ** 2
        R_i = ((X - p0_v[0]) / ai_v) ** 2 + ((Y - p0_v[1]) / ai_v) ** 2 + ((Z - p0_v[2]) / ai_v) ** 2
        self._Mb__mask = tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), self._Mb__rot_q, order=0)

        # Surface generation
        R_i = ((X - p0_v[0]) / rad_v) ** 2 + ((Y - p0_v[1]) / rad_v) ** 2 + ((Z - p0_v[2]) / rad_v) ** 2
        R_i = tomo_rotate(R_i, self._Mb__rot_q, mode='reflect')
        self._Mb__surf = iso_surface(R_i, 1)
        add_sfield_to_poly(self._Mb__surf, self._Mb__mask, 'mb_mask', dtype='int', interp='NN', mode='points')
        self._Mb__surf = poly_threshold(self._Mb__surf, 'mb_mask', mode='points', low_th=.5)

        # Outer layer
        R_o = ((X - p0_v[0]) / ao_v_p1) ** 2 + ((Y - p0_v[1]) / ao_v_p1) ** 2 + ((Z - p0_v[2]) / ao_v_p1) ** 2
        R_i = ((X - p0_v[0]) / ao_v_m1) ** 2 + ((Y - p0_v[1]) / ao_v_m1) ** 2 + ((Z - p0_v[2]) / ao_v_m1) ** 2
        G = tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), self._Mb__rot_q, order=0)
        # R = (X - p0_v[0])**2 + (Y - p0_v[1])**2 + (Z - p0_v[2])**2
        # G = tomo_rotate(np.logical_and(R >= ao_v_m1**2, R <= ao_v_p1**2), self._Mb__rot_q, order=0)

        # Inner layer
        R_o = ((X - p0_v[0]) / ai_v_p1) ** 2 + ((Y - p0_v[1]) / ai_v_p1) ** 2 + ((Z - p0_v[2]) / ai_v_p1) ** 2
        R_i = ((X - p0_v[0]) / ai_v_m1) ** 2 + ((Y - p0_v[1]) / ai_v_m1) ** 2 + ((Z - p0_v[2]) / ai_v_m1) ** 2
        G += tomo_rotate(np.logical_and(R_i >= 1, R_o <= 1), self._Mb__rot_q, order=0)
        # G += tomo_rotate(np.logical_and(R >= ai_v_m1**2, R_o <= ai_v_p1**2), self._Mb__rot_q, order=0)

        # Smoothing
        # TODO: is it required the density_norm() having lin_map()?
        # self._Mb__tomo = lin_map(density_norm(sp.ndimage.gaussian_filter(G.astype(float), s_v), inv=True), ub=0, lb=1)
        self._Mb__tomo = lin_map(-1 * sp.ndimage.gaussian_filter(G.astype(float), s_v), ub=0, lb=1)


class MbTorus(Mb):
    """
    Class for generating a membrane with Toroidal shape
    """

    def __init__(self, tomo_shape, v_size=1, center=(0, 0, 0), rot_q=(1, 0, 0, 0), thick=1, layer_s=1, rad_a=1, rad_b=1):
        """
        Constructor
        :param tomo_shape: reference tomogram shape (X, Y and Z dimensions)
        :param v_size: reference tomogram voxel size (default 1)
        :param center: ellipsoid center (VERY IMPORTANT: coordinates are not in voxels)
        :param rot_q: rotation expressed as quaternion with respect ellipsoid center (default [1, 0, 0, 0] no rotation)
        :param thick: membrane thickness (default 1)
        :param layer_s: Gaussian sigma for each layer
        :param rad_a: (default 1) torus radius
        :param rad_b: (default 1) torus tube radius
        """
        super(MbTorus, self).__init__(tomo_shape, v_size, center, rot_q, thick, layer_s)
        assert (rad_a > 0) and (rad_b > 0)
        self.__rad_a, self.__rad_b = float(rad_a), float(rad_b)
        self._Mb__build_tomos()

    def _Mb__build_tomos(self):

        # Input parsing
        t_v, s_v = .5 * (self._Mb__thick / self._Mb__v_size), self._Mb__layer_s / self._Mb__v_size
        rad_a_v, rad_b_v = self.__rad_a / self._Mb__v_size, self.__rad_b / self._Mb__v_size
        bo_v, bi_v = rad_b_v + t_v, rad_b_v - t_v
        bo_v_p1, bo_v_m1 = bo_v + 1, bo_v - 1
        bi_v_p1, bi_v_m1 = bi_v + 1, bi_v - 1
        p0_v = self._Mb__center / self._Mb__v_size

        # Generating the bilayer
        dx, dy, dz = float(self._Mb__tomo_shape[0]), float(self._Mb__tomo_shape[1]), float(self._Mb__tomo_shape[2])
        dx2, dy2, dz2 = math.floor(.5 * dx), math.floor(.5 * dy), math.floor(.5 * dz)
        p0_v[0] -= dx2
        p0_v[1] -= dy2
        p0_v[2] -= dz2
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(np.arange(x_l, x_h), np.arange(y_l, y_h), np.arange(z_l, z_h), indexing='xy')

        # Mask generation
        R_o = ((rad_a_v - np.sqrt((X-p0_v[0])**2 + (Y-p0_v[1])**2))**2 + (Z-p0_v[2])**2 - bo_v*bo_v) <= 1
        R_i = ((rad_a_v - np.sqrt((X-p0_v[0])**2 + (Y-p0_v[1])**2))**2 + (Z-p0_v[2])**2 - bi_v*bi_v) >= 1
        self._Mb__mask = tomo_rotate(np.logical_and(R_i, R_o), self._Mb__rot_q, order=0)

        # Surface generation
        R_i = ((rad_a_v - np.sqrt((X-p0_v[0])**2 + (Y-p0_v[1])**2))**2 + (Z-p0_v[2])**2 - rad_b_v*rad_b_v)
        R_i = tomo_rotate(R_i, self._Mb__rot_q, mode='reflect')
        self._Mb__surf = iso_surface(R_i, 1)
        add_sfield_to_poly(self._Mb__surf, self._Mb__mask, 'mb_mask', dtype='int', interp='NN', mode='points')
        self._Mb__surf = poly_threshold(self._Mb__surf, 'mb_mask', mode='points', low_th=.5)

        # Outer layer
        R_o = ((rad_a_v - np.sqrt((X-p0_v[0])**2 + (Y-p0_v[1])**2))**2 + (Z-p0_v[2])**2 - bo_v_p1*bo_v_p1) <= 1
        R_i = ((rad_a_v - np.sqrt((X-p0_v[0])**2 + (Y-p0_v[1])**2))**2 + (Z-p0_v[2])**2 - bo_v_m1*bo_v_m1) >= 1
        G = tomo_rotate(np.logical_and(R_i, R_o), self._Mb__rot_q, order=0)

        # Inner layer
        R_o = ((rad_a_v - np.sqrt((X-p0_v[0])**2 + (Y-p0_v[1])**2))**2 + (Z-p0_v[2])**2 - bi_v_p1 * bi_v_p1) <= 1
        R_i = ((rad_a_v - np.sqrt((X-p0_v[0])**2 + (Y-p0_v[1])**2))**2 + (Z-p0_v[2])**2 - bi_v_m1 * bi_v_m1) >= 1
        G += tomo_rotate(np.logical_and(R_i, R_o), self._Mb__rot_q, order=0)

        # Smoothing
        self._Mb__tomo = lin_map(density_norm(sp.ndimage.gaussian_filter(G.astype(float), s_v), inv=True), ub=0, lb=1)


class SetMembranes:
    """
    Class for modelling a set of membranes within a tomogram
    """

    def __init__(self, voi, v_size, gen_rnd_surfs, param_rg, thick_rg, layer_rg, occ, over_tolerance=0, bg_voi=None,
                 grow=0):
        """
        Construction
        :param voi: a 3D numpy array to define a VOI (Volume Of Interest) for membranes
        :param v_size: voxel size
        :param gen_rnd_surf: an of object that inherits from lrandom.SurfGen class to generate random instances with
                             membrane surface parameters, therefore the objects class determine the shape of the membranes
                             generated
        :param thick_rg: membrane thickness range (2-tuple)
        :param layer_s: lipid layer range (2-tuple)
        :param occ: occupancy threshold in percentage [0, 100]%
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0, in range [0,1))
        :param bg_voi: background VOI (Default None), if present membrane areas in this VOI will be removed and
        :param grow: number of voxel to grow the VOI
        not considered for overlapping. It must be binary with the same shape of voi
        """

        # Input parsing
        assert isinstance(voi, np.ndarray) and (voi.dtype == bool)
        assert issubclass(gen_rnd_surfs.__class__, SurfGen)
        assert hasattr(param_rg, '__len__') and (len(param_rg) == 3) and (param_rg[0] <= param_rg[1])
        assert hasattr(thick_rg, '__len__') and (len(thick_rg) == 2) and (thick_rg[0] <= thick_rg[1])
        assert hasattr(layer_rg, '__len__') and (len(layer_rg) == 2) and (layer_rg[0] <= layer_rg[1])
        assert (occ >= 0) and (occ <= 100)
        assert (over_tolerance >= 0) and (over_tolerance <= 100)
        if bg_voi is not None:
            assert isinstance(bg_voi, np.ndarray) and (bg_voi.dtype == bool)
            assert len(voi.shape) == len(bg_voi.shape) and voi.shape == bg_voi.shape
        assert grow >= 0

        # Variables assignment
        self.__voi = voi
        if bg_voi is not None:
            self.__bg_voi = bg_voi
        self.__vol = float(self.__voi.sum()) * v_size * v_size * v_size # without the float cast it may raise overflow warining in Windows
        self.__v_size = v_size
        self.__tomo, self.__gtruth = np.zeros(shape=voi.shape, dtype=np.float16), \
                                     np.zeros(shape=voi.shape, dtype=bool)
        self.__surfs, self.__app_vtp = vtk.vtkPolyData(), vtk.vtkAppendPolyData()
        self.__count_mbs = 0
        self.__gen_rnd_surfs = gen_rnd_surfs
        self.__param_rg, self.__thick_rg, self.__layer_rg = param_rg, thick_rg, layer_rg
        self.__occ, self.__over_tolerance = occ, over_tolerance
        self.__grow = grow

    def get_vol(self):
        return self.__vol

    def get_mb_occupancy(self):
        return self.__gtruth.sum() / np.prod(np.asarray(self.__voi.shape, dtype=float))

    def build_set(self, verbosity=False):
        """
        Build a set of ellipsoid membranes and insert them in a tomogram and a vtkPolyData object
        :param verbosity: if True (default False) the output message with info about the membranes generated is printed
        :return:
        """

        # Initialization
        count_mb = 1
        count_exp = 0

        # Network loop
        while self.get_mb_occupancy() < self.__occ:

            # Polymer initialization
            p0 = np.asarray((self.__voi.shape[0] * self.__v_size * random.random(),
                             self.__voi.shape[1] * self.__v_size * random.random(),
                             self.__voi.shape[2] * self.__v_size * random.random()))
            thick, layer_s = random.uniform(self.__thick_rg[0], self.__thick_rg[1]), \
                             random.uniform(self.__layer_rg[0], self.__layer_rg[1])
            rot_q = gen_rand_unit_quaternion()

            try:

                # Membrane generation according the predefined surface model
                if isinstance(self.__gen_rnd_surfs, EllipGen):
                    ellip_axes = self.__gen_rnd_surfs.gen_parameters_exp()
                    hold_mb = MbEllipsoid(self.__voi.shape, v_size=self.__v_size,
                                          center=p0, rot_q=rot_q, thick=thick, layer_s=layer_s,
                                          a=ellip_axes[0], b=ellip_axes[1], c=ellip_axes[2])
                    hold_rad = np.mean(ellip_axes)
                elif isinstance(self.__gen_rnd_surfs, SphGen):
                    rad = self.__gen_rnd_surfs.gen_parameters()
                    hold_mb = MbSphere(self.__voi.shape, v_size=self.__v_size,
                                       center=p0, rot_q=rot_q, thick=thick, layer_s=layer_s, rad=rad)
                    hold_rad = rad
                elif isinstance(self.__gen_rnd_surfs, TorGen):
                    tor_axes = self.__gen_rnd_surfs.gen_parameters()
                    hold_mb = MbTorus(self.__voi.shape, v_size=self.__v_size,
                                      center=p0, rot_q=rot_q, thick=thick, layer_s=layer_s,
                                      rad_a=tor_axes[0], rad_b=tor_axes[1])
                    hold_rad = np.mean(tor_axes)
                else:
                    print('ERROR: not valid random surface parameters generator: ' + str(self.__gen_rnd_surfs.__class__))
                    raise MbError

                # Background masking
                if self.__bg_voi is not None:
                    hold_mb.masking(self.__bg_voi)

                # Insert membrane
                self.insert_mb(hold_mb, merge='max', over_tolerance=self.__over_tolerance, check_vol=True, grow=self.__grow)
                if verbosity:
                    print('Membrane ' + str(count_mb) + ', total occupancy: ' + str(self.get_mb_occupancy()) +
                          ', volume: ' + str(hold_mb.get_vol()) + ', thickness: ' + str(hold_mb.get_thick()) +
                          ', layer_s: ' + str(hold_mb.get_layer_s()) + ', radius (avg): ' + str(hold_rad))
                count_mb += 1

            # Handling the exception raised when a membrane could not be generated properly
            except MbError:
                count_exp += 1
                # print('JOl')
                # print('Count: ' + str(count_exp))
                if count_exp == MAX_TRIES_MB:
                    print('WARNING: more than ' + str(MAX_TRIES_MB) + ' tries failed to insert a membrane!')
                    break
            count_exp = 0

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
        # return np.invert(self.__voi) * self.__tomo
        return self.__tomo

    def get_gtruth(self):
        """
        Get the ground truth within the VOI
        :return: an ndarray
        """
        # return self.__voi * self.__gtruth
        return self.__gtruth

    def get_vtp(self):
        """
        Get the set of membranes as a vtkPolyData with their surfaces
        :return: a vtkPolyData
        """
        return self.__surfs

    def get_num_mbs(self):
        """
        Get the number of membranes in the set
        :return: an integer with the number of membranes
        """
        return self.__count_mbs

    def check_overlap(self, mb, over_tolerance):
        """
        Determines if the membrane overlaps with any member within the membranes set
        :param mb: input Membrane to check for the overlapping
        :param over_tolerance: overlapping tolerance (percentage of membrane voxel overlapping)
        """
        mb_mask = mb.get_mask()
        #tomo_mb = np.zeros(shape=mb_mask.shape, dtype=bool)
        # mb.insert_density_svol(tomo_mb, merge='max')
        # tomo_over = np.logical_and(mb_mask, self.__gtruth)
        tomo_over = np.logical_and(mb_mask, np.logical_not(self.__voi))
        if 100. * (tomo_over.sum() / self.get_vol()) > over_tolerance:
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

    def insert_mb(self, mb, merge='min', over_tolerance=None, check_vol=True, grow=0):
        """
        Insert the membrane into the set (tomogram, vtkPolyData and Ground Truth)
        :param mb: input membrane (Mb) object
        :param merge: merging mode for density insertion, valid: 'min' (default), 'max', 'sum' and 'insert'
        :param over_tolerance: overlapping tolerance (percentage of membrane voxel overlapping), if None then disabled
        :param check_vol: if True (default) the input membrane volume is check and if equal to zero then the membrane
                          is not inserted and MbError is raised
        :param grow: number of voxel to grow the membrane tomogram to insert (default 0), only used in 'voi' mode
        :return: raises a ValueError if the membrane is not inserted
        """
        if check_vol and (mb.get_vol() <= 0):
            raise MbError
        if (over_tolerance is None) or (not self.check_overlap(mb, over_tolerance)):
            # Density tomogram insertion
            mb.insert_density_svol(self.__tomo, merge=merge, mode='tomo')
            # Ground Truth
            mb.insert_density_svol(self.__gtruth, merge='max', mode='mask')
            # VOI
            mb.insert_density_svol(self.__voi, merge='min', mode='voi', grow=grow)
            # Surfaces insertion
            self.__app_vtp.AddInputData(mb.get_vtp())
            self.__app_vtp.Update()
            self.__surfs = poly_scale(self.__app_vtp.GetOutput(), self.__v_size)
            self.__count_mbs += 1
        else:
            raise MbError
