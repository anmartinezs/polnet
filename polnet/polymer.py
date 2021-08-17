"""
Classes for modeling Self-Avoiding Worm-Like Chain (SAWLC) polymers.
A polymer is a linear sequence of monomers
"""

__author__ = 'Antonio Martinez-Sanchez'

import copy
import math
import random

import numpy as np

from .utils import *
from .affine import *
from abc import ABC, abstractmethod


class Monomer:
    """
    Class for a single monomer
    """

    def __init__(self, m_surf, diam):
        """
        Constructor
        :param m_surf: monomer surface (as vtkPolyData object)
        :param diam: monomer diameter
        """
        assert isinstance(m_surf, vtk.vtkPolyData)
        assert diam > 0
        self.__m_surf = vtk.vtkPolyData()
        self.__m_surf.DeepCopy(m_surf)
        self.__diam = diam
        self.__rot_angs = np.asarray((0., 0., 0.))
        self.__bounds = np.zeros(shape=6)
        self.compute_bounds()
        # Ordered transformation queue, each entry is a 2-tuple
        # (str in ['t', 'r'], transform value in [vector, quaternion])
        self.__trans = list()

    def get_vtp(self):
        return self.__m_surf

    def get_center_mass(self):
        """
        Computer and return the monomer center of mass
        :return: a numpy array
        """
        cm_flt = vtk.vtkCenterOfMass()
        cm_flt.SetInputData(self.__m_surf)
        cm_flt.Update()
        return np.asarray(cm_flt.GetCenter())

    def get_diameter(self):
        return self.__diam

    def compute_bounds(self):
        # Compute bounds
        arr = self.__m_surf.GetPoints().GetData()
        self.__bounds[0], self.__bounds[1] = arr.GetRange(0)
        self.__bounds[2], self.__bounds[3] = arr.GetRange(1)
        self.__bounds[4], self.__bounds[5] = arr.GetRange(2)

    def rotate_q(self, q):
        """
        Applies rotation rigid transformation around center from an input unit quaternion.
        :param q: input quaternion
        :return:
        """
        w, v_axis = quat_to_angle_axis(q[0], q[1], q[2], q[3])
        self.__m_surf = poly_rotate_wxyz(self.__m_surf, w, v_axis[0], v_axis[1], v_axis[2])
        self.compute_bounds()
        self.__trans.append(('r', q))

    def translate(self, t_v):
        """
        Applies rotation rigid transformation.
        :param t_v: translation vector (x, y, z)
        :return:
        """
        self.__m_surf = poly_translate(self.__m_surf, t_v)
        self.compute_bounds()
        self.__trans.append(('t', t_v))

    def point_in_bounds(self, point):
        x_over, y_over, z_over = True, True, True
        if (self.__bounds[0] > point[0]) or (self.__bounds[1] < point[0]):
            x_over = False
        if (self.__bounds[2] > point[1]) or (self.__bounds[3] < point[1]):
            y_over = False
        if (self.__bounds[4] > point[2]) or (self.__bounds[5] < point[2]):
            y_over = False
        return x_over and y_over and z_over

    def bound_in_bounds(self, bounds):
        """
        Check if the object's bound are at least partially in another bound
        :param bounds: input bound
        :return:
        """
        x_over, y_over, z_over = True, True, True
        if (self.__bounds[0] > bounds[1]) or (self.__bounds[1] < bounds[0]):
            x_over = False
        if (self.__bounds[2] > bounds[3]) or (self.__bounds[3] < bounds[2]):
            y_over = False
        if (self.__bounds[4] > bounds[5]) or (self.__bounds[5] < bounds[4]):
            y_over = False
        return x_over and y_over and z_over

    def overlap_voi(self, voi, v_size=1):
        """
        Determines if the monomer overlaps a VOI, that requires the next condition:
            - Any particle on the monomer surface is within the VOI
        :param voi: input VOI (Volume Of Interest), binary tomgram with True for VOI voxels
        :param v_size: voxel size, it must greater than 0 (default 1)
        :return: True if the monomer overlaps the VOI, False otherwise
        """

        # Initialization
        assert v_size > 0
        assert isinstance(voi, np.ndarray) and (voi.dtype == 'bool')
        nx, ny, nz = voi.shape
        v_size_i = 1. / v_size

        # Any particle on the monomer surface is within the VOI
        for i in range(self.__m_surf.GetNumberOfPoints()):
            pt = np.asarray(self.__m_surf.GetPoint(i)) * v_size_i
            x, y, z = np.round(pt).astype(int)
            if (x < nx) and (y < ny) and (z < nz) and (x >= 0) and (y >= 0) and (z >= 0):
                if not voi[x, y, z]:
                    return True
            else:
                return True

        return False

    def get_vol(self):
        """
        Get the polymer volume
        :param fast: if True (default) the volume monomer is only computed once
        :return: the computed volume
        """
        mass = vtk.vtkMassProperties()
        mass.SetInputData(self.__m_surf)
        return mass.GetVolume()

    def get_copy(self):
        """
        Get a copy of the current Monomer
        :return: a new instance of this monomer
        """
        return Monomer(self.__m_surf, self.__diam)

    def insert_density_svol(self, m_svol, tomo, v_size=1, merge='min'):
        """
        Insert a monomer subvolume into a tomogram
        :param m_svol: input monomer sub-volume
        :param tomo: tomogram where m_svol is added
        :param v_size: tomogram voxel size (default 1)
        :param merge: merging mode, valid: 'max' (default), 'min', 'sum' and 'insert'
        :return:
        """
        v_size_i = 1. / v_size
        tot_v = np.asarray((0., 0., 0.))
        hold_svol = m_svol
        for trans in self.__trans:
            if trans[0] == 't':
                tot_v += (trans[1] * v_size_i)
            elif trans[0] == 'r':
                hold_svol = tomo_rotate(hold_svol, trans[1], mode='constant', cval=hold_svol.max())
        insert_svol_tomo(hold_svol, tomo, tot_v, merge=merge)

    def overlap_mmer(self, mmer, over_tolerance=0):
        """
        Determines if the monomer overlaps with another
        :param mmer: input monomer to check overlap with self
        :return: True if overlapping, otherwise False
        """
        # Initialization
        selector = vtk.vtkSelectEnclosedPoints()
        selector.SetTolerance(VTK_RAY_TOLERANCE)
        selector.Initialize(self.get_vtp())

        dist = points_distance(self.get_center_mass(), mmer.get_center_mass())
        if dist <= self.get_diameter():
            poly_b = mmer.get_vtp()
            count, n_points = 0., poly_b.GetNumberOfPoints()
            n_points_if = 1. / float(n_points)
            for i in range(n_points):
                if selector.IsInsideSurface(poly_b.GetPoint(i)) > 0:
                    count += 1
                    over = count * n_points_if
                    if over > over_tolerance:
                        return True

        return False

    def overlap_net(self, net, over_tolerance=0):
        """
        Determines if the monomer overlaps with another momonmer in a network
        :param mmer: input monomer to check overlap with self
        :return: True if overlapping, otherwise False
        """
        # Initialization
        selector = vtk.vtkSelectEnclosedPoints()
        selector.SetTolerance(VTK_RAY_TOLERANCE)
        selector.Initialize(self.get_vtp())

        for pmer in net.get_pmers_list():
            for mmer in pmer.get_mmers_list():
                dist = points_distance(self.get_center_mass(), mmer.get_center_mass())
                if dist <= self.get_diameter():
                    poly_b = mmer.get_vtp()
                    count, n_points = 0., poly_b.GetNumberOfPoints()
                    n_points_if = 1. / float(n_points)
                    for i in range(n_points):
                        if selector.IsInsideSurface(poly_b.GetPoint(i)) > 0:
                            count += 1
                            over = count * n_points_if
                            if over > over_tolerance:
                                return True

        return False


class Polymer(ABC):
    """
    Abstract class for modeling a Polymer (a sequence of monomers)
    """

    def __init__(self, m_surf):
        """
        Constructor
        :param m_surf: monomer surface (as vtkPolyData object)
        """
        h_diam = poly_max_distance(m_surf)
        assert h_diam > 0
        self.__m_surf = m_surf
        self.__m_diam = h_diam
        self.__m = list()
        self.__p = None
        self.__u = None
        self.__t, self.__r, self.__q = list(), list(), list()
        self.__t_length = 0

    def get_vol(self, fast=True):
        """
        Get the polymer volume
        :param fast: if True (default) the volume monomer is only computed once
        :return: the computed volume
        """
        if fast:
            vol = 0
            if len(self.__m) == 0:
                return vol
            else:
                vol = self.__m[0].get_vol()
                for m in self.__m[1:]:
                    vol += m.get_vol()
                return vol
        else:
            vol = 0
            for m in self.__m:
                vol += m.compute_vol()
            return vol

    def get_total_len(self):
        return self.__t_length

    def get_num_monomers(self):
        return len(self.__m)

    def get_tail_point(self):
        """
        Get the central coordinate for the latest monomer
        :return: point coordinates as ndarray
        """
        return self.__r[-1]

    def get_vtp(self):
        """
        Get the polymer a skeleton, each momomer is point and lines conecting monomers
        :return: a vtkPolyData
        """

        app_flt = vtk.vtkAppendPolyData()

        # Polymers loop
        for m in self.__m:
            app_flt.AddInputData(m.get_vtp())
        app_flt.Update()

        return app_flt.GetOutput()

    def get_skel(self):
        """
        Get the polymer as a skeleton, each momomer is point and lines conecting monomers
        :return: a vtkPolyData
        """

        # Initialization
        poly, points = vtk.vtkPolyData(), vtk.vtkPoints()
        verts, lines = vtk.vtkCellArray(), vtk.vtkCellArray()

        # Monomers loop
        for i in range(1, len(self.__r)):
            id_p0, id_p1 = points.InsertNextPoint(self.__r[i - 1]), points.InsertNextPoint(self.__r[i])
            verts.InsertNextCell(1)
            verts.InsertCellPoint(id_p0)
            lines.InsertNextCell(2)
            lines.InsertCellPoint(id_p0)
            lines.InsertCellPoint(id_p1)

        # Construct poly
        poly.SetPoints(points)
        poly.SetVerts(verts)
        poly.SetLines(lines)

        return poly

    def add_monomer(self, r, t, q, m):
        """
        Add a new monomer surface to the polymer once affine transformation is known
        :param r: center point
        :param t: tangent vector
        :param q: unit quaternion for rotation
        :param m: monomer
        :return:
        """
        assert isinstance(r, np.ndarray) and isinstance(t, np.ndarray) and isinstance(q, np.ndarray) \
               and isinstance(m, Monomer)
        assert (len(r) == 3) and (len(t) == 3) and (len(q) == 4)
        self.__r.append(r)
        self.__t.append(t)
        self.__q.append(q)
        self.__m.append(m)
        # Update total length
        if self.get_num_monomers() <= 1:
            self.__t_length = 0
        else:
            self.__t_length += points_distance(self.__r[-1], self.__r[-2])

    def insert_density_svol(self, m_svol, tomo, v_size=1, merge='min'):
        """
        Insert a polymer as set of subvolumes into a tomogram
        :param m_svol: input monomer sub-volume reference
        :param tomo: tomogram where m_svol is added
        :param v_size: tomogram voxel size (default 1)
        :param merge: merging mode, valid: 'max' (default), 'min', 'sum' and 'insert'
        :return:
        """
        for mmer in self.__m:
            mmer.insert_density_svol(m_svol, tomo, v_size, merge=merge)

    @abstractmethod
    def set_reference(self):
        raise NotImplementedError

    @abstractmethod
    def gen_new_monomer(self):
        raise NotImplementedError

    def overlap_polymer(self, monomer, over_tolerance=0):
        """
        Determines if a monomer overlaps with other polymer's monomers
        :param monomer: input monomer
        :param over_tolerance: fraction of overlapping tolerance (default 0)
        :return: True if there is an overlapping and False otherwise
        """

        # Initialization
        selector = vtk.vtkSelectEnclosedPoints()
        selector.SetTolerance(VTK_RAY_TOLERANCE)
        selector.Initialize(monomer.get_vtp())

        # Polymer loop, no need to process monomer beyond diameter distance
        diam, center = monomer.get_diameter(), monomer.get_center_mass()
        for i in range(len(self.__m) - 1, 0, -1):
            hold_monomer = self.__m[i]
            dist = points_distance(center, hold_monomer.get_center_mass())
            if dist <= diam:
                poly_b = hold_monomer.get_vtp()
                count, n_points = 0., poly_b.GetNumberOfPoints()
                n_points_if = 1. / float(n_points)
                for i in range(n_points):
                    if selector.IsInsideSurface(poly_b.GetPoint(i)) > 0:
                        count += 1
                        over = count * n_points_if
                        if over > over_tolerance:
                            return True

        return False


class SAWLC(Polymer):
    """
    Class for fibers following model Self-Avoiding Worm-Like Chain (SAWLC)
    """

    def __init__(self, l_length, m_surf, p0=(0, 0, 0)):
        """
        Constructor
        :param l_lengh: link length
        :param m_surf: monomer surface (as vtkPolyData object)
        :param p0: starting point
        """
        super(SAWLC, self).__init__(m_surf)
        assert l_length > 0
        self.__l = l_length
        self.set_reference(p0)

    def set_reference(self, p0=(0., 0., 0)):
        """
        Initializes the chain with the specified point input point, if points were introduced before the are forgotten
        :param p0: starting point
        :return:
        """
        assert hasattr(p0, '__len__') and (len(p0) == 3)
        self.__p = np.asarray(p0)
        hold_monomer = Monomer(self._Polymer__m_surf, self._Polymer__m_diam)
        hold_q = gen_rand_unit_quaternion()
        # hold_q = np.asarray((1, 0., 0., 1.), dtype=np.float32)
        hold_monomer.rotate_q(hold_q)
        hold_monomer.translate(p0)
        self.add_monomer(p0, np.asarray((0., 0., 0.)), hold_q, hold_monomer)

    def gen_new_monomer(self, over_tolerance=0, voi=None, v_size=1):
        """
        Generates a new monomer for the polymer according to the specified random model
        TODO: According to current implementation only tangential and axial rotation angles are chosen from a uniform
              random distribution
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0)
        :param voi: VOI to define forbidden regions (default None, not applied)
        :param v_size: VOI voxel size, it must be greater than 0 (default 1)
        :return: a 4-tuple with monomer center point, associated tangent vector, rotated quaternion and monomer,
                 return None in case the generation has failed
        """

        # Translation
        t = gen_uni_s2_sample(np.asarray((0., 0., 0.)), self.__l)
        r = self._Polymer__r[-1] + t

        # Rotation
        q = gen_rand_unit_quaternion()
        # q = np.asarray((1, 0, 0, 1), dtype=np.float32)

        # Monomer
        # hold_m = self._Polymer__m[-1].get_copy()
        hold_m = Monomer(self._Polymer__m_surf, self._Polymer__m_diam)
        hold_m.rotate_q(q)
        hold_m.translate(r)

        # Check self-avoiding and forbidden regions
        if self.overlap_polymer(hold_m, over_tolerance=over_tolerance):
            return None
        elif voi is not None:
            if hold_m.overlap_voi(voi, v_size):
                return None

        return r, t, q, hold_m


class HelixFiber(Polymer):
    """
    Class for modelling a random helical flexible fiber
    """

    def __init__(self, l_length, m_surf, p_length, hp_length, mz_length, z_length_f=0, p0=(0, 0, 0), vz=(0, 0, 1)):
        """
        Constructor
        :param l_length: link length
        :param m_surf: monomer surface (as vtkPolyData object)
        :param p_length: persistence length
        :param hp_length: helix period length (distance required by azimuthal angle to cover 360deg)
        :param mz_length: monomer z-length
        :param z_length_f: z-length or helix elevation
        :param p0: starting point (default origin (0,0,0))
        :param vz: reference vector for z-axis (default (0, 0, 1)
        """
        super(HelixFiber, self).__init__(m_surf)
        assert (l_length > 0) and (p_length > 0) and (z_length_f >= 0) and (hp_length > 0) and (mz_length > 0)
        self.__l, self.__lp, self.__lz = l_length, p_length, l_length * z_length_f
        self.__hp, self.__mz_length = hp_length, mz_length
        self.__hp_astep = (360. * self.__mz_length) / self.__hp
        self.__compute_helical_parameters()
        assert hasattr(vz, '__len__') and (len(vz) == 3)
        self.__vz = np.asarray(vz, dtype=float)
        # Curve state member variables
        self.__ct, self.__za, self.__rq = 0., 0., np.asarray((1., 0., 0., 0.)) # z-aligned curve time (considering speed 1)
        self.set_reference(p0, vz)

    def set_reference(self, p0=(0., 0., 0.), vz=(0., 0., 1.)):
        """
        Initializes the chain with the specified point input point, if points were introduced before the are forgotten
        :param p0: starting point
        :param vz: z-axis reference vector for helicoidal parametrization
        :return:
        """
        assert hasattr(p0, '__len__') and (len(p0) == 3)
        self.__p = np.asarray(p0)
        self.__rq = gen_rand_unit_quaternion()
        hold_monomer = Monomer(self._Polymer__m_surf, self._Polymer__m_diam)
        vzr = gen_uni_s2_sample(np.asarray((0., 0., 0.)), 1.)
        M = vect_to_zmat(vzr, mode='passive')
        hold_q = rot_to_quat(M)
        vzr /= vector_module(vzr)
        vzr *= self.__mz_length
        hold_monomer.rotate_q(hold_q)
        hold_monomer.translate(p0)
        self.add_monomer(p0, vzr, hold_q, hold_monomer)

    def gen_new_monomer(self, over_tolerance=0, voi=None, v_size=1, net=None, branch=None):
        """
        Generates a new monomer according the flexible fiber model
        :param over_tolerance: fraction of overlapping tolerance for self avoiding (default 0)
        :param voi: VOI to define forbidden regions (default None, not applied)
        :param v_size: VOI voxel size, it must be greater than 0 (default 1)
        :param net: if not None (default) it contain a network of polymer that must be avoided
        :param branch: input branch from where the current mmer starts, is avoid network avoiding at the branch,
                       only valid in net is not None (default None).
        :return: a 4-tuple with monomer center point, associated tangent vector, rotated quaternion and monomer,
                 return None in case the generation has failed
        """

        hold_m = Monomer(self._Polymer__m_surf, self._Polymer__m_diam)

        # Rotation
        t = self._Polymer__t[-1] + self.__compute_tangent(self.__ct)
        t = t * (self.__mz_length / vector_module(t))
        self.__za = wrap_angle(self.__za + self.__hp_astep)
        q1 = angle_axis_to_quat(self.__za, t[0], t[1], t[2])
        M = vect_to_zmat(t, mode='passive')
        q = rot_to_quat(M)
        hold_m.rotate_q(quat_mult(q, q1))

        # Translation
        hold_r = self._Polymer__r[-1]
        self.__ct += self.__l
        r = hold_r + t
        hold_m.translate(r)

        # Avoid forbidden regions
        if voi is not None:
            if hold_m.overlap_voi(voi, v_size):
                return None
        # Self-avoiding and network avoiding
        if branch is None:
            if self.overlap_polymer(hold_m, over_tolerance=over_tolerance):
                return None
            if net is not None:
                if hold_m.overlap_net(net, over_tolerance=over_tolerance):
                    return None
        else:
            branch_dst = points_distance(branch.get_point(), hold_m.get_center_mass())
            if branch_dst > hold_m.get_diameter():
                if self.overlap_polymer(hold_m, over_tolerance=over_tolerance):
                    return None
                if net is not None:
                    if hold_m.overlap_net(net, over_tolerance=over_tolerance):
                        return None

        return r, t, q, hold_m

    def __compute_helical_parameters(self):
        """
        Private method (fill class member variables) to compute helical parameters (a and b) from persistence length
        :return:
        """

        # Compute curvature from persistence
        k = math.acos(math.exp(-self.__l / self.__lp)) / self.__l

        # Compute helix b parameter from z-length and persistence
        self.__b = self.__lz / self.__lp
        assert self.__b >= 0

        # Compute helix a parameter from k and b
        self.__a = (1 + math.sqrt(1-4.*k*self.__b*self.__b)) / (2. * k)
        assert self.__a > 0

    def __compute_tangent(self, t):
        """
        Computes curve (z-aligned axis) normalized tangent vector
        :param t: input parameter, time assuming that speed is 1
        :return: returns the normalized tangent vector (3 elements array)
        """
        sq = math.sqrt(self.__a * self.__a + self.__b * self.__b)
        # s = sq * t
        s = t
        t = (1. / sq) * np.asarray((self.__b, -self.__a * math.sin(s / sq), self.__a * math.cos(s / sq)))
        return rot_vect_quat(t, self.__rq)


class FiberUnit(ABC):
    """
    Abstract class to generate fiber unit (set of monomers)
    """

    @abstractmethod
    def get_vtp(self):
        raise NotImplementedError

    @abstractmethod
    def get_tomo(self):
        raise NotImplementedError


class FiberUnitSDimer(FiberUnit):
    """
    Class for modeling a fiber unit as dimer of two spheres
    """

    def __init__(self, sph_rad, v_size=1):
        """
        Constructor
        :param sph_rad: radius for spheres
        :param v_size: voxel size (default 1)
        """
        assert (sph_rad > 0) and (v_size > 0)
        self.__sph_rad, self.__v_size = float(sph_rad), float(v_size)
        self.__size = int(math.ceil(6. * (sph_rad / v_size)))
        if self.__size%2 != 0:
            self.__size += 1
        self.__tomo, self.__surf = None, None
        self.__gen_sdimer()

    def get_vtp(self):
        return self.__surf

    def get_tomo(self):
        return self.__tomo

    def __gen_sdimer(self):
        """
        Contains the procedure to generate the Dimer of spheres with the specified size by using logistic functions
        """

        # Input parsing
        sph_rad_v = self.__sph_rad / self.__v_size

        # Generating the grid
        self.__tomo = np.zeros(shape=(self.__size, self.__size, self.__size), dtype=np.float32)
        dx, dy, dz = float(self.__tomo.shape[0]), float(self.__tomo.shape[1]), float(self.__tomo.shape[2])
        dx2, dy2, dz2 = math.floor(.5 * dx), math.floor(.5 * dy), math.floor(.5 * dz)
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(np.arange(x_l, x_h), np.arange(y_l, y_h), np.arange(z_l, z_h), indexing='xy')
        X += .5
        Y += .5
        Z += .5
        # X, Y, Z = X.astype(np.float16), Y.astype(np.float16), X.astype(np.float16)

        # Generate the first unit
        Yh = Y + sph_rad_v
        R = np.abs(X * X + Yh * Yh + Z * Z)
        self.__tomo += 1. / (1. + np.exp(-R))

        # Generate the second unit
        Yh = Y - sph_rad_v
        R = np.abs(X * X + Yh * Yh + Z * Z)
        self.__tomo += 1. / (1. + np.exp(-R))

        # Generating the surfaces
        self.__tomo = lin_map(self.__tomo, lb=0, ub=1)
        self.__surf = iso_surface(self.__tomo, .75)
        self.__surf = poly_scale(self.__surf, self.__v_size)
        self.__surf = poly_translate(self.__surf, -.5 * self.__v_size * (np.asarray(self.__tomo.shape)-.5))


class MTUnit(FiberUnit):
    """
    Class for modelling a fiber unit for microtubules (MTs)
    """

    def __init__(self, sph_rad=40, mt_rad=100.5, n_units=13, v_size=1):
        """
        Constructor
        :param sph_rad: radius for spheres (default 40, approximate tubulin radius in A)
        :param mt_rad: microtubule radius (default 100.5, approximate microtubule radius in A)
        :param n_units: number of units (default 13, number of protofilaments that compund a MT)
        :param v_size: voxel size (default 1)
        """
        assert (sph_rad > 0) and (mt_rad > 0) and (n_units > 0) and (v_size > 0)
        self.__sph_rad, self.__mt_rad, self.__n_units, self.__v_size = float(sph_rad), float(mt_rad), int(n_units),\
                                                                       float(v_size)
        self.__size = int(math.ceil(6. * (sph_rad / v_size)))
        if self.__size % 2 != 0:
            self.__size += 1
        self.__tomo, self.__surf = None, None
        self.__gen_sdimer()

    def get_vtp(self):
        return self.__surf

    def get_tomo(self):
        return self.__tomo

    def __gen_sdimer(self):
        """
        Contains the procedure to generate the Dimer of spheres with the specified size by using logistic functions
        """

        # Input parsing
        sph_rad_v, mt_rad_v = self.__sph_rad / self.__v_size, self.__mt_rad / self.__v_size

        # Generating the grid
        self.__tomo = np.zeros(shape=(self.__size, self.__size, self.__size), dtype=np.float32)
        dx, dy, dz = float(self.__tomo.shape[0]), float(self.__tomo.shape[1]), float(self.__tomo.shape[2])
        dx2, dy2, dz2 = math.floor(.5 * dx), math.floor(.5 * dy), math.floor(.5 * dz)
        x_l, y_l, z_l = -dx2, -dy2, -dz2
        x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
        X, Y, Z = np.meshgrid(np.arange(x_l, x_h), np.arange(y_l, y_h), np.arange(z_l, z_h), indexing='xy')
        X += .5
        Y += .5
        Z += .5
        # X, Y, Z = X.astype(np.float16), Y.astype(np.float16), X.astype(np.float16)

        # Loop for generate the units
        Z2 = Z * Z
        ang_step = 2. * np.pi / self.__n_units
        ang = ang_step
        while ang <= 2. * np.pi:
            # Generate the unit
            x, y = mt_rad_v * math.cos(ang), mt_rad_v * math.sin(ang)
            Xh, Yh = X + x, Y + y
            R = np.abs(Xh * Xh + Yh * Yh + Z2)
            self.__tomo += 1. / (1. + np.exp(-R))
            ang += ang_step

        # Generating the surfaces
        self.__tomo = lin_map(self.__tomo, lb=0, ub=1)
        self.__surf = iso_surface(self.__tomo, .75)
        self.__surf = poly_scale(self.__surf, self.__v_size)
        self.__surf = poly_translate(self.__surf, -.5 * self.__v_size * (np.asarray(self.__tomo.shape) - .5))



