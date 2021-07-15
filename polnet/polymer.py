"""
Classes for modeling Self-Avoiding Worm-Like Chain (SAWLC) polymers.
A polymer is a linear sequence of monomers
"""

__author__ = 'Antonio Martinez-Sanchez'

import copy
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
        self.__trans_q = dict()
        self.__trans_q['t'] = list()
        self.__trans_q['r'] = list()

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
        self.__trans_q['r'].append(q)

    def translate(self, t_v):
        """
        Applies rotation rigid transformation.
        :param t_v: translation vector (x, y, z)
        :return:
        """
        self.__m_surf = poly_translate(self.__m_surf, t_v)
        self.compute_bounds()
        self.__trans_q['t'].append(t_v)

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

    def insert_density_svol(self, m_svol, tomo, v_size=1, merge='max'):
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
        for i in range(len(self.__trans_q['t'])):
            tot_v += (self.__trans_q['t'][i] * v_size_i)
            hold_svol = tomo_rotate(m_svol, self.__trans_q['r'][i])
        insert_svol_tomo(hold_svol, tomo, tot_v, merge=merge)


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

    def insert_density_svol(self, m_svol, tomo, v_size=1, merge='max'):
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


class SAWLC(Polymer):
    """
    Class for fibers following the model Self-Avoiding Worm-Like Chain (SAWLC)
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

    def set_reference(self, p0):
        """
        Initializes the chain with the specified point input point, if points were introduced before the are forgotten
        :param p0: starting point
        :return:
        """
        assert hasattr(p0, '__len__') and (len(p0) == 3)
        self.__p = np.asarray(p0)
        hold_monomer = Monomer(self.__m_surf, self.__m_diam)
        hold_q = gen_rand_unit_quaternion()
        # hold_q = np.asarray((1, 0., 0., 1.), dtype=np.float32)
        hold_monomer.rotate_q(hold_q)
        hold_monomer.translate(p0)
        self.add_monomer(p0, np.asarray((0., 0., 0.)), hold_q, hold_monomer)

    def gen_new_monomer(self, over_tolerance=0, voi=None, v_size=1):
        """
        Generates a new monomer for the polymer according to specified model
        TODO: According to current implementation only tangential and axial rotation angles are chosen from a unifrom
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

    def overlap_polymer(self, monomer, over_tolerance=0):
        """
        Determines if a monomer overlaps with polymer
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
        for i in range(len(self._Polymer__m) - 1, 0, -1):
            hold_monomer = self._Polymer__m[i]
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


class HelixFiber(Polymer):
    """
    Class for modelling a helical flexible fiber
    """

    def __init__(self, l_length, m_surf, p_length, z_length_f=0, z_ang=0, p0=(0, 0, 0), vz=(0, 0, 1)):
        """
        Constructor
        :param l_length: link length
        :param m_surf: monomer surface (as vtkPolyData object)
        :param p_length: persistence length
        :param z_length_f: z-length or helix elevation
        :param z_ang: azimutal (around reference z-axis) angle (default 0)
        :param p0: starting point (default origin (0,0,0))
        :param vz: reference vector for z-axis (default (0, 0, 1)
        """
        super(HelixFiber, self).__init__(m_surf)
        assert (l_length > 0) and (p_length > 0) and (z_length_f >= 0)
        self.__l, self.__lp, self.__lz = l_length, p_length, l_length * z_length_f
        assert z_ang >=0
        self.__compute_helical_parameters()
        self.set_reference(p0, vz)

    def gen_new_monomer(self, over_tolerance=0, voi=None, v_size=1):
        """
        Generates a new monomer according the flexible fiber model
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

