"""
Functionality of processing PolyData
"""

__author__ = 'Antonio Martinez-Sanchez'

import math

from polnet.affine import *
from polnet.utils import *
import vtk
import numpy as np
from scipy import stats

# CONSTANTS


# FUNCTIONS
from polnet.utils import wrap_angle


def find_point_on_poly(point, poly):
    """
    Find the closest point on a poly to a reference point
    :param point: input reference point
    :param poly: poly data where the closest output point has to be found
    :return: output point
    """
    assert hasattr(point, '__len__') and (len(point) == 3)

    point_tree = vtk.vtkKdTreePointLocator()
    point_tree.SetDataSet(poly)
    point_tree.BuildLocator()
    cpoint_id = point_tree.FindClosestPoint(point)
    normals = poly.GetPointData().GetNormals()
    return np.asarray(poly.GetPoint(cpoint_id)), np.asarray(normals.GetTuple(cpoint_id))


def gen_rand_quaternion_on_vector(vect):
    """
    Generates a unit quaternion which represents a random rotation around an input reference vector from Z-axis
    :param vect: reference vector
    :return: a quaternion which represents the rotation from Z-axis unit vector to be aligned to reference vector, plus
             a random rotation around the axis defined by the reference vector.
    """
    assert hasattr(vect, '__len__') and (len(vect) == 3)

    rnd_ang = 360. * np.random.random() - 180.
    q1 = angle_axis_to_quat(rnd_ang, vect[0], vect[1], vect[2])
    M = vect_to_zmat(np.asarray(vect), mode='passive')
    q = rot_to_quat(M)
    return quat_mult(q, q1)


def gen_uni_s2_sample_on_poly(center, rad, thick, poly):
    """
    Generates a coordinate from an approximately uniformly random distribution on the intersection between a holow
    sphere an a PolyData
    :param center: sphere center
    :param rad: sphere radius
    :param poly: input poly (vtkPolyData object)
    :param thick: hollow sphere thickness
    :return: the random coordinate generated or None if no intersection
    """
    assert hasattr(center, '__len__') and (len(center) == 3)
    assert (rad > 0) and (thick > 0)
    assert isinstance(poly, vtk.vtkPolyData)

    # Find poly points within rad+thick
    kdtree = vtk.vtkKdTreePointLocator()
    kdtree.SetDataSet(poly)
    kdtree.BuildLocator()
    pids = vtk.vtkIdList()
    kdtree.FindPointsWithinRadius(rad + .5*thick, center, pids)

    # save_vtp(vtp_inter, './out/hold_3.vtp')

    # Get a points randomly un util a point is found in the intersection
    min_dst, n_pts = rad - 0.5*thick, pids.GetNumberOfIds()
    for i in np.random.randint(0, n_pts, n_pts):
        pt = poly.GetPoint(pids.GetId(i))
        hold = center - pt
        if math.sqrt((hold * hold).sum()) > min_dst:
            return pt
    return None


def gen_uni_s2_sample_on_poly_inter(center, rad, poly, sph_res=360):
    """
    Generates a coordinate from an approximately uniformly random distribution on the intersection between an sphere
    an a PolyData
    @Deprecated: the usage vtkIntersectionPolyDataFilter makes this function too slow
    :param center: sphere center
    :param rad: sphere radius
    :param poly: input poly (vtkPolyData object)
    :param sph_res: resolution for generating the surfaces polydata for both (longitude and latitude resolution),
                    default is 360 (1 degree resolution)
    :return: the random coordinate generated or None if no intersection
    """
    assert hasattr(center, '__len__') and (len(center) == 3)
    assert rad > 0
    assert isinstance(poly, vtk.vtkPolyData)

    # Generate the sphere polydata
    vtp_source = vtk.vtkSphereSource()
    vtp_source.SetCenter(center[0], center[1], center[2])
    vtp_source.SetRadius(rad)
    vtp_source.SetPhiResolution(36)
    vtp_source.SetThetaResolution(36)
    vtp_source.Update()
    vtp_sphere = vtp_source.GetOutput()

    # # Debug
    # from polnet.lio import save_vtp
    # save_vtp(vtp_sphere, './out/hold_1.vtp')
    # save_vtp(poly, './out/hold_2.vtp')

    # Compute polydata objects intersection
    inter_flt = vtk.vtkIntersectionPolyDataFilter()
    inter_flt.SetInputDataObject(0, vtp_sphere)
    inter_flt.SetInputDataObject(1, poly)
    inter_flt.Update()
    vtp_inter = inter_flt.GetOutput()

    # save_vtp(vtp_inter, './out/hold_3.vtp')

    # Get a point randomly on intersection
    n_pts = vtp_inter.GetNumberOfPoints()
    if n_pts == 0:
        return None
    else:
        rnd_id = np.random.randint(0, vtp_inter.GetNumberOfPoints()-1, 1)
        return vtp_inter.GetPoint(rnd_id)


def poly_reverse_normals(poly):
    """
    Reverse the normals of an input polydata
    :param poly: input vtkPolyData object
    :return: a vtkPolyData object copy of the input but with the normals reversed
    """
    assert isinstance(poly, vtk.vtkPolyData)
    reverse = vtk.vtkReverseSense()
    reverse.SetInputData(poly)
    reverse.ReverseNormalsOn()
    reverse.Update()
    return reverse.GetOutput()


def poly_volume(poly):
    """
    Computes the volume of polydata
    :param poly: input vtkPolyData
    :return: the volume computed
    """
    assert isinstance(poly, vtk.vtkPolyData)
    mass = vtk.vtkMassProperties()
    mass.SetInputData(poly)
    return mass.GetVolume()


def add_sfield_to_poly(poly, sfield, name, dtype='float', interp='NN', mode='points'):
    """
    Add the values of a scalar field to a vtkPolyData object as point property
    :param poly: vtkPolyData objects where the scalar field values will be added
    :param sfield: input scalar field as ndarray
    :param name: string with name associated to the added property
    :param dtype: data type, valid 'float' or 'int'
    :param interp: interpolation mode, valid 'NN'-> nearest neighbour and 'trilin'-> trilinear
    :param mode: determines if the scalar field is either added to vtkPolyData points ('points', defualt) or
                 cells ('cells')
    """
    assert isinstance(sfield, np.ndarray)
    assert isinstance(name, str)
    assert (dtype == 'float') or (dtype == 'int')
    assert (interp == 'NN') or (interp == 'trilin')
    if interp == 'trilin':
        interp_func = trilin_interp
    else:
        interp_func = nn_iterp
    assert (mode == 'points') or (mode == 'cells')

    if mode == 'points':
        # Creating and adding the new property as a new array for PointData
        n_points = poly.GetNumberOfPoints()
        if dtype == 'int':
            arr = vtk.vtkIntArray()
        else:
            arr = vtk.vtkFloatArray()
        arr.SetName(name)
        arr.SetNumberOfComponents(1)
        arr.SetNumberOfValues(n_points)
        for i in range(n_points):
            x, y, z = poly.GetPoint(i)
            arr.SetValue(i, interp_func(x, y, z, sfield))
        poly.GetPointData().AddArray(arr)
    else:
        # Creating and adding the new property as a new array for CellData
        if dtype == 'int':
            arr = vtk.vtkIntArray()
        else:
            arr = vtk.vtkFloatArray()
        n_cells = poly.GetNumberOfCells()
        arr.SetName(name)
        arr.SetNumberOfComponents(1)
        arr.SetNumberOfValues(n_cells)
        for i in range(n_cells):
            cell = vtk.vtkGenericCell()
            poly.GetCell(i, cell)
            pts = cell.GetPoints()
            n_pts = pts.GetNumberOfPoints()
            if dtype == 'int':
                values = np.zeros(shape=n_pts, dtype=int)
            else:
                values = np.zeros(shape=n_pts, dtype=float)
            for j in range(n_pts):
                x, y, z = pts.GetPoint(j)
                values[j] = interp_func(x, y, z, sfield)
            arr.SetValue(i, stats.mode(values)[0][0])
        poly.GetCellData().AddArray(arr)