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

GTRUTH_VTP_LBLS = 'gt_labels'


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


def poly_surface_area(poly):
    """
    Computes the surface area of polydata
    :param poly: input vtkPolyData
    :return: the volume computed
    """
    assert isinstance(poly, vtk.vtkPolyData)
    mass = vtk.vtkMassProperties()
    mass.SetInputData(poly)
    return mass.GetSurfaceArea()


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


def merge_polys(poly_1, poly_2):
    """
    Merges two input poly_data in single one
    :param poly_1: input poly_data 1
    :param poly_2: input poly_data 2
    :return: an poly_data that merges the two inputs
    """
    assert isinstance(poly_1, vtk.vtkPolyData) and isinstance(poly_2, vtk.vtkPolyData)
    app_flt = vtk.vtkAppendPolyData()
    app_flt.AddInputData(poly_1)
    app_flt.AddInputData(poly_2)
    app_flt.Update()
    return app_flt.GetOutput()


def add_label_to_poly(poly, lbl, p_name, mode='cell'):
    """
    Add a label to all cells in a poly_data
    :param poly: input poly_data
    :param lbl: label (integer) value
    :param p_name: property name used for labels, if not exist in poly_dota is created
    :param mode: selected wheter the label is added to cells ('cell'), points ('point') or both ('both')
    """
    assert (mode == 'cell') or (mode == 'point') or (mode == 'both')
    assert isinstance(poly, vtk.vtkPolyData)
    lbl, p_name = int(lbl), str(p_name)

    if mode == 'cell' or mode == 'both':
        arr = vtk.vtkIntArray()
        n_cells = poly.GetNumberOfCells()
        arr.SetName(p_name)
        arr.SetNumberOfComponents(1)
        arr.SetNumberOfValues(n_cells)
        for i in range(n_cells):
            arr.SetValue(i, lbl)
        poly.GetCellData().AddArray(arr)
    elif mode == 'cell' or mode == 'both':
        arr = vtk.vtkIntArray()
        n_points = poly.GetNumberOfPoints()
        arr.SetName(p_name)
        arr.SetNumberOfComponents(1)
        arr.SetNumberOfValues(n_points)
        for i in range(n_points):
            arr.SetValue(i, lbl)
        poly.GetPointData().AddArray(arr)


def points_to_poly_spheres(points, rad):
    """
    From an array of coordinates generates a poly_data associating a sphere centered a each point
    :param points: array or list n points with shape [n, 3]
    :param rad: sphere radius
    :return: an output poly_data
    """
    assert hasattr(points, '__len__') and (len(points) > 0) and (len(points[0]) == 3)
    rad = float(rad)
    app_flt = vtk.vtkAppendPolyData()

    for i in range(len(points)):
        center = points[i]
        vtp_source = vtk.vtkSphereSource()
        vtp_source.SetCenter(center[0], center[1], center[2])
        vtp_source.SetRadius(rad)
        vtp_source.Update()
        app_flt.AddInputData(vtp_source.GetOutput())

    app_flt.Update()
    return app_flt.GetOutput()


def poly_max_distance(vtp):
    """
    Computes the maximum distance in vtkPolyData
    :param vtp: input vtkPolyData
    :return: the maximum distance as real value
    """
    if vtp.GetNumberOfPoints() <= 1:
        return 0
    else:
        mx = 0
        for i in range(0, vtp.GetNumberOfPoints() - 1):
            ref_p = np.asarray(vtp.GetPoint(i))
            for j in range(i + 1, vtp.GetNumberOfPoints()):
                hold_p = np.asarray(vtp.GetPoint(j))
                hold_mx = points_distance(ref_p, hold_p)
                if hold_mx > mx:
                    mx = hold_mx
        return mx

def poly_diam(vtp):
    """
    Computes the diameter of a polydata, approximated to two times the maximumd point distance to its center of mass
    :param vtp: input vtkPolyData
    :return: the maximum distance as real value
    """
    if vtp.GetNumberOfPoints() <= 1:
        return 0
    else:
        mx = 0
        ref_p = poly_center_mass(vtp)
        for i in range(0, vtp.GetNumberOfPoints()):
            hold_p = np.asarray(vtp.GetPoint(i))
            hold_mx = points_distance(ref_p, hold_p)
            if hold_mx > mx:
                mx = hold_mx
        return mx


def poly_point_min_dst(poly, point, chull=False):
    """
    Compute the minimum distance from a point to a poly
    :param poly: input poly
    :param point: input point
    :param chull: computation mode, if True (default False) the convex hull surface is firstly extracted to avoid poly holes,
                'otherwise the minimum distance is directly computed
    :return: the minimum distance found
    """

    if poly.GetNumberOfPoints() <= 0:
        return 0
    else:
        mn = np.finfo(float).max
        if chull:
            poly = convex_hull_surface(poly)
        ref_p = np.asarray(point, dtype=float)
        for j in range(0, poly.GetNumberOfPoints()):
            hold_p = np.asarray(poly.GetPoint(j))
            hold_mn = points_distance(ref_p, hold_p)
            if hold_mn < mn:
                mn = hold_mn
        return mn


def poly_center_mass(poly):
    """
    Computes the center of mass of polydata
    :param poly: input poly
    :return: center of mass coordinates
    """
    cm_flt = vtk.vtkCenterOfMass()
    cm_flt.SetInputData(poly)
    cm_flt.Update()
    return np.asarray(cm_flt.GetCenter())


def convex_hull_surface(poly):
    """
    Extract the convex full surface of a polydata
    :param poly: input polydata
    :return: convex hull surface
    """
    convexHull = vtk.vtkDelaunay3D()
    convexHull.SetInputData(poly)
    outerSurface = vtk.vtkGeometryFilter()
    convexHull.Update()
    outerSurface.SetInputData(convexHull.GetOutput())
    outerSurface.Update()

    return outerSurface.GetOutput()


def poly_decimate(poly, dec):
    """
    Decimate a vtkPolyData
    :param poly: input vtkPolyData
    :param dec: Specify the desired reduction in the total number of polygons, default None (not applied)
               (e.g., if TargetReduction is set to 0.9,
               this filter will try to reduce the data set to 10% of its original size).
    :return: the input poly filtered
    """
    tr_dec = vtk.vtkDecimatePro()
    tr_dec.SetInputData(poly)
    tr_dec.SetTargetReduction(dec)
    tr_dec.Update()
    return tr_dec.GetOutput()