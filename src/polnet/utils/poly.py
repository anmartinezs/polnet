"""VTK PolyData processing utilities.

Contains functions for isosurface extraction, scalar field
assignment, thresholding, and geometric queries on vtkPolyData
objects used to represent membrane surfaces and protein models.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import math

import numpy as np
import vtk
from scipy import stats
from vtkmodules.util import numpy_support

from ..logging_conf import _LOGGER as logger
from .affine import (
    angle_axis_to_quat,
    quat_mult,
    rot_to_quat,
    vect_to_zmat,
)
from .utils import (
    nn_iterp,
    points_distance,
    trilin_interp,
)

GTRUTH_VTP_LBLS = "gt_labels"


def find_point_on_poly(point, poly):
    """Find the closest point on a vtkPolyData to a reference point.

    Args:
        point (array-like): Reference 3-D coordinate.
        poly (vtk.vtkPolyData): Surface to search.

    Returns:
        tuple[numpy.ndarray, numpy.ndarray]: A pair
            (closest_point, normal) both as 3-element arrays.

    Raises:
        ValueError: If point does not have exactly 3 elements.
    """
    if not hasattr(point, "__len__") or len(point) != 3:
        raise ValueError("point must have exactly 3 elements.")

    point_tree = vtk.vtkKdTreePointLocator()
    point_tree.SetDataSet(poly)
    point_tree.BuildLocator()
    cpoint_id = point_tree.FindClosestPoint(point)
    normals = poly.GetPointData().GetNormals()
    return np.asarray(poly.GetPoint(cpoint_id)), np.asarray(
        normals.GetTuple(cpoint_id)
    )


def gen_rand_quaternion_on_vector(vect):
    """Generate a unit quaternion aligning Z-axis to a vector.

    Combines a rotation from the Z-axis to vect with a uniform
    random spin around vect, producing an orientation compatible
    with placing a molecule on a membrane normal.

    Args:
        vect (array-like): Reference direction as a 3-element
            vector.

    Returns:
        numpy.ndarray: A unit quaternion [w, x, y, z].

    Raises:
        ValueError: If vect does not have exactly 3 elements.
    """
    if not hasattr(vect, "__len__") or len(vect) != 3:
        raise ValueError("vect must have exactly 3 elements.")

    rnd_ang = 360.0 * np.random.random() - 180.0
    q1 = angle_axis_to_quat(rnd_ang, vect[0], vect[1], vect[2])
    M = vect_to_zmat(np.asarray(vect), mode="passive")
    q = rot_to_quat(M)
    return quat_mult(q, q1)


def gen_uni_s2_sample_on_poly(center, rad, thick, poly):
    """Sample a random point on a hollow-sphere / PolyData intersection.

    Finds all PolyData points within a hollow sphere shell
    [rad - thick/2, rad + thick/2] centred at center and returns
    one uniformly at random.

    Args:
        center (array-like): Sphere centre as a 3-element vector.
        rad (float): Sphere radius in voxels.
        thick (float): Shell thickness in voxels.
        poly (vtk.vtkPolyData): Surface to sample from.

    Returns:
        tuple | None: A 3-D coordinate on the surface, or None if
            no intersection is found.

    Raises:
        ValueError: If center has != 3 elements, or rad/thick <= 0.
        TypeError: If poly is not a vtkPolyData instance.
    """
    if not hasattr(center, "__len__") or len(center) != 3:
        raise ValueError("center must have exactly 3 elements.")
    if rad <= 0 or thick <= 0:
        raise ValueError("rad and thick must be > 0.")
    if not isinstance(poly, vtk.vtkPolyData):
        raise TypeError("poly must be a vtkPolyData instance.")

    # Find poly points within rad+thick
    kdtree = vtk.vtkKdTreePointLocator()
    kdtree.SetDataSet(poly)
    kdtree.BuildLocator()
    pids = vtk.vtkIdList()
    kdtree.FindPointsWithinRadius(rad + 0.5 * thick, center, pids)

    min_dst, n_pts = rad - 0.5 * thick, pids.GetNumberOfIds()
    for i in np.random.randint(0, n_pts, n_pts):
        pt = poly.GetPoint(pids.GetId(i))
        hold = center - pt
        if math.sqrt((hold * hold).sum()) > min_dst:
            return pt
    return None


def poly_reverse_normals(poly):
    """Reverse the outward normals of a vtkPolyData surface.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        vtk.vtkPolyData: Copy of the input with normals flipped.

    Raises:
        TypeError: If poly is not a vtkPolyData instance.
    """
    if not isinstance(poly, vtk.vtkPolyData):
        raise TypeError("poly must be a vtkPolyData instance.")
    reverse = vtk.vtkReverseSense()
    reverse.SetInputData(poly)
    reverse.ReverseNormalsOn()
    reverse.Update()
    return reverse.GetOutput()


def poly_volume(poly):
    """Compute the enclosed volume of a closed vtkPolyData surface.

    Args:
        poly (vtk.vtkPolyData): Closed polygon dataset.

    Returns:
        float: Enclosed volume in voxel units cubed.

    Raises:
        TypeError: If poly is not a vtkPolyData instance.
    """
    if not isinstance(poly, vtk.vtkPolyData):
        raise TypeError("poly must be a vtkPolyData instance.")
    mass = vtk.vtkMassProperties()
    mass.SetInputData(poly)
    return mass.GetVolume()


def poly_surface_area(poly):
    """Compute the surface area of a vtkPolyData mesh.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        float: Surface area in voxel units squared.

    Raises:
        TypeError: If poly is not a vtkPolyData instance.
    """
    if not isinstance(poly, vtk.vtkPolyData):
        raise TypeError("poly must be a vtkPolyData instance.")
    mass = vtk.vtkMassProperties()
    mass.SetInputData(poly)
    return mass.GetSurfaceArea()


def add_sfield_to_poly(
    poly, sfield, name, dtype="float", interp="NN", mode="points"
):
    """Sample a 3-D scalar field onto vtkPolyData point/cell data.

    Each point (or cell centroid) is mapped to its scalar-field
    value via the chosen interpolation method and stored as a
    named array on the poly.

    Args:
        poly (vtk.vtkPolyData): Dataset to annotate.
        sfield (numpy.ndarray): Input 3-D scalar field.
        name (str): Name for the new scalar array.
        dtype (str): Output data type: 'float' (default) or
            'int'.
        interp (str): Interpolation strategy: 'NN' (default,
            nearest-neighbour) or 'trilin' (trilinear).
        mode (str): Whether to annotate 'points' (default) or
            'cells'.

    Raises:
        TypeError: If sfield is not a numpy array or name is not
            a string.
        ValueError: If dtype, interp, or mode is invalid.
    """
    if not isinstance(sfield, np.ndarray):
        raise TypeError("sfield must be a numpy array.")
    if not isinstance(name, str):
        raise TypeError("name must be a string.")
    if dtype not in ("float", "int"):
        raise ValueError("dtype must be 'float' or 'int'.")
    if interp not in ("NN", "trilin"):
        raise ValueError("interp must be 'NN' or 'trilin'.")
    if interp == "trilin":
        interp_func = trilin_interp
    else:
        interp_func = nn_iterp
    if mode not in ("points", "cells"):
        raise ValueError("mode must be 'points' or 'cells'.")

    cast = int if dtype == "int" else float

    if mode == "points":
        n_points = poly.GetNumberOfPoints()
        if dtype == "int":
            arr = vtk.vtkIntArray()
        else:
            arr = vtk.vtkFloatArray()
        arr.SetName(name)
        arr.SetNumberOfComponents(1)
        arr.SetNumberOfValues(n_points)
        for i in range(n_points):
            x, y, z = poly.GetPoint(i)
            arr.SetValue(i, cast(interp_func(x, y, z, sfield)))
        poly.GetPointData().AddArray(arr)
    else:
        if dtype == "int":
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
            if dtype == "int":
                values = np.zeros(shape=n_pts, dtype=int)
            else:
                values = np.zeros(shape=n_pts, dtype=float)
            for j in range(n_pts):
                x, y, z = pts.GetPoint(j)
                values[j] = interp_func(x, y, z, sfield)
            arr.SetValue(i, cast(stats.mode(values, keepdims=False).mode))
        poly.GetCellData().AddArray(arr)


def poly_mask(poly: vtk.vtkPolyData, mask: np.ndarray) -> vtk.vtkPolyData:
    """Remove vtkPolyData cells whose vertices fall outside a mask.

    Any cell with at least one vertex at a False (or out-of-bounds)
    voxel is removed.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.
        mask (numpy.ndarray): 3-D boolean mask array.

    Returns:
        vtk.vtkPolyData: Filtered dataset with masked cells
            removed.

    Raises:
        TypeError: If mask dtype is not bool.
    """
    if mask.dtype != bool:
        raise TypeError("mask must be boolean.")
    m_x, m_y, m_z = mask.shape
    del_cell_ids = vtk.vtkIdTypeArray()
    for i in range(poly.GetNumberOfCells()):
        cell = vtk.vtkGenericCell()
        poly.GetCell(i, cell)
        pts = cell.GetPoints()
        for j in range(pts.GetNumberOfPoints()):
            x, y, z = pts.GetPoint(j)
            x, y, z = int(round(x)), int(round(y)), int(round(z))
            if (
                x < 0
                or x >= m_x
                or y < 0
                or y >= m_y
                or z < 0
                or z >= m_z
                or not mask[x, y, z]
            ):
                del_cell_ids.InsertNextValue(i)
                break
    rm_filter = vtk.vtkRemovePolyData()
    rm_filter.AddInputData(poly)
    rm_filter.SetCellIds(del_cell_ids)
    rm_filter.Update()
    return rm_filter.GetOutput()


def merge_polys(poly_1, poly_2):
    """Merge two vtkPolyData objects into a single dataset.

    Args:
        poly_1 (vtk.vtkPolyData): First polygon dataset.
        poly_2 (vtk.vtkPolyData): Second polygon dataset.

    Returns:
        vtk.vtkPolyData: Combined dataset.

    Raises:
        TypeError: If either argument is not a vtkPolyData.
    """
    if not isinstance(poly_1, vtk.vtkPolyData) or not isinstance(
        poly_2, vtk.vtkPolyData
    ):
        raise TypeError("Both inputs must be vtkPolyData.")
    app_flt = vtk.vtkAppendPolyData()
    app_flt.AddInputData(poly_1)
    app_flt.AddInputData(poly_2)
    app_flt.Update()
    return app_flt.GetOutput()


def add_label_to_poly(poly, lbl, p_name, mode="cell"):
    """Assign a uniform integer label to all elements of a poly.

    Creates (or overwrites) a named integer scalar array on the
    poly's cell data, point data, or both.

    Args:
        poly (vtk.vtkPolyData): Dataset to annotate.
        lbl (int): Label value to assign.
        p_name (str): Name of the scalar array.
        mode (str): Target: 'cell' (default), 'point', or 'both'.

    Raises:
        ValueError: If mode is not 'cell', 'point', or 'both'.
        TypeError: If poly is not a vtkPolyData instance.
    """
    if mode not in ("cell", "point", "both"):
        raise ValueError("mode must be 'cell', 'point' or 'both'.")
    if not isinstance(poly, vtk.vtkPolyData):
        raise TypeError("poly must be a vtkPolyData instance.")
    lbl, p_name = int(lbl), str(p_name)

    if mode == "cell" or mode == "both":
        arr = vtk.vtkIntArray()
        n_cells = poly.GetNumberOfCells()
        arr.SetName(p_name)
        arr.SetNumberOfComponents(1)
        arr.SetNumberOfValues(n_cells)
        for i in range(n_cells):
            arr.SetValue(i, lbl)
        poly.GetCellData().AddArray(arr)
    if mode == "point" or mode == "both":
        arr = vtk.vtkIntArray()
        n_points = poly.GetNumberOfPoints()
        arr.SetName(p_name)
        arr.SetNumberOfComponents(1)
        arr.SetNumberOfValues(n_points)
        for i in range(n_points):
            arr.SetValue(i, lbl)
        poly.GetPointData().AddArray(arr)


def points_to_poly_spheres(points, rad):
    """Build a vtkPolyData with one sphere surface per input point.

    Useful for visualising particle positions as spherical glyphs.

    Args:
        points (array-like): Sequence of N 3-D coordinates with
            shape (N, 3).
        rad (float): Sphere radius in voxels.

    Returns:
        vtk.vtkPolyData: Merged polygon dataset containing one
            sphere per point.

    Raises:
        ValueError: If points is empty or coordinates are not
            3-element vectors.
    """
    if (
        not hasattr(points, "__len__")
        or len(points) == 0
        or len(points[0]) != 3
    ):
        raise ValueError(
            "points must be a non-empty sequence " "of 3-element coordinates."
        )
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
    """Compute the maximum pairwise vertex distance in a poly.

    This is an O(N²) exhaustive search; prefer :func:`poly_diam`
    for large meshes.

    Args:
        vtp (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        float: Maximum pairwise Euclidean distance; 0 if <= 1
            vertex.
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
    """Approximate the diameter of a vtkPolyData mesh.

    The diameter is estimated as twice the maximum vertex distance
    from the centre of mass, providing an O(N) approximation
    versus the exact O(N²) computation in :func:`poly_max_distance`.

    Args:
        vtp (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        float: Approximate diameter; 0 if <= 1 vertex.
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
        return 2.0 * mx


def poly_point_min_dst(poly, point, chull=False):
    """Compute the minimum distance from a 3-D point to a poly.

    Args:
        poly (vtk.vtkPolyData): Reference surface.
        point (array-like): Query 3-D coordinate.
        chull (bool): If True, extract the convex hull first to
            avoid artefacts from surface holes (default False).

    Returns:
        float: Minimum vertex-to-point distance; 0 if poly is
            empty.
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
    """Compute the centre of mass of a vtkPolyData mesh.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        numpy.ndarray: Centre of mass as a 3-element array.
    """
    cm_flt = vtk.vtkCenterOfMass()
    cm_flt.SetInputData(poly)
    cm_flt.Update()
    return np.asarray(cm_flt.GetCenter())


def convex_hull_surface(poly):
    """Extract the convex-hull surface of a vtkPolyData.

    Uses Delaunay3D followed by vtkGeometryFilter to obtain the
    outer surface.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        vtk.vtkPolyData: Convex-hull surface mesh.
    """
    convexHull = vtk.vtkDelaunay3D()
    convexHull.SetInputData(poly)
    outerSurface = vtk.vtkGeometryFilter()
    convexHull.Update()
    outerSurface.SetInputData(convexHull.GetOutput())
    outerSurface.Update()

    return outerSurface.GetOutput()


def poly_decimate(poly, dec):
    """Reduce the polygon count of a vtkPolyData mesh.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.
        dec (float): Target reduction fraction in [0, 1). For
            example 0.9 reduces the mesh to roughly 10 % of its
            original polygon count.

    Returns:
        vtk.vtkPolyData: Decimated polygon dataset.
    """
    tr_dec = vtk.vtkDecimatePro()
    tr_dec.SetInputData(poly)
    tr_dec.SetTargetReduction(dec)
    tr_dec.Update()
    return tr_dec.GetOutput()


def image_to_vti(img: np.ndarray) -> vtk.vtkImageData:
    """Convert a 2-D or 3-D numpy array to a vtkImageData object.

    Args:
        img (numpy.ndarray): Input 2-D or 3-D array.

    Returns:
        vtk.vtkImageData: Resulting VTK image with float scalars.

    Raises:
        TypeError: If img is not a numpy array.
        ValueError: If img is not 2-D or 3-D.
    """
    if not isinstance(img, np.ndarray):
        raise TypeError("img must be a numpy array.")
    n_D = len(img.shape)
    if n_D not in (2, 3):
        raise ValueError("img must be 2-D or 3-D.")
    data_type = vtk.VTK_FLOAT
    shape = img.shape

    flat_data_array = img.flatten()
    vtk_data = numpy_support.numpy_to_vtk(
        num_array=flat_data_array, deep=True, array_type=data_type
    )

    img = vtk.vtkImageData()
    img.GetPointData().SetScalars(vtk_data)
    if n_D == 2:
        img.SetDimensions(shape[0], shape[1], 1)
    else:
        img.SetDimensions(shape[0], shape[1], shape[2])

    return img


def save_vti(img: vtk.vtkImageData, fname: str):
    """Write a vtkImageData object to a .vti or .vtk file.

    Args:
        img (vtk.vtkImageData): Input VTK image.
        fname (str): Output file path ending with '.vtk' or
            '.vti'.

    Raises:
        TypeError: If img is not a vtkImageData or fname is not a
            string.
        ValueError: If fname does not end with '.vtk' or '.vti'.
    """
    if not isinstance(img, vtk.vtkImageData):
        raise TypeError("img must be a vtkImageData instance.")
    if not isinstance(fname, str):
        raise TypeError("fname must be a string.")
    if not (fname.endswith(".vtk") or fname.endswith(".vti")):
        raise ValueError("fname must end with '.vtk' or '.vti'.")

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(fname)
    writer.SetInputData(img)
    if writer.Write() != 1:
        logger.warning("Failed to write .vti file: %s", fname)
