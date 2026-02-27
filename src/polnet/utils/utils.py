"""Common mathematical and volumetric utility functions.

Provides density normalisation, uniform sampling on spheres,
sub-volume insertion into tomograms, and other shared operations
used across the Polnet simulation modules.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import errno
import math
import os
import shutil
import stat

import numpy as np
import skimage
import vtk
from vtkmodules.util import numpy_support

from ..logging_conf import _LOGGER as logger

PI_2 = 2 * np.pi
VTK_RAY_TOLERANCE = 0.000001  # 0.001


def gen_uni_s2_sample(center, rad):
    """Generate a coordinate drawn uniformly from a sphere surface.

    Args:
        center (array-like): Sphere centre as a 3-element vector.
        rad (float): Sphere radius.

    Returns:
        numpy.ndarray: A 3-D coordinate on the sphere surface.
    """
    X = np.random.randn(1, 3)[0]
    norm = rad / np.linalg.norm(X)
    X *= norm
    return X + center


def points_distance(a, b):
    """Compute the Euclidean distance between two points.

    Args:
        a (numpy.ndarray): First point as a 3-element array.
        b (numpy.ndarray): Second point as a 3-element array.

    Returns:
        float: The Euclidean distance d(a, b).
    """
    hold = b - a
    return math.sqrt((hold * hold).sum())


def iso_surface(tomo, th, flp=None, closed=False, normals=None):
    """Extract an isosurface from a 3-D volume via Marching Cubes.

    Args:
        tomo (numpy.ndarray): Input 3-D volume.
        th (float): Isovalue threshold.
        flp (int, optional): Axis to flip before extraction
            (0, 1, or 2). Defaults to None.
        closed (bool): If True, pads the volume to guarantee a
            closed surface (valid only for binary inputs).
            Defaults to False.
        normals (str, optional): Normal orientation — 'inwards',
            'outwards', or None (default, no reorientation).

    Returns:
        vtk.vtkPolyData: A triangle-only polygon dataset.

    Raises:
        RuntimeError: If closed is True but the output surface
            is not closed.
    """

    march = vtk.vtkMarchingCubes()
    tomo_vtk = numpy_to_vti(tomo)
    if closed:
        padder = vtk.vtkImageConstantPad()
        padder.SetInputData(tomo_vtk)
        padder.SetConstant(0)
        padder.SetOutputWholeExtent(
            -1, tomo.shape[0], -1, tomo.shape[1], -1, tomo.shape[2]
        )
        padder.Update()
        tomo_vtk = padder.GetOutput()

    if flp is not None:
        flp_i = int(flp)
        if (flp_i >= 0) and (flp_i <= 3):
            fliper = vtk.vtkImageFlip()
            fliper.SetFilteredAxis(flp_i)
            fliper.SetInputData(tomo_vtk)
            fliper.Update()
            tomo_vtk = fliper.GetOutput()

    march.SetInputData(tomo_vtk)
    march.SetValue(0, th)
    march.Update()
    hold_poly = march.GetOutput()

    hold_poly = poly_filter_triangles(hold_poly)

    if normals is not None:
        orienter = vtk.vtkPolyDataNormals()
        orienter.SetInputData(hold_poly)
        orienter.AutoOrientNormalsOn()
        if normals == "inwards":
            orienter.FlipNormalsOn()
        orienter.Update()
        hold_poly = orienter.GetOutput()

    if closed and (not is_closed_surface(hold_poly)):
        raise RuntimeError

    return hold_poly


def is_closed_surface(poly):
    """Check whether a vtkPolyData surface is closed (watertight).

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        bool: True if the surface is closed, False otherwise.
    """
    selector = vtk.vtkSelectEnclosedPoints()
    selector.CheckSurfaceOn()
    selector.SetSurfaceData(poly)
    if selector.GetCheckSurface() > 0:
        return True
    else:
        return False


def poly_filter_triangles(poly):
    """Filter a vtkPolyData to retain only triangle cells.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.

    Returns:
        vtk.vtkPolyData: Copy containing only triangles.
    """
    cut_tr = vtk.vtkTriangleFilter()
    cut_tr.SetInputData(poly)
    cut_tr.PassVertsOff()
    cut_tr.PassLinesOff()
    cut_tr.Update()
    return cut_tr.GetOutput()


def numpy_to_vti(array, spacing=None):
    """Convert a 3-D numpy array to a vtkImageData object.

    Args:
        array (numpy.ndarray): Input 3-D array.
        spacing (list[float] | None): Voxel spacing for each axis
            (default ``[1, 1, 1]``).

    Returns:
        vtk.vtkImageData: The resulting VTK image object.
    """
    if spacing is None:
        spacing = [1, 1, 1]
    if not isinstance(array, np.ndarray):
        raise TypeError("Input array must be a 3-D numpy array.")

    array_1d = numpy_support.numpy_to_vtk(
        num_array=np.reshape(array, -1, order="F"),
        deep=True,
        array_type=vtk.VTK_FLOAT,
    )

    nx, ny, nz = array.shape
    image = vtk.vtkImageData()
    image.SetSpacing(spacing)
    image.SetDimensions(nx, ny, nz)
    image.AllocateScalars(vtk.VTK_FLOAT, 1)
    image.GetPointData().SetScalars(array_1d)

    return image


def get_sub_copy(tomo, sub_pt, sub_shape):
    """Extract a sub-volume centred on a point from a tomogram.

    Out-of-bounds regions are zero-filled.

    Args:
        tomo (numpy.ndarray): Input 3-D tomogram.
        sub_pt (array-like): Centre voxel of the sub-volume
            (x, y, z).
        sub_shape (tuple[int, int, int]): Output sub-volume shape;
            all dimensions must be even.

    Returns:
        numpy.ndarray: Copy of the sub-volume with the requested
            shape.
    """

    nx, ny, nz = sub_shape[0], sub_shape[1], sub_shape[2]
    mx, my, mz = tomo.shape[0], tomo.shape[1], tomo.shape[2]
    mx1, my1, mz1 = mx - 1, my - 1, mz - 1
    hl_x, hl_y, hl_z = int(nx * 0.5), int(ny * 0.5), int(nz * 0.5)
    x, y, z = (
        int(round(sub_pt[0])),
        int(round(sub_pt[1])),
        int(round(sub_pt[2])),
    )

    off_l_x, off_l_y, off_l_z = x - hl_x, y - hl_y, z - hl_z
    off_h_x, off_h_y, off_h_z = x + hl_x, y + hl_y, z + hl_z
    dif_l_x, dif_l_y, dif_l_z = 0, 0, 0
    dif_h_x, dif_h_y, dif_h_z = nx, ny, nz
    if off_l_x < 0:
        dif_l_x = abs(off_l_x)
        off_l_x = 0
    if off_l_y < 0:
        dif_l_y = abs(off_l_y)
        off_l_y = 0
    if off_l_z < 0:
        dif_l_z = abs(off_l_z)
        off_l_z = 0
    if off_h_x >= mx:
        dif_h_x = nx - off_h_x + mx1
        off_h_x = mx1
    if off_h_y >= my:
        dif_h_y = ny - off_h_y + my1
        off_h_y = my1
    if off_h_z >= mz:
        dif_h_z = nz - off_h_z + mz1
        off_h_z = mz1

    # Make the subvolume copy
    hold_sv = np.zeros(shape=sub_shape, dtype=tomo.dtype)
    hold_sv[dif_l_x:dif_h_x, dif_l_y:dif_h_y, dif_l_z:dif_h_z] = tomo[
        off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z
    ]

    return hold_sv


def insert_svol_tomo(svol, tomo, sub_pt, merge="max"):
    """Stamp a sub-volume into a target tomogram at a given centre.

    Voxels that fall outside the tomogram bounds are silently
    clipped.

    Args:
        svol (numpy.ndarray): Input sub-volume to stamp.
        tomo (numpy.ndarray): Target tomogram modified in place.
        sub_pt (array-like): Centre voxel of the insertion point
            (x, y, z).
        merge (str): Blending strategy — 'max' (default), 'min',
            'sum', 'insert', or 'and'.
    """

    sub_shape = svol.shape
    nx, ny, nz = sub_shape[0], sub_shape[1], sub_shape[2]
    mx, my, mz = tomo.shape[0], tomo.shape[1], tomo.shape[2]
    mx1, my1, mz1 = mx - 1, my - 1, mz - 1
    hl_x, hl_y, hl_z = int(nx * 0.5), int(ny * 0.5), int(nz * 0.5)
    x, y, z = (
        int(round(sub_pt[0])),
        int(round(sub_pt[1])),
        int(round(sub_pt[2])),
    )

    # Compute bounding restriction
    off_l_x, off_l_y, off_l_z = x - hl_x, y - hl_y, z - hl_z
    off_h_x, off_h_y, off_h_z = x + hl_x, y + hl_y, z + hl_z
    dif_l_x, dif_l_y, dif_l_z = 0, 0, 0
    dif_h_x, dif_h_y, dif_h_z = nx, ny, nz
    if off_l_x < 0:
        dif_l_x = abs(off_l_x)
        off_l_x = 0
    if off_l_y < 0:
        dif_l_y = abs(off_l_y)
        off_l_y = 0
    if off_l_z < 0:
        dif_l_z = abs(off_l_z)
        off_l_z = 0
    if off_h_x >= mx:
        dif_h_x = nx - off_h_x + mx1
        off_h_x = mx1
    if off_h_y >= my:
        dif_h_y = ny - off_h_y + my1
        off_h_y = my1
    if off_h_z >= mz:
        dif_h_z = nz - off_h_z + mz1
        off_h_z = mz1
    if off_l_x > off_h_x:
        off_h_x = off_l_x
    if off_l_y > off_h_y:
        off_h_y = off_l_y
    if off_l_z > off_h_z:
        off_h_z = off_l_z
    if dif_l_x > dif_h_x:
        dif_h_x = dif_l_x
    if dif_l_y > dif_h_y:
        dif_h_y = dif_l_y
    if dif_l_z > dif_h_z:
        dif_h_z = dif_l_z
    sz_svol = [dif_h_x - dif_l_x, dif_h_y - dif_l_y, dif_h_z - dif_l_z]
    sz_off = [off_h_x - off_l_x, off_h_y - off_l_y, off_h_z - off_l_z]
    if (sz_svol[0] > sz_off[0]) and (sz_svol[0] > 1):
        dif_h_x -= 1
    if (sz_svol[1] > sz_off[1]) and (sz_svol[1] > 1):
        dif_h_y -= 1
    if (sz_svol[2] > sz_off[2]) and (sz_svol[2] > 1):
        dif_h_z -= 1

    # Modify the input tomogram
    if merge == "insert":
        tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z] = svol[
            dif_l_x:dif_h_x, dif_l_y:dif_h_y, dif_l_z:dif_h_z
        ]
    elif merge == "sum":
        tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z] += svol[
            dif_l_x:dif_h_x, dif_l_y:dif_h_y, dif_l_z:dif_h_z
        ]
    elif merge == "min":
        tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z] = np.minimum(
            svol[dif_l_x:dif_h_x, dif_l_y:dif_h_y, dif_l_z:dif_h_z],
            tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z],
        )
    elif merge == "max":
        tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z] = np.maximum(
            svol[dif_l_x:dif_h_x, dif_l_y:dif_h_y, dif_l_z:dif_h_z],
            tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z],
        )
    elif merge == "and":
        tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z] = (
            np.logical_and(
                svol[dif_l_x:dif_h_x, dif_l_y:dif_h_y, dif_l_z:dif_h_z],
                tomo[off_l_x:off_h_x, off_l_y:off_h_y, off_l_z:off_h_z],
            )
        )


# Applies a linear mapping to the input array for getting an array in the specified range
def lin_map(array, lb=0, ub=1):
    """Linearly rescale an array to the range [lb, ub].

    Args:
        array (numpy.ndarray): Input array to remap.
        lb (float): Lower output bound (default 0).
        ub (float): Upper output bound (default 1).

    Returns:
        numpy.ndarray: The remapped array with values in [lb, ub].
    """
    a = np.max(array)
    i = np.min(array)
    den = a - i
    if den == 0:
        return array
    m = (ub - lb) / den
    c = ub - m * a
    return m * array + c


def wrap_angle(ang, deg=True):
    """Wrap an angle to the range (-180, 180] or (-pi, pi].

    Args:
        ang (float | numpy.ndarray): Input angle or array of
            angles.
        deg (bool): If True (default), ang is in degrees;
            otherwise in radians.

    Returns:
        float | numpy.ndarray: Wrapped angle(s) in the same units
            as the input.
    """
    if deg:
        phase = ((-ang + 180.0) % (2.0 * 180.0) - 180.0) * -1.0
    else:
        phase = ((-ang + np.pi) % (2.0 * np.pi) - np.pi) * -1.0
    return phase


def point_to_poly(point, normal=None, n_name="n_normal"):
    """Convert a single 3-D point into a one-vertex vtkPolyData.

    Args:
        point (array-like): 3-tuple with point coordinates.
        normal (array-like, optional): 3-tuple with the surface
            normal stored as a point-data array. Defaults to None.
        n_name (str): Name for the normal array
            (default 'n_normal').

    Returns:
        vtk.vtkPolyData: A vtkPolyData containing a single vertex.
    """
    poly = vtk.vtkPolyData()
    p_points = vtk.vtkPoints()
    p_cells = vtk.vtkCellArray()
    p_points.InsertNextPoint(point)
    p_cells.InsertNextCell(1)
    p_cells.InsertCellPoint(0)
    poly.SetPoints(p_points)
    poly.SetVerts(p_cells)
    if normal is not None:
        p_norm = vtk.vtkFloatArray()
        p_norm.SetName(n_name)
        p_norm.SetNumberOfComponents(3)
        p_norm.InsertTuple(0, normal)
        poly.GetPointData().AddArray(p_norm)
    return poly


def density_norm(tomo, mask=None, inv=True):
    """Normalise a tomogram to zero mean and unit standard deviation.

    Normalisation formula: (I(x,y,z) - mean) / std.

    Args:
        tomo (numpy.ndarray): Input 3-D tomogram.
        mask (numpy.ndarray, optional): Boolean mask selecting the
            region used for statistics. None (default) uses the
            whole tomogram.
        inv (bool): If True (default), values are inverted
            (multiplied by -1) before normalisation.

    Returns:
        numpy.ndarray: Normalised tomogram as float32.
    """
    if mask is None:
        mask = np.ones(shape=tomo.shape, dtype=bool)

    if inv:
        hold_tomo = -1.0 * tomo
    else:
        hold_tomo = tomo

    stat_tomo = hold_tomo[mask > 0]
    mn, st = stat_tomo.mean(), stat_tomo.std()

    tomo_out = np.zeros(shape=tomo.shape, dtype=np.float32)
    if st > 0:
        tomo_out = (hold_tomo - mn) / st
    else:
        logger.warning(
            "density_norm: std=%s, returning zero-filled tomogram", st
        )

    return tomo_out


def trilin_interp(x, y, z, tomogram):
    """Trilinear interpolation of a scalar field at a sub-voxel point.

    Args:
        x (float): X coordinate (within tomogram bounds).
        y (float): Y coordinate (within tomogram bounds).
        z (float): Z coordinate (within tomogram bounds).
        tomogram (numpy.ndarray): Input 3-D scalar field.

    Returns:
        float: Interpolated value at (x, y, z).

    Raises:
        TypeError: If tomogram is not a 3-D numpy array.
        ValueError: If any coordinate is out of bounds.
    """

    if not isinstance(tomogram, np.ndarray) or tomogram.ndim != 3:
        raise TypeError("tomogram must be a 3-D numpy array.")
    xc = int(math.ceil(x))
    yc = int(math.ceil(y))
    zc = int(math.ceil(z))
    xf = int(math.floor(x))
    yf = int(math.floor(y))
    zf = int(math.floor(z))
    if (
        xc >= tomogram.shape[0]
        or yc >= tomogram.shape[1]
        or zc >= tomogram.shape[2]
        or xf < 0
        or yf < 0
        or zf < 0
    ):
        raise ValueError("Coordinates are out of tomogram bounds.")

    v000 = float(tomogram[xf, yf, zf])
    v100 = float(tomogram[xc, yf, zf])
    v010 = float(tomogram[xf, yc, zf])
    v001 = float(tomogram[xf, yf, zc])
    v101 = float(tomogram[xc, yf, zc])
    v011 = float(tomogram[xf, yc, zc])
    v110 = float(tomogram[xc, yc, zf])
    v111 = float(tomogram[xc, yc, zc])

    xn = x - xf
    yn = y - yf
    zn = z - zf
    x1 = 1 - xn
    y1 = 1 - yn
    z1 = 1 - zn

    return (
        (v000 * x1 * y1 * z1)
        + (v100 * xn * y1 * z1)
        + (v010 * x1 * yn * z1)
        + (v001 * x1 * y1 * zn)
        + (v101 * xn * y1 * zn)
        + (v011 * x1 * yn * zn)
        + (v110 * xn * yn * z1)
        + (v111 * xn * yn * zn)
    )


def nn_iterp(x, y, z, tomogram):
    """Nearest-neighbour interpolation of a scalar field at a point.

    Args:
        x (float): X coordinate (within tomogram bounds).
        y (float): Y coordinate (within tomogram bounds).
        z (float): Z coordinate (within tomogram bounds).
        tomogram (numpy.ndarray): Input 3-D scalar field.

    Returns:
        float: Scalar value at the nearest voxel to (x, y, z).

    Raises:
        TypeError: If tomogram is not a 3-D numpy array.
        ValueError: If any coordinate is out of bounds.
    """

    if not isinstance(tomogram, np.ndarray) or tomogram.ndim != 3:
        raise TypeError("tomogram must be a 3-D numpy array.")
    xc = int(math.ceil(x))
    yc = int(math.ceil(y))
    zc = int(math.ceil(z))
    xf = int(math.floor(x))
    yf = int(math.floor(y))
    zf = int(math.floor(z))
    if (
        xc >= tomogram.shape[0]
        or yc >= tomogram.shape[1]
        or zc >= tomogram.shape[2]
        or xf < 0
        or yf < 0
        or zf < 0
    ):
        raise ValueError("Coordinates are out of tomogram bounds.")

    point = np.asarray((x, y, z))
    X, Y, Z = np.meshgrid(
        range(xf, xc + 1), range(yf, yc + 1), range(zf, zc + 1), indexing="ij"
    )
    X, Y, Z = X.flatten(), Y.flatten(), Z.flatten()
    min_point = np.asarray((X[0], Y[0], Z[0]))
    hold = point - min_point
    min_dist = np.sqrt((hold * hold).sum())
    for i in range(1, len(X)):
        hold_point = np.asarray((X[i], Y[i], Z[i]))
        hold = point - hold_point
        hold_dist = np.sqrt((hold * hold).sum())
        if hold_dist < min_dist:
            min_point = hold_point
            min_dist = hold_dist

    return tomogram[min_point[0], min_point[1], min_point[2]]


def poly_threshold(poly, p_name, mode="points", low_th=None, hi_th=None):
    """Threshold a vtkPolyData by the values of a named property.

    Args:
        poly (vtk.vtkPolyData): Input polygon dataset.
        p_name (str): Name of the scalar property to threshold on.
        mode (str): Whether the property belongs to 'points'
            (default) or 'cells'.
        low_th (float, optional): Lower threshold bound. Defaults
            to the property minimum.
        hi_th (float, optional): Upper threshold bound. Defaults
            to the property maximum.

    Returns:
        vtk.vtkPolyData: Thresholded polygon dataset.

    Raises:
        ValueError: If mode is invalid or p_name is not found.
    """

    prop = None
    if mode not in ("points", "cells"):
        raise ValueError("mode must be 'points' or 'cells'.")
    if mode == "points":
        n_arrays = poly.GetPointData().GetNumberOfArrays()
        for i in range(n_arrays):
            if p_name == poly.GetPointData().GetArrayName(i):
                prop = poly.GetPointData().GetArray(p_name)
                break
    else:
        n_arrays = poly.GetCellData().GetNumberOfArrays()
        for i in range(n_arrays):
            if p_name == poly.GetCellData().GetArrayName(i):
                prop = poly.GetCellData().GetArray(p_name)
                break
    if prop is None:
        raise ValueError(f"Property '{p_name}' not found in poly.")
    rg_low, rg_hi = None, None
    if (low_th is None) or (hi_th is None):
        rg_low, rg_hi = prop.GetRange()
    if low_th is None:
        low_th = rg_low
    if hi_th is None:
        hi_th = rg_hi

    th_flt = vtk.vtkThreshold()
    th_flt.SetInputData(poly)
    if mode == "cells":
        th_flt.SetInputArrayToProcess(
            0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_CELLS, p_name
        )
    else:
        th_flt.SetInputArrayToProcess(
            0, 0, 0, vtk.vtkDataObject.FIELD_ASSOCIATION_POINTS, p_name
        )
    th_flt.SetLowerThreshold(low_th)
    th_flt.SetUpperThreshold(hi_th)
    th_flt.AllScalarsOff()
    th_flt.Update()

    surf_flt = vtk.vtkDataSetSurfaceFilter()
    surf_flt.SetInputData(th_flt.GetOutput())
    surf_flt.Update()

    return surf_flt.GetOutput()


def gen_six_connectivity_mask():
    """Generate a 6-connectivity structuring element.

    Returns:
        numpy.ndarray: A (3, 3, 3) boolean array with the six
            face-connected neighbours set to True.
    """
    mask = np.zeros(shape=(3, 3, 3), dtype=bool)
    mask[1, 1, 0] = True
    mask[0, 1, 1] = True
    mask[1, 1, 1] = True
    mask[2, 1, 1] = True
    mask[1, 0, 1] = True
    mask[1, 2, 1] = True
    mask[1, 1, 2] = True
    return mask


def clean_dir(dir_path):
    """Remove all contents of a directory while preserving it.

    Args:
        dir_path (str | Path): Path to the directory to clean.
    """
    for root, dirs, files in os.walk(dir_path):
        for f in files:
            f_name = os.path.join(root, f)
            try:
                os.unlink(f_name)
            except OSError as err:
                logger.debug("Error deleting file %s: %s", f_name, err)
                if err.errno != errno.EBUSY:
                    os.chmod(f_name, stat.S_IWRITE)
                    os.remove(f_name)
        for d in dirs:
            d_name = os.path.join(root, d)
            try:
                shutil.rmtree(d_name)
            except OSError as err:
                logger.debug("Error deleting directory %s: %s", d_name, err)
                if err.errno != errno.EBUSY:
                    os.chmod(d_name, stat.S_IWRITE)
                    os.remove(d_name)


def vol_cube(vol, off=0):
    """Embed a 3-D volume in a cubic array padded to its max side.

    Args:
        vol (numpy.ndarray): Input 3-D array.
        off (int): Extra padding voxels added symmetrically on all
            sides (default 0).

    Returns:
        numpy.ndarray: Cubic array with side length
            max(vol.shape) + off.

    Raises:
        TypeError: If vol is not a 3-D numpy array.
        ValueError: If off is negative.
    """
    if not isinstance(vol, np.ndarray) or vol.ndim != 3:
        raise TypeError("vol must be a 3-D numpy array.")
    if off < 0:
        raise ValueError("off must be >= 0.")
    dim_max = np.argmax(vol.shape)
    cube_dim = vol.shape[dim_max]
    out_vol = np.zeros(shape=(cube_dim, cube_dim, cube_dim), dtype=vol.dtype)
    if dim_max == 0:
        off_ly, off_lz = (cube_dim - vol.shape[1]) // 2, (
            cube_dim - vol.shape[2]
        ) // 2
        off_hy, off_hz = off_ly + vol.shape[1], off_lz + vol.shape[2]
        out_vol[:, off_ly:off_hy, off_lz:off_hz] = vol
    elif dim_max == 1:
        off_lx, off_lz = (cube_dim - vol.shape[0]) // 2, (
            cube_dim - vol.shape[2]
        ) // 2
        off_hx, off_hz = off_lx + vol.shape[0], off_lz + vol.shape[2]
        out_vol[off_lx:off_hx, :, off_lz:off_hz] = vol
    else:
        off_lx, off_ly = (cube_dim - vol.shape[0]) // 2, (
            cube_dim - vol.shape[1]
        ) // 2
        off_hx, off_hy = off_lx + vol.shape[0], off_ly + vol.shape[1]
        out_vol[off_lx:off_hx, off_ly:off_hy, :] = vol
    if off > 0:
        off_2 = off // 2
        off_vol = np.zeros(
            shape=(cube_dim + off, cube_dim + off, cube_dim + off),
            dtype=vol.dtype,
        )
        off_vol[
            off_2 : off_2 + cube_dim,
            off_2 : off_2 + cube_dim,
            off_2 : off_2 + cube_dim,
        ] = out_vol
        out_vol = off_vol
    return out_vol


def gen_sphere_mask(shape, radius, center=None):
    """Generate a binary spherical mask as a 3-D boolean array.

    Args:
        shape (array-like): Output array shape as (nx, ny, nz).
        radius (float): Sphere radius in voxels.
        center (array-like, optional): Sphere centre coordinates.
            Defaults to the array centre.

    Returns:
        numpy.ndarray: Boolean array with True inside the sphere.

    Raises:
        ValueError: If shape or center have invalid dimensions or
            non-positive values.
    """
    if not hasattr(shape, "__len__") or len(shape) != 3:
        raise ValueError("shape must have exactly 3 elements.")
    if shape[0] <= 0 or shape[1] <= 0 or shape[2] <= 0:
        raise ValueError("All shape dimensions must be > 0.")
    if center is None:
        center = 0.5 * np.asarray(shape)
    else:
        if not hasattr(center, "__len__") or len(center) != 3:
            raise ValueError("center must have exactly 3 elements.")
        if center[0] <= 0 or center[1] <= 0 or center[2] <= 0:
            raise ValueError("All center values must be > 0.")
        center = np.asarray(center)

    x_mat, y_mat, z_mat = np.meshgrid(
        range(shape[0]), range(shape[1]), range(shape[2]), indexing="ij"
    )
    x_mat = x_mat.astype(np.float32) - center[0]
    y_mat = y_mat.astype(np.float32) - center[1]
    z_mat = z_mat.astype(np.float32) - center[2]
    return x_mat * x_mat + y_mat * y_mat + z_mat * z_mat <= radius * radius


def tomo_crop_non_zeros(tomo):
    """Crop a tomogram to the tight bounding box of non-zero voxels.

    Args:
        tomo (numpy.ndarray): Input 3-D tomogram.

    Returns:
        numpy.ndarray: Cropped tomogram containing only non-zero
            regions.
    """
    coords = np.argwhere(tomo > 0)
    x_min, y_min, z_min = coords.min(axis=0)
    x_max, y_max, z_max = coords.max(axis=0)
    return tomo[x_min : x_max + 1, y_min : y_max + 1, z_min : z_max + 1]


def connectivity_analysis(tomo, th, connectivity=None):
    """Remove connected components smaller than a size threshold.

    Args:
        tomo (numpy.ndarray): Input 3-D binary (or positive)
            tomogram.
        th (int): Minimum voxel count for a region to be kept.
        connectivity (int, optional): Passed to
            skimage.measure.label. Defaults to None (full).

    Returns:
        numpy.ndarray: Integer array where voxel values equal the
            size of their connected region (0 for removed).
    """
    tomo_mb, num_lbls = skimage.measure.label(
        tomo > 0, connectivity=connectivity, return_num=True
    )
    tomo_sz = np.zeros(shape=tomo_mb.shape, dtype=np.int32)
    for lbl in range(1, num_lbls + 1):
        ids = tomo_mb == lbl
        feat_sz = ids.sum()
        if feat_sz >= th:
            tomo_sz[ids] = feat_sz
    return tomo_sz


MAX_TRIES_EXP = int(1e6)


def gen_bounded_exp(mean, lb, ub):
    """Draw a random number from a bounded exponential distribution.

    Args:
        mean (float): Mean of the underlying exponential
            distribution (1/lambda).
        lb (float): Lower bound for acceptance.
        ub (float): Upper bound for acceptance.

    Returns:
        float: A sampled value in [lb, ub].

    Raises:
        RuntimeError: If a valid sample is not drawn within
            MAX_TRIES_EXP attempts.
    """
    hold = np.random.exponential(scale=mean)
    count = 1
    while ((hold < lb) or (hold > ub)) and (count < MAX_TRIES_EXP):
        hold = np.random.exponential(scale=mean)
        count += 1
    if count >= MAX_TRIES_EXP:
        raise RuntimeError
    else:
        return hold
