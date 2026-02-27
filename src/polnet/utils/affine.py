"""Affine transformation utilities for 3-D cryo-ET volumes.

Provides rotation matrices, quaternion operations, and coordinate
transformations used throughout the Polnet simulation pipeline.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import math
import random

import numpy as np
import vtk
from scipy import ndimage as spnd
from scipy.spatial.transform import Rotation as spR

from .utils import wrap_angle

PI_2 = 2 * np.pi
LOW_VALUE = 0.000001
# Super-fibonacci constants
PHI = np.sqrt(2.0)
PSI = 1.533751168755204288118041


def vector_module(v):
    r"""Compute the Euclidean module (L2 norm) of a vector.

    Args:
        v (numpy.ndarray): One-dimensional input vector.

    Returns:
        float: The computed module \|v\|.
    """
    return math.sqrt((v * v).sum())


def poly_rotate_wxyz(in_vtp, w, x, y, z):
    """Apply a rigid rotation to a vtkPolyData from angle-axis input.

    Args:
        in_vtp (vtk.vtkPolyData): Input polygon dataset.
        w (float): Rotation angle in degrees.
        x (float): Rotation axis X component.
        y (float): Rotation axis Y component.
        z (float): Rotation axis Z component.

    Returns:
        vtk.vtkPolyData: The rotated polygon dataset.
    """
    rot_tf = vtk.vtkTransform()
    rot_tf.RotateWXYZ(w, x, y, z)
    tf_rot = vtk.vtkTransformPolyDataFilter()
    tf_rot.SetInputData(in_vtp)
    tf_rot.SetTransform(rot_tf)
    tf_rot.Update()
    return tf_rot.GetOutput()


def poly_translate(in_vtp, t):
    """Apply a rigid translation to a vtkPolyData.

    Args:
        in_vtp (vtk.vtkPolyData): Input polygon dataset.
        t (array-like): Translation vector with 3 elements [tx, ty, tz].

    Returns:
        vtk.vtkPolyData: The translated polygon dataset.

    Raises:
        ValueError: If t does not have exactly 3 elements.
    """

    if not hasattr(t, "__len__") or len(t) != 3:
        raise ValueError("Translation vector t must have 3 elements.")

    box_tr = vtk.vtkTransform()
    box_tr.Translate(t[0], t[1], t[2])
    tr_box = vtk.vtkTransformPolyDataFilter()
    tr_box.SetInputData(in_vtp)
    tr_box.SetTransform(box_tr)
    tr_box.Update()
    return tr_box.GetOutput()


def poly_scale(in_vtp, s):
    """Apply a uniform scaling transformation to a vtkPolyData.

    Args:
        in_vtp (vtk.vtkPolyData): Input polygon dataset.
        s (float): Uniform scaling factor applied to all axes.

    Returns:
        vtk.vtkPolyData: The scaled polygon dataset.
    """
    box_sc = vtk.vtkTransform()
    box_sc.Scale(s, s, s)
    tr_box = vtk.vtkTransformPolyDataFilter()
    tr_box.SetInputData(in_vtp)
    tr_box.SetTransform(box_sc)
    tr_box.Update()
    return tr_box.GetOutput()


def gen_rand_unit_quaternion():
    """Generate a uniformly random unit quaternion on SO(3).

    Uses the Kuffner (2004) method for uniform SO(3) sampling.

    Reference:
        Kuffner, J.J. (2004). Effective sampling and distance metrics
        for 3D rigid body path planning. ICRA 2004.

    Returns:
        numpy.ndarray: Unit quaternion of shape (4,) as
            (w, x, y, z) = (cos(θ/2), vx·sin(θ/2), vy·sin(θ/2),
            vz·sin(θ/2)), where (vx, vy, vz) is the rotation axis
            and θ is the angle.
    """
    s = random.random()
    sigma_1, sigma_2 = math.sqrt(1 - s), math.sqrt(s)
    theta_1, theta_2 = PI_2 * random.random(), PI_2 * random.random()
    w = math.cos(theta_2) * sigma_2
    x = math.sin(theta_1) * sigma_1
    y = math.cos(theta_1) * sigma_1
    z = math.sin(theta_2) * sigma_2
    return np.asarray((w, x, y, z))


def quat_to_angle_axis(qw, qx, qy, qz, deg=True):
    """Convert a unit quaternion to angle-axis representation.

    Args:
        qw (float): Quaternion real part w.
        qx (float): Quaternion imaginary part x.
        qy (float): Quaternion imaginary part y.
        qz (float): Quaternion imaginary part z.
        deg (bool): If True (default), the returned angle is in
            degrees; otherwise in radians.

    Returns:
        tuple[float, numpy.ndarray]: A 2-tuple (angle, axis) where
            angle is a scalar and axis is a unit vector of shape (3,).
    """
    x = np.asarray((qx, qy, qz))
    norm = vector_module(x)
    if norm <= 0:
        return 0, np.asarray((0.0, 0.0, 0.0))
    if deg:
        ang_rad = 2.0 * math.atan2(norm, qw)
        ang = math.degrees(ang_rad)
        return ang, x / math.sin(0.5 * ang_rad)
    else:
        ang = 2.0 * math.atan2(norm, qw)
        return ang, x / math.sin(0.5 * ang)


def angle_axis_to_quat(ang, x, y, z, deg=True):
    """Convert an angle-axis representation to a unit quaternion.

    Args:
        ang (float): Rotation angle.
        x (float): Rotation axis X component.
        y (float): Rotation axis Y component.
        z (float): Rotation axis Z component.
        deg (bool): If True (default), ang is in degrees;
            otherwise in radians.

    Returns:
        numpy.ndarray: Unit quaternion of shape (4,) as
            (w, x, y, z).
    """
    ax = np.asarray((x, y, z), dtype=float)
    ax /= vector_module(ax)
    hold_ang = 0.5 * wrap_angle(ang, deg=deg)
    if deg:
        hold_ang = math.radians(hold_ang)
    ca, sa = math.cos(hold_ang), math.sin(hold_ang)
    return np.asarray((ca, ax[0] * sa, ax[1] * sa, ax[2] * sa), dtype=float)


def rot_vect_quat(v, q):
    """Rotate a vector by a quaternion using Rodrigues' formula.

    Args:
        v (array-like): Input 3-D vector.
        q (array-like): Unit quaternion of shape (4,) as
            (w, x, y, z).

    Returns:
        array-like: The rotated vector with the same type as v.
    """

    ang, ax = quat_to_angle_axis(q[0], q[1], q[2], q[3], deg=False)

    k = np.asarray((ax[0], ax[1], ax[2]))
    mod_k = math.sqrt((k * k).sum())
    if mod_k > 0:
        k /= mod_k
        vc = np.asarray(v)

        # Rodrigues formula
        cos_ang, sin_ang = math.cos(ang), math.sin(ang)
        return (
            vc * cos_ang
            + np.cross(k, vc) * sin_ang
            + k * np.dot(k, vc) * (1.0 - cos_ang)
        )
    else:
        return v


def quat_to_mat(q):
    """Convert a unit quaternion to a 3x3 rotation matrix.

    Args:
        q (array-like): Quaternion of shape (4,) as
            (q0, q1, q2, q3).

    Returns:
        numpy.ndarray: A (3, 3) rotation matrix that converts a
            point from the local reference frame to the global
            reference frame.
    """
    q0, q1, q2, q3 = q[0], q[1], q[2], q[3]

    r00 = 2 * (q0 * q0 + q1 * q1) - 1
    r01 = 2 * (q1 * q2 - q0 * q3)
    r02 = 2 * (q1 * q3 + q0 * q2)

    r10 = 2 * (q1 * q2 + q0 * q3)
    r11 = 2 * (q0 * q0 + q2 * q2) - 1
    r12 = 2 * (q2 * q3 - q0 * q1)

    r20 = 2 * (q1 * q3 - q0 * q2)
    r21 = 2 * (q2 * q3 + q0 * q1)
    r22 = 2 * (q0 * q0 + q3 * q3) - 1

    return np.array([[r00, r01, r02], [r10, r11, r12], [r20, r21, r22]])


def tomo_rotate(
    tomo,
    q,
    center=None,
    active=True,
    order=3,
    mode="constant",
    cval=0.0,
    prefilter=True,
):
    """Apply a quaternion-encoded rotation to a 3-D volume.

    Args:
        tomo (numpy.ndarray): Input 3-D volume.
        q (array-like): Unit quaternion of shape (4,) encoding
            the rotation.
        center (numpy.ndarray, optional): Centre of rotation as a
            3-element array. Defaults to the tomogram centre.
        active (bool): If True (default), applies an active
            rotation; otherwise a passive (inverse) rotation.
        order (int): Interpolation spline order for
            scipy.ndimage.affine_transform (default 3).
        mode (str): Border mode for
            scipy.ndimage.affine_transform (default 'constant').
        cval (float): Fill value for mode='constant'
            (default 0.0).
        prefilter (bool): Pre-filter flag for
            scipy.ndimage.affine_transform (default True).

    Returns:
        numpy.ndarray: The rotated volume with the same shape
            as tomo.

    Raises:
        TypeError: If tomo is not a 3-D numpy array.
        ValueError: If q does not have 4 elements or center is
            not a 3-element array.
    """

    if not isinstance(tomo, np.ndarray) or tomo.ndim != 3:
        raise TypeError("tomo must be a 3-D numpy array.")
    if not hasattr(q, "__len__") or len(q) != 4:
        raise ValueError("Quaternion q must have 4 elements.")
    if center is None:
        center = 0.5 * (np.asarray(tomo.shape, dtype=np.float32) - 1)
    else:
        if not isinstance(center, np.ndarray) or len(center) != 3:
            raise ValueError(
                "center must be a numpy array with " "3 elements."
            )

    q /= vector_module(q)
    R = quat_to_mat(q)
    if not active:
        R = R.T

    matrix = R.T
    offset = center - matrix @ center
    return spnd.affine_transform(
        tomo,
        matrix,
        offset=offset,
        order=order,
        mode=mode,
        cval=cval,
        prefilter=prefilter,
        output_shape=tomo.shape,
    )


def vect_rotate(vect, q, active=True):
    """Apply a quaternion-encoded rotation to a 3-D vector.

    Args:
        vect (numpy.ndarray): Input vector of shape (3,).
        q (array-like): Unit quaternion of shape (4,) as
            (w, x, y, z).
        active (bool): If True (default), applies an active
            rotation; otherwise a passive (inverse) rotation.

    Returns:
        numpy.ndarray: The rotated vector of shape (3,).

    Raises:
        TypeError: If vect is not a numpy array with 3 elements.
        ValueError: If q does not have exactly 4 elements.
    """

    if not isinstance(vect, np.ndarray) or len(vect) != 3:
        raise TypeError("vect must be a numpy array with " "3 elements.")
    if not hasattr(q, "__len__") or len(q) != 4:
        raise ValueError("Quaternion q must have 4 elements.")

    q /= vector_module(q)
    R = quat_to_mat(q)
    if not active:
        R = R.T

    return np.matmul(R, vect)


def quat_mult(q1, q2):
    """Multiply two quaternions (Hamilton product q1 * q2).

    Args:
        q1 (array-like): First quaternion of shape (4,) as
            (w, x, y, z).
        q2 (array-like): Second quaternion of shape (4,) as
            (w, x, y, z).

    Returns:
        numpy.ndarray: Resulting quaternion of shape (4,) as
            float64.
    """
    w0, x0, y0, z0 = q1
    w1, x1, y1, z1 = q2
    hold_q = np.array(
        [
            -x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
            x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
            -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
            x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0,
        ],
        dtype=np.float64,
    )
    return hold_q


def quat_two_vectors(a, b):
    """Compute the shortest-arc quaternion rotating vector a onto b.

    Args:
        a (array-like): Origin 3-D vector.
        b (array-like): Destination 3-D vector.

    Returns:
        numpy.ndarray: Unit quaternion of shape (4,) encoding the
            rotation from a to b.
    """
    u, v = a / vector_module(a), b / vector_module(b)
    if vector_module(u - v) < LOW_VALUE:
        ax = ortho_vector(u)
        ax_n = ax / vector_module(ax)
        return np.asarray((0.0, ax_n[0], ax_n[1], ax_n[2]))
    else:
        half = u + v
        half /= vector_module(half)
        ang, ax = np.dot(u, v), np.cross(u, v)
        return np.asarray((ang, ax[0], ax[1], ax[2]))


def ortho_vector(v):
    """Return any vector orthogonal to the input.

    Cross-multiplies v with the most orthogonal standard basis
    vector to avoid numerical degeneracy.

    Args:
        v (array-like): Input 3-D vector.

    Returns:
        numpy.ndarray: A vector perpendicular to v (not
            necessarily unit-length).
    """
    # Handle input is a basis
    if (v[0] < LOW_VALUE) and (v[1] < LOW_VALUE):
        return np.asarray((1.0, 0.0, 0.0))
    elif (v[0] < LOW_VALUE) and (v[2] < LOW_VALUE):
        return np.asarray((1.0, 0.0, 0.0))
    elif (v[1] < LOW_VALUE) and (v[2] < LOW_VALUE):
        return np.asarray((0.0, 1.0, 0.0))
    # Search for the most orthogonal basis vector
    v_abs = np.abs(v)
    if v_abs[0] < v_abs[1]:
        if v_abs[0] < v_abs[2]:
            other = np.asarray((1.0, 0.0, 0.0))
        else:
            other = np.asarray((0.0, 0.0, 1.0))
    else:
        if v_abs[1] < v_abs[2]:
            other = np.asarray((0.0, 1.0, 0.0))
        else:
            other = np.asarray((0.0, 0.0, 1.0))
    return np.cross(v, other)


def vect_to_zmat(v_in, mode="active"):
    """Compute the rotation matrix that maps the Z-axis onto v_in.

    Angles are computed in the extrinsic ZYZ system and converted
    to the Relion intrinsic ZY'Z'' convention.

    Args:
        v_in (array-like): Target 3-D unit vector.
        mode (str): Either 'active' (default) or 'passive';
            passive returns the transpose.

    Returns:
        numpy.ndarray: A (3, 3) rotation matrix.
    """

    n = v_in / vector_module(v_in)

    # Computing angles in Extrinsic ZYZ system
    alpha = np.arccos(n[2])
    beta = np.arctan2(n[1], n[0])

    # Transform to Relion system (intrinsic ZY'Z'' where rho is free)
    rot, tilt, psi = (
        0.0,
        wrap_angle(math.degrees(alpha), deg=True),
        wrap_angle(180.0 - math.degrees(beta), deg=True),
    )

    M = rot_mat_zyz(rot, tilt, psi, deg=True)

    if mode == "passive":
        M = M.T

    return M


def rot_mat_zyz(rot, tilt, psi, deg=True):
    """Build a rotation matrix from ZY'Z'' Euler angles.

    Compatible with Relion and XMIPP conventions. Translated from
    the Relion C++ source (euler.cpp).

    Args:
        rot (float): First Euler angle Z.
        tilt (float): Second Euler angle Y'.
        psi (float): Third Euler angle Z''.
        deg (bool): If True (default), angles are in degrees;
            otherwise in radians.

    Returns:
        numpy.matrix: A (3, 3) rotation matrix (float32).
    """

    # XMIPP doc
    if deg:
        rot, tilt, psi = (
            math.radians(rot),
            math.radians(tilt),
            math.radians(psi),
        )
    mt = np.zeros(shape=(3, 3), dtype=np.float32)
    ca, sa = math.cos(rot), math.sin(rot)
    cb, sb = math.cos(tilt), math.sin(tilt)
    cg, sg = math.cos(psi), math.sin(psi)
    cc, cs = cb * ca, cb * sa
    sc, ss = sb * ca, sb * sa

    # XMIPP doc inverted
    mt[0][0] = cg * cc - sg * sa
    mt[1][0] = cg * cs + sg * ca
    mt[2][0] = -cg * sb
    mt[0][1] = -sg * cc - cg * sa
    mt[1][1] = -sg * cs + cg * ca
    mt[2][1] = sg * sb
    mt[0][2] = sc
    mt[1][2] = ss
    mt[2][2] = cb

    return np.asmatrix(mt)


def rot_to_quat(rot):
    """Convert a rotation matrix to a unit quaternion.

    Args:
        rot (numpy.ndarray | numpy.matrix): A (3, 3) rotation
            matrix.

    Returns:
        numpy.ndarray: Unit quaternion of shape (4,) as
            (w, x, y, z) float64.
    """
    r = spR.from_matrix(rot)
    hold_q = r.as_quat()
    return np.asarray(
        (hold_q[3], hold_q[0], hold_q[1], hold_q[2]), dtype=float
    )


def tomo_shift(tomo, shift):
    """Shift a 3-D tomogram in Fourier space.

    Applies a sub-voxel-accurate rigid translation by multiplying
    the FFT of the volume with a linear phase ramp.

    Args:
        tomo (numpy.ndarray): Input 3-D volume.
        shift (array-like): Shift in voxels for each of the three
            dimensions.

    Returns:
        numpy.ndarray: The shifted volume as a real-valued float
            array.

    Raises:
        TypeError: If tomo is not a 3-D numpy array.
        ValueError: If shift does not have exactly 3 elements.
    """

    if not isinstance(tomo, np.ndarray) or tomo.ndim != 3:
        raise TypeError("tomo must be a 3-D numpy array.")
    if not hasattr(shift, "__len__") or len(shift) != 3:
        raise ValueError("shift must have exactly 3 elements.")
    dx, dy, dz = (
        float(tomo.shape[0]),
        float(tomo.shape[1]),
        float(tomo.shape[2]),
    )
    dx2, dy2, dz2 = (
        math.floor(0.5 * dx),
        math.floor(0.5 * dy),
        math.floor(0.5 * dz),
    )
    if isinstance(shift, np.ndarray):
        delta = np.copy(shift)
    else:
        delta = np.asarray(shift, dtype=np.float32)
    dim = np.asarray((dx, dy, dz), dtype=np.float32)

    x_l, y_l, z_l = -dx2, -dy2, -dz2
    x_h, y_h, z_h = -dx2 + dx, -dy2 + dy, -dz2 + dz
    X, Y, Z = np.meshgrid(
        np.arange(x_l, x_h),
        np.arange(y_l, y_h),
        np.arange(z_l, z_h),
        indexing="ij",
    )

    # Check for trivial dimensions
    ids = np.where(dim <= 1)[0]
    delta[ids] = 0

    # Shift grid in Fourier space
    delta[0], delta[1], delta[2] = delta[0] / dx, delta[1] / dy, delta[2] / dz
    X = np.fft.ifftshift(delta[0] * X + delta[1] * Y + delta[2] * Z)
    del Y, Z

    # Tomogram shifting in Fourier space
    j = 1j
    img = np.fft.fftn(tomo)
    return np.real(np.fft.ifftn(img * np.exp(-2.0 * np.pi * j * X)))


def uniform_sampling_so3(n):
    """Generate n deterministic uniform samples on SO(3).

    Uses the Super-Fibonacci spiral method.

    Reference:
        Alexa, M. (2022). Super-Fibonacci Spirals: Fast,
        Low-Discrepancy Sampling of SO(3). CVPR 2022,
        doi:10.1109/CVPR52688.2022.00811.

    Args:
        n (int): Number of samples to generate.

    Returns:
        numpy.ndarray: Array of shape (n, 4) containing unit
            quaternions.
    """

    Q = np.empty(shape=(n, 4), dtype=float)

    for i in range(n):
        s = i + 0.5
        r = np.sqrt(s / n)
        R = np.sqrt(1.0 - s / n)
        alpha = 2.0 * np.pi * s / PHI
        beta = 2.0 * np.pi * s / PSI
        Q[i, 0] = r * np.sin(alpha)
        Q[i, 1] = r * np.cos(alpha)
        Q[i, 2] = R * np.sin(beta)
        Q[i, 3] = R * np.cos(beta)

    return Q
