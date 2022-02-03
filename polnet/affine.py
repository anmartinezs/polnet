"""
Functionality for Affine Transformations
"""

__author__ = 'Antonio Martinez-Sanchez'

import vtk
import math
import random
import numpy as np
from scipy import ndimage as spnd
from polnet.utils import wrap_angle
from scipy.spatial.transform import Rotation as spR

# CONSTANTS

PI_2 = 2 * np.pi
LOW_VALUE = 0.000001

# FUNCTIONS


def vector_module(v):
    """
    Computes the module of a vector (|v|)
    :param v: input vector (one dimensional numpy array)
    :return: the computed module
    """
    return math.sqrt((v*v).sum())


def poly_rotate_wxyz(in_vtp, w, x, y, z):
    """
    Applies rigid rotation to a vtkPolyData from angle-axis input form
    :param in_vtp: input vtkPolyData
    :param w: angle
    :param x: axis - x
    :param y: axis - y
    :param z: axis - z
    :return:
    """
    rot_tf = vtk.vtkTransform()
    rot_tf.RotateWXYZ(w, x, y, z)
    tf_rot = vtk.vtkTransformPolyDataFilter()
    tf_rot.SetInputData(in_vtp)
    tf_rot.SetTransform(rot_tf)
    tf_rot.Update()
    return tf_rot.GetOutput()


def poly_translate(in_vtp, t):
    """
    Applies rigid rotation to a vtkPolyData
    :param in_vtp: input vtkPolyData
    :param t: translation vector
    :return: the transformed vtkPolyData
    """

    # Input parsing
    assert hasattr(t, '__len__') and (len(t) == 3)

    # Translation
    box_tr = vtk.vtkTransform()
    box_tr.Translate(t[0], t[1], t[2])
    tr_box = vtk.vtkTransformPolyDataFilter()
    tr_box.SetInputData(in_vtp)
    tr_box.SetTransform(box_tr)
    tr_box.Update()
    return tr_box.GetOutput()


def poly_scale(in_vtp, s):
    """
    Applies scaling transformation to a vtkPolyData
    :param in_vtp: input vtkPolyData
    :param s: scaling factor
    :return: the transformed vtkPolyData
    """
    # Translation
    box_sc = vtk.vtkTransform()
    box_sc.Scale(s, s, s)
    tr_box = vtk.vtkTransformPolyDataFilter()
    tr_box.SetInputData(in_vtp)
    tr_box.SetTransform(box_sc)
    tr_box.Update()
    return tr_box.GetOutput()


def gen_rand_unit_quaternion():
    """
    Generates a random sample for a unit quaternion
    KUFFNER, James J. Effective sampling and distance metrics for 3D rigid body path planning.
    In IEEE International Conference on Robotics and Automation, 2004. Proceedings. ICRA'04. 2004. IEEE, 2004. p. 3993-3998.
    :return: a unit quaternion in the shape Q=(w,x,y,z)=(cos(θ/2),vx*sin(θ/2),vy*sin(θ/2),vz*sin(θ/2)),
             where (vx,vy,vz) is the axis and \theta is the angle
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
    """
    Given a quaternion, q = s + X with s=qw and X=[qx, qy, qz], returns its angle axis representation
    :param qw: quaternion w value
    :param qx: quaternion x value
    :param qy: quaternion y value
    :param qz: quaternion z value
    :param deg: if True (default) the output is given in degrees, otherwise in radians
    :return: an angle in degrees and axis (3-array)
    """
    x = np.asarray((qx, qy, qz))
    norm = vector_module(x)
    if norm <= 0:
        return 0, np.asarray((0., 0., 0.))
    if deg:
        ang_rad = (2. * math.atan2(norm, qw))
        ang = math.degrees(ang_rad)
        return ang, x / math.sin(.5 * ang_rad)
    else:
        ang = 2. *  math.atan2(norm, qw)
        return ang, x / math.sin(.5 * ang)


def angle_axis_to_quat(ang, x, y, z, deg=True):
    """
    Given a angle and axis [x, y, z] computes its corresponding unit quaternion
    :param ang: angle
    :param x: axis x value
    :param y: axis y value
    :param z: axis z value
    :param deg: if True (default) the input angle is given in degrees, otherwise in radians
    :return: a unit quaternion (4-array)
    """
    ax = np.asarray((x, y, z), dtype=float)
    ax /= vector_module(ax)
    hold_ang = .5 * wrap_angle(ang, deg=deg)
    if deg:
        hold_ang = math.radians(hold_ang)
    ca, sa = math.cos(hold_ang), math.sin(hold_ang)
    return np.asarray((ca, ax[0]*sa, ax[1]*sa, ax[2]*sa), dtype=float)


def rot_vect_quat(v, q):
    """
    Rotate a vector by a quaternion using Rodrigues formula
    :param v: input vector
    :param q: input quaternion (angle, axis)
    :return: output rotated vector
    """

    # Get axis angle representation
    ang, ax = quat_to_angle_axis(q[0], q[1], q[2], q[3])

    # Make sure axis is a unit vector
    k = np.asarray((ax[0], ax[1], ax[2]))
    mod_k = math.sqrt((k * k).sum())
    if mod_k > 0:
        k /= mod_k
        vc = np.asarray(v)

        # Rodrigues formula
        cos_ang, sin_ang = math.cos(ang), math.sin(ang)
        return vc * cos_ang + np.cross(k, vc) * sin_ang + k * np.dot(k, vc) * (1.-cos_ang)
    else:
        return v


def quat_to_mat(q):
    """
    Covert a quaternion into a full three-dimensional rotation matrix
    :param q: A 4 element array representing the quaternion (q0,q1,q2,q3)
    :return: A 3x3 element matrix representing the full 3D rotation matrix.
             This rotation matrix converts a point in the local reference
             frame to a point in the global reference frame.
    """
    # Extract the values from Q
    q0, q1, q2, q3 = q[0], q[1], q[2], q[3]

    # First row of the rotation matrix
    r00 = 2 * (q0 * q0 + q1 * q1) - 1
    r01 = 2 * (q1 * q2 - q0 * q3)
    r02 = 2 * (q1 * q3 + q0 * q2)

    # Second row of the rotation matrix
    r10 = 2 * (q1 * q2 + q0 * q3)
    r11 = 2 * (q0 * q0 + q2 * q2) - 1
    r12 = 2 * (q2 * q3 - q0 * q1)

    # Third row of the rotation matrix
    r20 = 2 * (q1 * q3 - q0 * q2)
    r21 = 2 * (q2 * q3 + q0 * q1)
    r22 = 2 * (q0 * q0 + q3 * q3) - 1

    # 3x3 rotation matrix
    return np.array([[r00, r01, r02], [r10, r11, r12], [r20, r21, r22]])


def tomo_rotate(tomo, q, center=None, active=True, order=3, mode='constant', cval=0.0, prefilter=True):
    """
    Applies the rotation defined in a quaternion to a tomogram
    :param tomo: input tomogram as 3D numpy array
    :param q: quaternion encoding the rotation
    :param center: center for rotation, if None (default) the tomogram center is used
    :param active: rotation mode, if True (default) active otherwise passive
    :param order: see doc for scipy.ndimage.interpolation.map_coordinate
    :param mode: see doc for scipy.ndimage.interpolation.map_coordinate
    :param cval: see doc for scipy.ndimage.interpolation.map_coordinate
    :param prefilter: see doc for scipy.ndimage.interpolation.map_coordinate
    :return: the rotated tomogram
    """

    # Input parsing
    assert isinstance(tomo, np.ndarray) and (len(tomo.shape) == 3)
    assert hasattr(q, '__len__') and (len(q) == 4)
    if center is None:
        # center = np.round(.5 * np.asarray(tomo.shape, dtype=np.float32))
        center = .5 * np.asarray(tomo.shape, dtype=np.float32)
    else:
        assert isinstance(tomo, np.ndarray) and (len(center) == 3)

    # Getting rotation matrix
    q /= vector_module(q)
    R = quat_to_mat(q)
    if not active:
        R = R.T

        # Compute grid
    X, Y, Z = np.meshgrid(np.arange(tomo.shape[0]), np.arange(tomo.shape[1]), np.arange(tomo.shape[2]),
                          indexing='ij')

    # From indices to coordinates
    X, Y, Z = (X - center[0]).astype(np.float32), (Y - center[1]).astype(np.float32), (Z - center[2]).astype(np.float32)

    # Grid rotation
    Xr = X * R[0, 0] + Y * R[1, 0] + Z * R[2, 0]
    Yr = X * R[0, 1] + Y * R[1, 1] + Z * R[2, 1]
    Zr = X * R[0, 2] + Y * R[1, 2] + Z * R[2, 2]

    # From coordinates to indices
    X, Y, Z = Xr + center[0], Yr + center[1], Zr + center[2]

    # Re-mapping (interpolation)
    ts = tomo.size
    inds = np.zeros(shape=(3, ts), dtype=np.float32)
    inds[0, :], inds[1, :], inds[2, :] = X.reshape(ts), Y.reshape(ts), Z.reshape(ts)
    tomo_r = spnd.interpolation.map_coordinates(tomo, inds, order=order, mode=mode, cval=cval, prefilter=prefilter)

    return tomo_r.reshape(tomo.shape)


def vect_rotate(vect, q, active=True):
    """
    Applies the rotation defined in a quaternion to a vector
    :param tomo: input vector as 1D numpy array with length 3
    :param q: quaternion encoding the rotation
    :param active: rotation mode, if True (default) active otherwise passive
    :return: the rotated vector as 1D numpy array
    """

    # Input parsing
    assert isinstance(vect, np.ndarray) and (len(vect) == 3)
    assert hasattr(q, '__len__') and (len(q) == 4)

    # Getting rotation matrix
    q /= vector_module(q)
    R = quat_to_mat(q)
    if not active:
        R = R.T

    # Applying rotation matrix
    return np.matmul(R, vect)


def quat_mult(q1, q2):
    """
    Multiply two quaternions
    :param q1: first input quaternion
    :param q2: second input quaternion
    """
    w0, x0, y0, z0 = q1
    w1, x1, y1, z1 = q2
    hold_q = np.array([-x1 * x0 - y1 * y0 - z1 * z0 + w1 * w0,
                       x1 * w0 + y1 * z0 - z1 * y0 + w1 * x0,
                       -x1 * z0 + y1 * w0 + z1 * x0 + w1 * y0,
                       x1 * y0 - y1 * x0 + z1 * w0 + w1 * z0], dtype=np.float64)
    return hold_q


def quat_two_vectors(a, b):
    """
    Computes the quaternion to rotate from one vector to another
    :param a: origin vector
    :param b: destination vector
    :return: the quaternion thar rotates vector a to be aligned with b
    """
    u, v = a / vector_module(a), b / vector_module(b)
    if vector_module(u - v) < LOW_VALUE:
        ax = ortho_vector(u)
        return np.asarray(0, ax / vector_module(ax))
    else:
        half = u + v
        half /= vector_module(half)
        ang, ax = np.dot(u, v), np.cross(u, v)
        return np.asarray((ang, ax[0], ax[1], ax[2]))


def ortho_vector(v):
    """
    Computes any orthogonal vector to an input. This implementation uses the cross product with the most orthogonal
    basis vector.
    :param v: input vector
    :return: output orthogonal vector
    """
    # Handle input is a basis
    if (v[0] < LOW_VALUE) and (v[1] < LOW_VALUE):
        return np.asarray((1., 0., 0.))
    elif (v[0] < LOW_VALUE) and (v[2] < LOW_VALUE):
        return np.asarray((1., 0., 0.))
    elif (v[1] < LOW_VALUE) and (v[2] < LOW_VALUE):
        return np.asarray((0., 1., 0.))
    # Search for the most orthogonal basis vector
    v_abs = np.abs(v)
    if v_abs[0] < v_abs[1]:
        if v_abs[0] < v_abs[2]:
            other = np.asarray((1., 0., 0.))
        else:
            other = np.asarray((0., 0., 1.))
    else:
        if v_abs[1] < v_abs[2]:
            other = np.asarray((0., 1., 0.))
        else:
            other = np.asarray((0., 0., 1.))
    return np.cross(v, other)


def vect_to_zmat(v_in, mode='active'):
    """
    Computes the matrix to rotate unit Z-axis vector to a given vector
    :param v_in: input vector
    :param mode: either 'active' (default) or 'pasive'
    :returns: a rotation matrix as numpy ndarray (shape=3x3)
    """

    # Normalization
    n = v_in / vector_module(v_in)


    # Computing angles in Extrinsic ZYZ system
    alpha = np.arccos(n[2])
    beta = np.arctan2(n[1], n[0])

    # Transform to Relion system (intrinsic ZY'Z'' where rho is free)
    rot, tilt, psi = 0., wrap_angle(math.degrees(alpha), deg=True), \
                     wrap_angle(180.-math.degrees(beta), deg=True)

    # Matrix computation
    M = rot_mat_zyz(rot, tilt, psi, deg=True)

    # By default is active, invert if passive
    if mode == 'passive':
        M = M.T

    return M


def rot_mat_zyz(rot, tilt, psi, deg=True):
    """
    Creates 3D rotation matrix according to ZY'Z'' Euler angles convention (Relion and Xmipp compatible)
    This is a direct translation from code https://github.com/jjcorreao/relion/blob/master/relion-1.3/src/euler.cpp
    :param rot: first input Euler angle (Z)
    :param tilt: second input Euler angle (X')
    :paramn psi: third input Euler angle (Z'')
    :param deg: if True (default) angles are computed in degrees, otherwise radians
    :return:  a rotation matrix as numpy ndarray (shape=3x3)
    """

    # XMIPP doc
    if deg:
        rot, tilt, psi = math.radians(rot), math.radians(tilt), math.radians(psi)
    mt = np.zeros(shape=(3, 3), dtype=np.float32)
    ca, sa = math.cos(rot), math.sin(rot)
    cb, sb = math.cos(tilt), math.sin(tilt)
    cg, sg = math.cos(psi), math.sin(psi)
    cc, cs = cb*ca, cb*sa
    sc, ss = sb*ca, sb*sa

    # XMIPP doc inverted
    mt[0][0] = cg*cc - sg*sa
    mt[1][0] = cg*cs + sg*ca
    mt[2][0] = -cg*sb
    mt[0][1] = -sg*cc - cg*sa
    mt[1][1] = -sg*cs + cg*ca
    mt[2][1] = sg*sb
    mt[0][2] = sc
    mt[1][2] = ss
    mt[2][2] = cb

    return np.mat(mt)


def rot_to_quat(rot):
    """
    Computes quaternion from an input rotation matrix
    :param rot: rotation numpy matrix
    :return: a quaternion as 4 dirmension array (real part, 3-tuple imaginary part)
    """
    r = spR.from_matrix(rot)
    hold_q = r.as_quat()
    return np.asarray((hold_q[3], hold_q[0], hold_q[1], hold_q[2]), dtype=float)


def tomo_shift(tomo, shift):
    """
    Tomogram shift in Fourier space
    :param tomo: the tomo numpy.ndarray (it must be 3D) to shift
    :param shift: 3-tuple with the shifting for every dimension
    """

    # Input parsing
    assert isinstance(tomo, np.ndarray) and (len(tomo.shape) == 3)
    assert hasattr(shift, '__len__') and (len(shift) == 3)
    dx, dy, dz = float(tomo.shape[0]), float(tomo.shape[1]), float(tomo.shape[2])
    dx2, dy2, dz2 = math.floor(.5*dx), math.floor(.5*dy), math.floor(.5*dz)
    if isinstance(shift, np.ndarray):
        delta = np.copy(shift)
    else:
        delta = np.asarray(shift, dtype=np.float32)
    dim = np.asarray((dx, dy, dz), dtype=np.float32)

    # Generating the grid
    x_l, y_l, z_l = -dx2, -dy2, -dz2
    x_h, y_h, z_h = -dx2+dx, -dy2+dy, -dz2+dz
    X, Y, Z = np.meshgrid(np.arange(x_l, x_h), np.arange(y_l, y_h), np.arange(z_l, z_h), indexing='xy')

    # Check for trivial dimensions
    ids = np.where(dim <= 1)[0]
    delta[ids] = 0

    # Shift grid in Fourier space
    delta[0], delta[1], delta[2] = delta[0]/dx, delta[1]/dy, delta[2]/dz
    X = np.fft.ifftshift(delta[0]*X + delta[1]*Y + delta[2]*Z)
    del Y, Z

    # Tomogram shifting in Fourier space
    j = np.complex(0, 1)
    img = np.fft.fftn(tomo)
    return np.real(np.fft.ifftn(img * np.exp(-2.*np.pi*j*X)))