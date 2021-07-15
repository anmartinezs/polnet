"""
Functionality for Affine Transformations
"""

__author__ = 'Antonio Martinez-Sanchez'

import vtk
import math
import random
import numpy as np
from scipy import ndimage as spnd

# CONSTANTS

PI_2 = 2 * np.pi

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
    x = math.sin(theta_2) * sigma_1
    y = math.cos(theta_2) * sigma_1
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
    if deg:
        ang_rad = (2. * math.atan2(norm, qw))
        ang = (180./np.pi) * ang_rad
        return ang, x / math.sin(.5 * ang_rad)
    else:
        ang = 2. *  math.atan2(norm, qw)
        return ang, x / math.sin(.5 * ang)


def rot_vect_quat(v, q):
    """
    Rotate a vector by a quaternion using Rodrigues formula
    :param v: input vector
    :param q: input quaternion (angle, axis)
    :return: output rotated vector
    """

    # Make sure axis is a unit vector
    k = np.asarray((q[1], q[2], q[3]))
    mod_k = math.sqrt((k * k).sum())
    assert mod_k > 0
    k /= mod_k
    vc = np.asarray(v)

    # Rodrigues formula
    cos_ang, sin_ang = math.cos(q[0]), math.sin(q[0])
    return vc * cos_ang + np.cross(k, vc) * sin_ang + k * np.dot(k, vc) * (1.-cos_ang)


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
    assert len(q) == 4
    if center is None:
        # center = np.round(.5 * np.asarray(tomo.shape, dtype=np.float32))
        center = .5 * np.asarray(tomo.shape, dtype=np.float32)
    else:
        assert isinstance(tomo, np.ndarray) and (len(center.shape) == 3)

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


