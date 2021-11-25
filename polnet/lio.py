"""
I/O functions

"""

__author__ = 'Antonio Martinez-Sanchez'

import vtk
import mrcfile
import numpy as np


def load_mrc(fname, mmap=False):
    """
    Load an input MRC tomogram as ndarray
    :param fname: the input MRC
    :param mmap: if True (default False) the data are read as a memory map
    :return: a ndarray (or memmap is mmap=True)
    """
    if mmap:
        mrc = mrcfile.mmap(fname, permissive=True)
    else:
        mrc = mrcfile.open(fname, permissive=True)
        # return np.swapaxes(mrc.data, 0, 2)
    return mrc.data


def write_mrc(tomo, fname, v_size=1, dtype=None):
    """
    Saves a tomo (3D dataset) as MRC file
    :param tomo: tomo to save as ndarray
    :param fname: output file path
    :param v_size: voxel size (default 1)
    :param dtype: data type (default None, then the dtype of tomo is considered)
    :return:
    """
    with mrcfile.new(fname, overwrite=True) as mrc:
        if dtype is None:
            # mrc.set_data(np.swapaxes(tomo, 0, 2))
            mrc.set_data(tomo)
        else:
            # mrc.set_data(np.swapaxes(tomo, 0, 2).astype(dtype))
            mrc.set_data(tomo.astype(dtype))
        mrc.voxel_size.flags.writeable = True
        mrc.voxel_size = (v_size, v_size, v_size)
        mrc.set_volume()
        # mrc.header.ispg = 401


def save_vtp(poly, fname):
    """
    Store data vtkPolyData as a .vtp file
    :param poly: input vtkPolyData to store
    :param fname: output path file
    :return:
    """

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(fname)
    writer.SetInputData(poly)
    if writer.Write() != 1:
        raise IOError


def save_vti(image, fname):
    """
    Store data vtkPolyData as a .vti file
    :param image: input image as numpy array
    :param fname: output path file
    :return:
    """

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(fname)
    writer.SetInputData(image)
    if writer.Write() != 1:
        raise IOError


def load_poly(fname):
    """
    Load data vtkPolyData object from a file
    :param fname: input .vtp file
    :return: the vtkPolyData object loaded
    """

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(fname)
    reader.Update()

    return reader.GetOutput()
