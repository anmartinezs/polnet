"""I/O helper functions for MRC, VTP, and CSV formats.

Wraps mrcfile, VTK, and pandas to provide a unified read/write
interface for volumetric and mesh data used by Polnet.

:author: Antonio Martinez-Sanchez
:maintainer: Juan Diego Gallego Nicolás
"""

import csv

import mrcfile
import numpy as np
import pandas as pd
import vtk
from vtkmodules.util.numpy_support import numpy_to_vtk


def load_mrc(fname, mmap=False, no_saxes=True):
    """Load an MRC tomogram file into a numpy array.

    Args:
        fname (str | Path): Path to the input MRC file.
        mmap (bool): If True, data are read as a memory-mapped
            array (default False).
        no_saxes (bool): If True (default), swaps axes 0 and 2
            to cancel the X/Y swap introduced by mrcfile.

    Returns:
        numpy.ndarray: The tomogram data; a memmap if mmap is
            True.
    """
    if mmap:
        with mrcfile.mmap(str(fname), permissive=True, mode="r") as mrc:
            data = mrc.data.copy()
    else:
        with mrcfile.open(str(fname), permissive=True, mode="r") as mrc:
            data = mrc.data.copy()
    if no_saxes:
        return np.swapaxes(data, 0, 2)
    return data


def write_mrc(tomo, fname, v_size=1, dtype=None, no_saxes=True):
    """Save a 3-D array as an MRC file.

    Args:
        tomo (numpy.ndarray): The 3-D tomogram to save.
        fname (str | Path): Output file path.
        v_size (float): Isotropic voxel size in Angstroms
            (default 1).
        dtype: Data type to cast to before saving; None (default)
            preserves the input dtype.
        no_saxes (bool): If True (default), swaps axes 0 and 2
            before writing to cancel the X/Y swap.
    """
    with mrcfile.new(str(fname), overwrite=True) as mrc:
        if dtype is None:
            if no_saxes:
                mrc.set_data(np.swapaxes(tomo, 0, 2))
            else:
                mrc.set_data(tomo)
        else:
            if no_saxes:
                mrc.set_data(np.swapaxes(tomo, 0, 2).astype(dtype))
            else:
                mrc.set_data(tomo.astype(dtype))
        mrc.voxel_size.flags.writeable = True
        mrc.voxel_size = (v_size, v_size, v_size)
        mrc.set_volume()


def read_mrc_v_size(fname):
    """Read the voxel size from the header of an MRC file.

    Args:
        fname (str | Path): Path to the MRC file.

    Returns:
        tuple[float, float, float]: Voxel size in Angstroms as
            (x, y, z).
    """
    with mrcfile.mmap(str(fname)) as mrc:
        return (mrc.voxel_size["x"], mrc.voxel_size["y"], mrc.voxel_size["z"])


def save_vtp(poly, fname):
    """Store a vtkPolyData object to a .vtp file.

    Args:
        poly (vtk.vtkPolyData): The polygon dataset to store.
        fname (str | Path): Output file path.

    Raises:
        IOError: If the VTK writer reports a failure.
    """

    writer = vtk.vtkXMLPolyDataWriter()
    writer.SetFileName(str(fname))
    writer.SetInputData(poly)
    if writer.Write() != 1:
        raise IOError


def save_vti(image, fname):
    """Store a vtkImageData object to a .vti file.

    Args:
        image (vtk.vtkImageData): The image dataset to store.
        fname (str | Path): Output file path.

    Raises:
        IOError: If the VTK writer reports a failure.
    """

    writer = vtk.vtkXMLImageDataWriter()
    writer.SetFileName(str(fname))
    writer.SetInputData(image)
    if writer.Write() != 1:
        raise IOError


def load_poly(fname):
    """Load a vtkPolyData object from a .vtp file.

    Args:
        fname (str | Path): Path to the input .vtp file.

    Returns:
        vtk.vtkPolyData: The loaded polygon dataset.
    """

    reader = vtk.vtkXMLPolyDataReader()
    reader.SetFileName(str(fname))
    reader.Update()

    return reader.GetOutput()


def load_csv_into_tomo_tables(in_csv_file):
    """Load a per-particle CSV file and organise rows by tomogram.

    Args:
        in_csv_file (str | Path): Path to the tab-separated CSV
            with columns: Density Micrographs, PolyData, Tomo3D,
            Type, Label, Code, Polymer, X, Y, Z, Q1, Q2, Q3, Q4.

    Returns:
        dict: Nested dictionary keyed by Tomo3D path.  Each value
            is a dict mapping column names to lists of row values.
    """
    tables_df = pd.read_csv(
        in_csv_file,
        delimiter="\t",
        names=[
            "Density Micrographs",
            "PolyData",
            "Tomo3D",
            "Type",
            "Label",
            "Code",
            "Polymer",
            "X",
            "Y",
            "Z",
            "Q1",
            "Q2",
            "Q3",
            "Q4",
        ],
        header=0,
    )
    den_tomos = set(tables_df["Tomo3D"].tolist())
    tables_dic = dict().fromkeys(den_tomos)
    for key in tables_dic:
        tables_dic[key] = dict().fromkeys(tables_df.columns.tolist())
        for kkey in tables_dic[key]:
            tables_dic[key][kkey] = list()
    for row in tables_df.iterrows():
        key = row[1]["Tomo3D"]
        for item, value in row[1].items():
            tables_dic[key][item].append(value)
    return tables_dic


def write_table(table, out_file):
    """Write a column-oriented table dictionary to a TSV file.

    Args:
        table (dict): Mapping from column names to lists of row
            values.
        out_file (str | Path): Output file path.
    """
    with open(out_file, "w", newline="", encoding="utf-8") as csv_file:
        fieldnames = list(table.keys())
        writer = csv.DictWriter(
            csv_file, fieldnames=fieldnames, delimiter="\t"
        )
        writer.writeheader()
        for row in range(len(table[fieldnames[0]])):
            dic_row = dict().fromkeys(fieldnames)
            for key in fieldnames:
                dic_row[key] = table[key][row]
            writer.writerow(dic_row)


def numpy_to_vti(tomo: np.ndarray, dtype=vtk.VTK_FLOAT) -> vtk.vtkImageData:
    """Convert a 3-D numpy array to a vtkImageData object.

    Args:
        tomo (numpy.ndarray): Input 3-D array.
        dtype (int): VTK scalar type constant
            (default vtk.VTK_FLOAT).

    Returns:
        vtk.vtkImageData: The resulting VTK image object.

    Raises:
        TypeError: If tomo is not a 3-D array.
    """
    if tomo.ndim != 3:
        raise TypeError("tomo must be a 3-D numpy array.")

    vtk_data = numpy_to_vtk(
        num_array=tomo.flatten(), deep=True, array_type=dtype
    )
    img = vtk.vtkImageData()
    img.GetPointData().SetScalars(vtk_data)
    img.SetDimensions(tomo.shape[0], tomo.shape[1], tomo.shape[2])
    return img
