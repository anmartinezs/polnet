"""
I/O functions

"""

__author__ = 'Antonio Martinez-Sanchez'

import vtk
import csv
import mrcfile
import numpy as np
import pandas as pd


def load_mrc(fname, mmap=False, no_saxes=True):
    """
    Load an input MRC tomogram as ndarray
    :param fname: the input MRC
    :param mmap: if True (default False) the data are read as a memory map
    :param no_saxes: if True (default) then X and Y axes are swaped to cancel the swaping made by mrcfile package
    :return: a ndarray (or memmap is mmap=True)
    """
    if mmap:
        mrc = mrcfile.mmap(fname, permissive=True, mode='r+')
    else:
        mrc = mrcfile.open(fname, permissive=True, mode='r+')
    if no_saxes:
        return np.swapaxes(mrc.data, 0, 2)
    return mrc.data


def write_mrc(tomo, fname, v_size=1, dtype=None, no_saxes=True):
    """
    Saves a tomo (3D dataset) as MRC file
    :param tomo: tomo to save as ndarray
    :param fname: output file path
    :param v_size: voxel size (default 1)
    :param dtype: data type (default None, then the dtype of tomo is considered)
    :param no_saxes: if True (default) then X and Y axes are swaped to cancel the swaping made by mrcfile package
    :return:
    """
    with mrcfile.new(fname, overwrite=True) as mrc:
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


def load_csv_into_tomo_tables(in_csv_file):
    """
    Load a CSV file as a dictionary of tables, one for each density
    :param in_csv_file: input CSV file path
    :return: a dictionary where each density path is an entry for a table, each table contains all particles of single
             density
    """
    tables_df = pd.read_csv(in_csv_file, delimiter='\t', names=['Density Micrographs', 'PolyData', 'Tomo3D', 'Type',
                                                                 'Label', 'Code', 'Polymer', 'X', 'Y', 'Z',
                                                                 'Q1', 'Q2', 'Q3', 'Q4'], header=0)
    den_tomos = set(tables_df['Tomo3D'].tolist())
    tables_dic = dict().fromkeys(den_tomos)
    for key in tables_dic:
        tables_dic[key] = dict().fromkeys(tables_df.columns.tolist())
        for kkey in tables_dic[key]:
            tables_dic[key][kkey] = list()
    for row in tables_df.iterrows():
        key = row[1]['Tomo3D']
        for item, value in row[1].items():
            tables_dic[key][item].append(value)
    return tables_dic


def write_table(table, out_file):
    """
    Store a table in a CSV file
    :param table: input table dictionary
    :param out_file: path for the output file
    :return:
    """
    with open(out_file, 'w', newline='') as csv_file:
        fieldnames = list(table.keys())
        writer = csv.DictWriter(csv_file, fieldnames=fieldnames, delimiter='\t')
        writer.writeheader()
        for row in range(len(table[fieldnames[0]])):
            dic_row = dict().fromkeys(fieldnames)
            for key in fieldnames:
                dic_row[key] = table[key][row]
            writer.writerow(dic_row)