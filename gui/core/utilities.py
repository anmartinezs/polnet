import wget
import os
from polnet import lio, polymer, utils
import numpy as np
from .pdbtomrc import pdb_2_mrc
import skimage
from scipy.ndimage import rotate, affine_transform
from .tk_utilities import *


def check_path(destination_path, name, ext):
    """
    Check correct path
    :param destination_path: destination path to check
    :param name: if path dont have name 
    :param ext: check correct extension
    """
    if destination_path:   
        if not os.path.basename(destination_path):
            destination_path = os.path.join(destination_path, f"{name}{ext}")
        if not destination_path.lower().endswith(ext):
            if destination_path.lower().endswith("/"):
                destination_path += name + ext 
            else:
                destination_path += ext

    return destination_path
    

def download_pdb(code,destination_path):  
    """
    Donwload pdb from protein data bank
    :param code: protein name 
    :param destination_path: download in destination path
    :return: path to visualize it
    """
    pdb_url = f"https://files.rcsb.org/download/{code}.pdb"
    destination_path = check_path(destination_path, code, ".pdb")

    print(destination_path)
    try:
        if os.path.exists(destination_path):
            os.remove(destination_path)  
        wget.download(pdb_url, destination_path)
        print(f"Archivo PDB {code} descargado y guardado en: {destination_path}")
    except Exception as e:
        print(f"Error al descargar el archivo PDB {code}: {str(e)}")

    return destination_path

                
def create_new_pdb(upload_value):
    """
    Create new file to a file uploaded
    :param upload_value: dictionary with upload widget values
    :return: new path created
    """
    print(upload_value)
    name = upload_value['name']
    content = upload_value['content']
    
    # New path 
    new_path = f"../data/templates/pdbs/{name}"
    
    with open(new_path, 'wb') as new_file:
        new_file.write(content)
    return new_path
    

def convert_to_mrc(file_paths, destination_path, apix, res, offset, het, selected_atoms, model):
    """
    convert pdb to mrc using pdb2mrc
    :param files: files selected
    :param file_destination_widget: FileChooser widget for choosing destination folder
    :param apix_widget: BoundedFloatText widget for Pixels per Ångström
    :param offset_widget: BoundedFloatText widget for Additional Offset
    """
    for path in file_paths:
        p = os.path.basename(path)
        iden = os.path.splitext(p)[0]
        output = destination_path + iden + ".mrc"
        pdb_2_mrc(path, output, apix, res, offset, het, selected_atoms, model)
        window_convert_to_mrc(output)
          

def create_axis_mrc(mb_size, zaxis_rad, vsize, files_path, out_dir, scale_factor):   
    """
    Create an mrc axis for membrane proteins
    :param mb_size: membrane size in amstrongs
    :param zaxis_rad: rad z axis in amstrongs
    :param vsize: voxel size
    :param files_path: files to create axis
    :param out_dir: path where save the files
    :param scale_factor: to resize the tomo box size
    """
    mb_size_vx, zaxis_rad_vx = mb_size / vsize, zaxis_rad / vsize
    for in_mrc in files_path:
        tomo = lio.load_mrc(in_mrc)
        tomo= skimage.transform.rescale(tomo, scale=scale_factor, order=0, anti_aliasing=True)
        # Constructing the reference subvolume
        center = .5 * (np.asarray(tomo.shape, dtype=np.float32) - 1)
        X, Y, Z = np.meshgrid(np.arange(tomo.shape[0]), np.arange(tomo.shape[1]), np.arange(tomo.shape[2]), indexing='ij')
        X, Y, Z = (X - center[0]).astype(np.float32), (Y - center[1]).astype(np.float32), (Z - center[2]).astype(np.float32)
        mask_mb = (Z <= 0) * (Z >= -mb_size_vx)
        zaxis_mask = X*X + Y*Y <= zaxis_rad
        ref_svol = np.zeros(shape=tomo.shape, dtype=np.int16)
        ref_svol[mask_mb] = 1
        ref_svol[zaxis_mask] = 1
    
        # First command: guess minimal dimension
        hold_out_mrc = out_dir + os.path.splitext(os.path.split(in_mrc)[1])[0] + '_mb_ref.mrc'
        lio.write_mrc(ref_svol, hold_out_mrc, v_size=vsize)
        window_convert_to_mrc(hold_out_mrc)



def insert_maxis(protein_array, eje_array, center, path, name):
    """
    Insert tomo in tomo axis
    :param protein_array: input subvolume (or subtomogram)
    :param eje_array: input tomogram that is going to be modified
    :param center: subvolume center point
    :param path: path to save te final tomogram
    :param angles: vector to rotate the protein
    :return: path 
    """
    utils.insert_svol_tomo(protein_array, eje_array, center, merge='sum')
    name = name + "_mbz_align_shitf"
    output = check_path(path, name, ".mrc")
    print("ESte es el path", output)
    lio.write_mrc(eje_array, output, v_size=10)
    return output


def write_mmolecules(mmer_id, path, outpath, mmer_iso, pmer_l, pmer_l_max, pmer_occ, pmer_over_tol, pmer_reverse_normal, is_membrane):
    """
    Create a pns or pms file with mmolecules content
    :parama mmer_id: identifier
    :param path: path to density .mrc
    :param outpath: path to save the file
    :param mmer_iso: isosurface threshold
    :param pmer_l: 
    :param pmer_occ: cluster max length
    :param pmer_over_tol: overlapping tolerance
    :param is_membrane: indicate protein type
    :param proteins_list: add the path if standard protein
    :param mb_protein_list: add the path if membrane protein
    :return: file_path to save it in a array
    """
    # File content
    file_content = f"""MMER_ID = {mmer_id}
MMER_SVOL = {path}
MMER_ISO = {mmer_iso}
PMER_L = {pmer_l}
PMER_L_MAX = {pmer_l_max}
PMER_OCC = {pmer_occ}
PMER_OVER_TOL = {pmer_over_tol}
"""

    # File name
    nombre_archivo = os.path.basename(path)
    iden = os.path.splitext(nombre_archivo)[0]
    partes = iden.split("_")
    primer_segmento = partes[0] if partes else iden
    
    if is_membrane:
        file_content += f"""PMER_REVERSE_NORMALS = {pmer_reverse_normal}"""
        file_extension = ".pms"
        name = "mb_"+ primer_segmento + "_10A"
        file_path = check_path(outpath, name, file_extension) 
    else:
        file_extension = ".pns"
        name = primer_segmento + "_10A" 
        file_path = check_path(outpath, name, file_extension) 
    # Save content
    
    with open(file_path, "w") as file:
        file.write(file_content)

    window_create_mmolecules(file_path)


def write_membrane(mb_type, mb_occ, mb_thick_rg, mb_layer_s_rg, mb_max_ecc, mb_over_tol, mb_min_rad, mb_den_cf_rg, outpath):
    """
    Create a mbs file 
    :parama mb_type: membrane geometetry
    :param mb_occ: ocupancy
    :param mb_thick_rg: membrane tickness or distance between both layers
    :param mb_layer_s_rg: membrane layers tickness
    :param mb_max_ecc: maximun ellipsoid eccentricity
    :param mb_over_tol: overlapping tolerance
    :param mb_min_rad: minimum membrane sphere radius
    :param mb_den_cf_rg: density factor 
    :parama outpath: path to save the file
    """
    # Crear el contenido del archivo de texto
    file_content = f"""MB_TYPE = {mb_type}
MB_OCC = {mb_occ}
MB_THICK_RG = {mb_thick_rg} # A
MB_LAYER_S_RG = {mb_layer_s_rg} # A
MB_MAX_ECC = {mb_max_ecc}
MB_OVER_TOL = {mb_over_tol} # In percentage
MB_MIN_RAD = {mb_min_rad}
MB_DEN_CF_RG = {mb_den_cf_rg}
"""

    # Guardar el contenido en un archivo de texto
    file_path = check_path(outpath, mb_type, ".mbs")
    
    with open(file_path, "w") as file:
        file.write(file_content)

    window_create_membranes(file_path)


def write_helix(hlix_type, hlix_mmer_rad, hlix_pmer_l, hlix_pmer_occ, hlix_min_p_len, hlix_hp_len, hlix_mz_len, hlix_mz_len_f, hlix_over_tol, hlix_min_nmmer, a_bprop, a_max_p_branch, mt_rad, mt_nunits, outpath):
    """
    Create hns file
    :param hlix_type: type of helix.
    :param hlix_mmer_rad: helix mmer radius.
    :param hlix_pmer_l: helix pmer length.
    :param hlix_pmer_occ: helix pmer occupancy.
    :param hlix_min_p_len : helix minimum p length.
    :param hlix_hp_len: helix hp length.
    :param hlix_mz_len: helix mz length.
    :param hlix_mz_len_f: helix mz length factor.
    :param hlix_over_tol: helix over-tolerance.
    :param hlix_min_nmmer: helix minimum nmmer.
    :param a_bprop: helix a bprop.
    :param a_max_p_branch: helix a max p branch.
    :param mt_rad: microtubule radius.
    :param mt_nunits: Microtubule number of units.
    :param outpath: path to save file
    """
    # File content
    file_content = f'''HLIX_TYPE = {hlix_type}
HLIX_MMER_RAD = {hlix_mmer_rad}
HLIX_PMER_L = {hlix_pmer_l}
HLIX_PMER_OCC = {hlix_pmer_occ}
HLIX_MIN_P_LEN = {hlix_min_p_len}
HLIX_HP_LEN = {hlix_hp_len}
HLIX_MZ_LEN = {hlix_mz_len}
HLIX_MZ_LEN_F = {hlix_mz_len_f}
HLIX_OVER_TOL = {hlix_over_tol}
HLIX_MIN_NMMER = {hlix_min_nmmer}
'''
    if "actin" == hlix_type:
        file_content += f'''
A_BPROP  = {a_bprop}
A_MAX_P_BRANCH  = {a_max_p_branch}
'''
    else:
        file_content+= f'''
MT_RAD = {mt_rad}
MT_NUNITS = {mt_nunits}
'''
    file_path = check_path(outpath, hlix_type, ".hns")
    
    # Escribir el contenido en un archivo de texto
    with open(file_path, 'w') as file:
        file.write(file_content)

    window_create_helix(file_path)


def create_actin_poly_data(hlix_mmer_rad, v_size):
    """
    Create acti poly data
    :param hlix_mmer_dad: monomer radius
    :param v_size: voxel size
    :return: vtk poly data
    """
    actin_filaments = polymer.FiberUnitSDimer(sph_rad=hlix_mmer_rad, v_size=v_size)  
    vtk_poly_data = actin_filaments.get_vtp()
    return vtk_poly_data


def create_mt_poly_data(hlix_mmer_rad, mt_rad, mt_units, v_size):
    """
    Create a mt poly data
    :param hlix_mmer_rad: monomer radius
    :param mt_rad: mt radius
    :param mt_units: mt ring number of monomers
    :param v_size: voxel size
    :return: vtk poly data
    """
    microtubes = polymer.MTUnit(sph_rad=hlix_mmer_rad, mt_rad= mt_rad, n_units=mt_units, v_size=v_size)  
    vtk_poly_data = microtubes.get_vtp()
    return vtk_poly_data


def create_poly_data(file_path, v_size):
    """
    Create the poly data
    :param file_path: path to the file
    :param v_size: pixel size
    :return: vtk poly data 
    """
    if not os.path.exists(file_path):
        raise FileNotFoundError(f"File not found: {file_path}")

    with open(file_path, 'r') as file:
        file_content = file.read()

    lines = file_content.split('\n')

    lines = [line for line in lines if line.strip()]

    print("Number of lines:", len(lines))  
    hlix_type = lines[0].split('=')[1].strip()
    hlix_mmer_rad = float(lines[1].split('=')[1].strip())
    print(hlix_mmer_rad)
    # Add more parameter parsing as needed
    if hlix_type == "actin":
        vtk_poly_data = create_actin_poly_data(hlix_mmer_rad, v_size)
        return vtk_poly_data
    else:
        mt_rad = float(lines[10].split('=')[1].strip())
        mt_units = int(lines[11].split('=')[1].strip())
        print(mt_rad, " " , mt_units)
        vtk_poly_data = create_mt_poly_data(hlix_mmer_rad, mt_rad, mt_units, v_size)
        print(vtk_poly_data)
        return vtk_poly_data


def files_selected(files_path, files_label):
    """
    Change file upload label
    :files_path: list with files
    :files_label: label to update
    """
    num_files = len(files_path) 
    if num_files < 3:
        text = ", ".join(files_path) if num_files > 0 else "No files selected"
    else:
        text = ", ".join(files_path[:2]) + ", ..."
    
    files_label.value = text


def check_dir(dir_path, path):
    """
    Check correct out dir
    :param dir_path: path selected by user
    :param path: default path 
    :return: correct path
    """
    if dir_path == None:
        dir_path = path
    return dir_path


def check_prop_list(check_box, prop_list):
    """
    Check widget prop list
    :param check_box: boolean to see it is active
    :param prop_list: to get the value
    :return: prop list value
    """
    if check_box:
        return prop_list
    else:
        return None
        
    
    
