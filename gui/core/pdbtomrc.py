import numpy as np
import math

import scipy.ndimage

from polnet import lio, utils, affine
import skimage

# By default, sigma_factor = 1/(π * (2/log2)½) ≈ 0.187 which makes the Fourier transform(FT) of the distribution fall to
# fall half maximum value at wavenumber 1 / resolution. Resolution here is 1A
# Other plausible chices are 0.225, 0.356 and 0.425 (see https://www.cgl.ucsf.edu/chimera/docs/UsersGuide/midas/molmap.html)
SIGMA_FACTOR = 0.187

def add_cloud_gauss(tomo, coords, g_std, value):
    """
    Add a translated Gaussian at the specified coordinates
    :param tomo: input tomogram where the Gaussian are going to be added
    :param coords: list of coordinates for the Gaussians
    :param g_std: Gaussian standard deviation
    :param value: Atomic number as mult factor
    :return: a np.ndarray with the added Guassian density
    """

    import time

    hold_time = time.time()
    # Mark Gaussian centers
    hold_tomo = np.zeros(shape=tomo.shape, dtype=np.float32)
    for coord in coords:
        x, y, z = int(round(coord[0])), int(round(coord[1])), int(round(coord[2]))
        # Detect limits
        if 0 <= x < tomo.shape[0] and 0 <= y < tomo.shape[1] and 0 <= z < tomo.shape[2]:
            hold_tomo[x, y, z] = 1
    print('\t\t-Time placing atoms:', time.time() - hold_time, 'secs')

    hold_time = time.time()
    # Building Guassian model
    nx, ny, nz = (tomo.shape[0] - 1) * .5, (tomo.shape[1] - 1) * .5, (tomo.shape[2] - 1) * .5
    if (nx % 1) == 0:
        arr_x = np.concatenate((np.arange(-nx, 0, 1), np.arange(0, nx + 1, 1)))
    else:
        if nx < 1:
            arr_x = np.arange(0, 1)
        else:
            nx = math.ceil(nx)
            arr_x = np.concatenate((np.arange(-nx, 0, 1), np.arange(0, nx, 1)))
    if (ny % 1) == 0:
        arr_y = np.concatenate((np.arange(-ny, 0, 1), np.arange(0, ny + 1, 1)))
    else:
        if ny < 1:
            arr_y = np.arange(0, 1)
        else:
            ny = math.ceil(ny)
            arr_y = np.concatenate((np.arange(-ny, 0, 1), np.arange(0, ny, 1)))
    if (nz % 1) == 0:
        arr_z = np.concatenate((np.arange(-nz, 0, 1), np.arange(0, nz + 1, 1)))
    else:
        if nz < 1:
            arr_z = np.arange(0, 1)
        else:
            nz = math.ceil(nz)
            arr_z = np.concatenate((np.arange(-nz, 0, 1), np.arange(0, nz, 1)))
    [X, Y, Z] = np.meshgrid(arr_x, arr_y, arr_z, indexing='ij')
    X = X.astype(np.float32, copy=False)
    Y = Y.astype(np.float32, copy=False)
    Z = Z.astype(np.float32, copy=False)
    R = np.sqrt(X * X + Y * Y + Z * Z)
    del X
    del Y
    del Z
    gauss = (1/(g_std*g_std*g_std*math.sqrt(2*np.pi))) * np.exp(-R / (2. * g_std * g_std))
    print('\t\t-Time building Gaussian model:', time.time() - hold_time, 'secs')

    hold_time = time.time()
    # Convolution
    tomo_conv = np.fft.fftshift(np.real(np.fft.ifftn(np.fft.fftn(hold_tomo) * np.fft.fftn(gauss))))
    # import scipy
    # tomo_conv = scipy.fft.fftshift(scipy.real(scipy.fft.ifftn(scipy.fft.fftn(hold_tomo) * scipy.fft.fftn(gauss))))

    # Add influence
    tomo_conv *= value
    print('\t\t-Time convolution:', time.time() - hold_time, 'secs')
    # Adding the Gaussian to the input tomogram
    return tomo_conv


def pdb_2_mrc(file_name, output_file, apix, res, offset, het, selected_atoms, model):
    """
    Converts a PDB file to an MRC file.

    :param file_name: Input PDB file name.
    :param output_file: Output MRC file name.
    :param apix: Voxel size per Ångström, scale of the tomogram resolution.
    :param res: Tomogram resolution in Ångströms.
    :param offset: Additional offset for the tomogram size.
    :param het: Flag to include HETATM atoms (True) or not (False).
    :param selected_atoms: List of protein to include in the tomogram (None to include all).
    :param model: Model number in the PDB file to use (None to include all models).
    """

    print('Processing PDB:', file_name)
    
    # HO is a hydrogen attached to an oxygen. 'W' is water (infrequently found)
    atomdefs={'H':(1.0,1.00794),'HO':(1.0,1.00794),'C':(6.0,12.0107),'A':(7.0,14.00674),'N':(7.0,14.00674),'O':(8.0,15.9994),'P':(15.0,30.973761),'K':(19.0,39.0983),'S':(16.0,32.066),'W':(18.0,1.00794*2.0+15.9994),'AU':(79.0,196.96655) }
    transmap=str.maketrans("", "", "0123456789")
    
    try:
        infile = open(file_name, "r")
    except:
        raise IOError("%s is an invalid file name" % file_name)

    amin = [1.0e20, 1.0e20, 1.0e20]
    amax = [-1.0e20, -1.0e20, -1.0e20]
    atoms = {}
    coords = []

    # parse the pdb file and pull out relevant atoms
    stm = False
    for line in infile:
        # only for one atom
        if model is not None:
            if line[:5] == "MODEL":
                if int(line.split()[1]) == model:
                    stm = True
            if stm and line[:6] == "ENDMDL":
                break
            if not stm:
                continue
            
        # process all models
        if line[:4] == 'ATOM' or (line[:6] == 'HETATM' and het):
            if selected_atoms and not (line[21] in selected_atoms):
                continue
            try:
                a = line[12:14].strip().translate(transmap)
                x = float(line[30:38])
                y = float(line[38:46])
                z = float(line[46:54])
                coords.append([x, y, z])
            except:
                print("PDB Parse error:\n%s\n'%s','%s','%s'  '%s','%s','%s'\n" % (
                    line, line[12:14], line[6:11], line[22:26], line[30:38], line[38:46], line[46:54]))
                print(a, x, y, z)

            if a in atoms:
                atoms[a].append([x, y, z])
            else:
                atoms[a] = [[x, y, z]]
        
            amin = [min(x, amin[0]), min(y, amin[1]), min(z, amin[2])]
            amax = [max(x, amax[0]), max(y, amax[1]), max(z, amax[2])]


    infile.close()
   
    # Prepare coords
    coords = np.array(coords, dtype=np.float32)
    amax = np.array(amax, dtype=np.float32)
    amin = np.array(amin, dtype=np.float32)
    
    # Calculate center of mass
    cm = np.mean(coords, axis=0)  
    print('\tCenter of mass:', cm)

    # Calculate box size
    print('\tOffset:', offset)
    axis_size = np.ceil((amax - amin + 2 * offset)).astype(np.int32) 
    max_dimension = np.max(axis_size)
    box_size = [max_dimension, max_dimension, max_dimension]
    print("\tInitial box size:", box_size)

    tomo_final = np.zeros((box_size[0], box_size[1], box_size[2]), dtype=np.float32)
    
    # Modify coordinates
    for key, coord_list in atoms.items():
        # Use scalar dimensions
        print("\tAdding atoms:", str(key))
        tomo_size = np.zeros((box_size[0], box_size[1], box_size[2]), dtype=np.float32)
        c = np.array(coord_list, dtype=np.float32) 
        c -= cm
        c[:,0] += box_size[0] / 2 - 1
        c[:,1] += box_size[1] / 2 - 1
        c[:,2] += box_size[2] / 2 - 1
        value = atomdefs[key][0]
        tomo_size = add_cloud_gauss(tomo_size, c, SIGMA_FACTOR, value)
        tomo_final += tomo_size
        
    # Apply transform
    print('\tRescaling (scale=' + str(1/apix) + ')...')
    tomo_final = skimage.transform.rescale(tomo_final, scale=1/apix, order=0, anti_aliasing = True)
    if res > apix:
        f_res = res / apix
        print('\tLow-pass filter by (f_res=' + str(SIGMA_FACTOR * f_res) + ')...')
        tomo_final = scipy.ndimage.gaussian_filter(tomo_final, SIGMA_FACTOR*f_res)
    tomo_final = tomo_final.astype(np.float32)
    print("\tFinal box size: ", tomo_final.shape)
    
    # #tipificar
    # tomo_final -= tomo_final.mean()
    # tomo_final /= tomo_final.std()
    lio.write_mrc(tomo_final, output_file, apix)
    print("Finished script-----  ")
