"""
Script for generating tomograms simulating a suspension with different proteins distributed according SAWLC network
model
    Input:
        - Number of tomograms to simulate
        - Tomogram dimensions parameter
        - Tomogram maximum occupancy
        - Cytosolic protein networks:
            + List of proteins files
        - 3D reconstruction paramaters
    Output:
        - The simulated density maps
        - The 3D reconstructed tomograms
        - Micrograph stacks
        - Polydata files
        - STAR file mapping particle coordinates and orientations with tomograms
"""

__author__ = 'Antonio Martinez-Sanchez'

import shutil

from polnet.utils import *
from polnet import lio
from polnet import tem
from polnet import poly as pp
from polnet.network import PGenHelixFiber, NetSAWLC
from polnet.stomo import MmerFile, SynthTomo, SetTomos

##### Input parameters

# Common tomogram settings
ROOT_PATH = ''
NTOMOS = 2
VOI_SHAPE = (400, 400, 150) # (924, 924, 300) # vx
VOI_OFF = 4 # vx
VOI_VSIZE = 2.2 # A/vx
GTRUTH_VTP_LBLS = 'gt_labels'
GTRUTH_POINTS_RAD = 35 # nm

# Proteins list
PROTEINS_LIST = ['', '']

# Reconstruction tomograms
TILT_ANGS = np.arange(-60, 60, 2)
DETECTOR_SNR = 0.1

# OUTPUT FILES
OUT_DIR = ROOT_PATH + ''
TEM_DIR = OUT_DIR + '/tem'
TOMOS_DIR = OUT_DIR + '/tem'

##### Main procedure

set_stomos = SetTomos()

# Preparing intermediate directories
clean_dir(TEM_DIR)
clean_dir(TOMOS_DIR)

# Loop for tomograms
for tomod_id in range(NTOMOS):

    # Generate the VOI and tomogram density
    voi = np.ones(shape=VOI_SHAPE, dtype=bool)
    voi[VOI_OFF:VOI_SHAPE[0] - VOI_OFF, VOI_OFF:VOI_SHAPE[1] - VOI_OFF, VOI_OFF:VOI_SHAPE[2] - VOI_OFF] = True
    tomo_den = np.zeros(shape=voi.shape, dtype=np.float32)
    synth_tomo = SynthTomo()

    # Loop for the list of input proteins loop
    for p_id, p_file in enumerate(PROTEINS_LIST):

        # Loading the protein
        protein = MmerFile(ROOT_PATH + '/' + p_file)

        # Genrate the SAWLC network associated to the input protein
        # Polymer parameters
        model = lio.load_mrc(protein.get_mmer_svol())
        model_mask = model < protein.get_iso()
        model[model_mask] = 0
        model_surf = pp.iso_surface(model, protein.get_iso(), closed=False, normals=None)
        center = .5 * np.asarray(model.shape, dtype=float)
        # Monomer centering
        model_surf = pp.poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = pp.poly_scale(model_surf, VOI_VSIZE)
        surf_diam = pp.poly_max_distance(model_surf)
        pol_l_generator = PGenHelixFiber()

        # Network generation
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, protein.get_pmer_l() * surf_diam, model_surf, protein.get_pmer_l_max(),
                             pol_l_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(), poly=None,
                             svol=model < protein.get_iso())
        net_sawlc.build_network()
        voi = net_sawlc.get_voi()

        # Density tomogram updating
        net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min')
        net_sawlc.insert_density_svol(model, tomo_den, VOI_VSIZE, merge='max')
        hold_vtp = net_sawlc.get_vtp()
        pp.add_label_to_poly(hold_vtp, p_id, p_name=GTRUTH_VTP_LBLS)
        poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
        synth_tomo.add_network(net_sawlc, 'Protein', p_id, protein.get_mmer_id())

    # Storing simulated density results
    tomo_den_out = TOMOS_DIR + '/tomo_den_' + p_id + '.mrc'
    lio.write_mrc(tomo_den, tomo_den_out, v_size=VOI_VSIZE)
    synth_tomo.set_den(tomo_den_out)
    poly_den_out = TOMOS_DIR + '/poly_den_' + p_id + '.mrc'
    pp.save_vtp(poly_vtp, poly_den_out)
    synth_tomo.set_poly(poly_den_out)

    # TEM for 3D reconstructions
    temic = tem.TEM(TEM_DIR)
    vol = lio.load_mrc(tomo_den_out)
    temic.gen_tilt_series_imod(vol, TILT_ANGS)
    temic.add_detector_noise(DETECTOR_SNR)
    temic.invert_mics_den()
    temic.recon3D_imod()
    out_mics, out_tomo_rec = TOMOS_DIR + '/tomo_mics_' + p_id + '.mrc', TOMOS_DIR + '/tomo_rec_' + p_id + '.mrc'
    shutil.copyfile(TEM_DIR + '/out_micrographs.mrc', out_mics)
    shutil.copyfile(TEM_DIR + '/out_rec3d.mrc', out_tomo_rec)
    synth_tomo.set_mics(out_mics)
    synth_tomo.set_tomo(out_tomo_rec)

    # Update the set
    set_stomos.add_tomo(synth_tomo)

# Storing tomograms CSV file
set_stomos.save_csv(OUT_DIR + '/tomos_motif_list.csv')