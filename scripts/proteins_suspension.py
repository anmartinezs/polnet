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

import random

from polnet.utils import *
from polnet import lio
from polnet import tem
from polnet import poly as pp
from polnet.network import PGenHelixFiber, NetSAWLC
from polnet.stomo import MmerFile, SynthTomo, SetTomos

##### Input parameters

# Common tomogram settings
ROOT_PATH = '/fs/pool/pool-lucic2/antonio/polnet/riboprot/synth' # '/home/antonio/workspace/synth_tomo/riboprot'
NTOMOS = 10 # 12
VOI_SHAPE = (512, 512, 256) # (1856, 1856, 236) # (1856, 1856, 464) # (400, 400, 150) # (400, 400, 464) # (924, 924, 300) # vx
VOI_OFFS = ((4,508), (4,508), (4,252)) # ((4,1852), (4,1852), (4,232)) # vx
VOI_VSIZE = 10 # 2.2 # A/vx
GTRUTH_VTP_LBLS = 'gt_labels'
GTRUTH_POINTS_RAD = 35 # nm

# Proteins list
PROTEINS_LIST = ['in_shrec/5MRC_10A.pns', 'in_shrec/4V94_10A.pns', 'in_shrec/4CR2_10A.pns', 'in_shrec/1QVR_10A.pns',
                 'in_shrec/1BXN_10A.pns', 'in_shrec/3CF3_10A.pns', 'in_shrec/1U6G_10A.pns', 'in_shrec/3D2F_10A.pns',
                 'in_shrec/2CG9_10A.pns', 'in_shrec/3H84_10A.pns', 'in_shrec/3GL1_10A.pns', 'in_shrec/3QM1_10A.pns',
                 'in_shrec/1S3X_10A.pns']
                # ['in_2.2A/4V4R_2.2A.pns', 'in_2.2A/4V4R_50S_2.2A.pns', 'in_2.2A/4V4R_30S_2.2A.pns', 'in_2.2A/3J9I_2.2A.pns',
                #  'in_2.2A/3IPM_2.2A.pns', 'in_2.2A/6UTJ_2.2A.pns',
                #  'in_2.2A/5MRC_2.2A.pns', 'in_2.2A/4V94_2.2A.pns', 'in_2.2A/4CR2_2.2A.pns', 'in_2.2A/1QVR_2.2A.pns',
                #  'in_2.2A/1BXN_2.2A.pns', 'in_2.2A/3CF3_2.2A.pns', 'in_2.2A/1U6G_2.2A.pns', 'in_2.2A/3D2F_2.2A.pns',
                #  'in_2.2A/2CG9_2.2A.pns', 'in_2.2A/3H84_2.2A.pns', 'in_2.2A/3GL1_2.2A.pns', 'in_2.2A/3QM1_2.2A.pns',
                #  'in_2.2A/1S3X_2.2A.pns']
SURF_DEC = None # 0.9 # Target reduction factor for surface decimation (defatul None)

# Reconstruction tomograms
TILT_ANGS = range(-60, 61, 3) # np.arange(-60, 60, 3) # at MPI-B IMOD only works for ranges
DETECTOR_SNR = [.15, .25] # 0.2 # None
MALIGN_MN = 1
MALIGN_MX = 1.5
MALIGN_SG = 0.2

# OUTPUT FILES
OUT_DIR = ROOT_PATH + '/out_shrec_snr' # '/out_multi_v1' # '/out_v2'
TEM_DIR = OUT_DIR + '/tem'
TOMOS_DIR = OUT_DIR + '/tomos'

##### Main procedure

set_stomos = SetTomos()

# Preparing intermediate directories
clean_dir(TEM_DIR)
clean_dir(TOMOS_DIR)

# Loop for tomograms
for tomod_id in range(NTOMOS):

    # Generate the VOI and tomogram density
    voi = np.zeros(shape=VOI_SHAPE, dtype=bool)
    voi[VOI_OFFS[0][0]:VOI_OFFS[0][1], VOI_OFFS[1][0]:VOI_OFFS[1][1], VOI_OFFS[2][0]:VOI_OFFS[2][1]] = True
    tomo_den = np.zeros(shape=voi.shape, dtype=np.float32)
    synth_tomo = SynthTomo()
    poly_vtp = None

    # Loop for the list of input proteins loop
    for p_id, p_file in enumerate(PROTEINS_LIST):

        print('PROCESSING FILE:', p_file)

        # Loading the protein
        protein = MmerFile(ROOT_PATH + '/' + p_file)

        # Genrate the SAWLC network associated to the input protein
        # Polymer parameters
        model = lio.load_mrc(protein.get_mmer_svol())
        model_mask = model < protein.get_iso()
        model[model_mask] = 0
        model_surf = pp.iso_surface(model, protein.get_iso(), closed=False, normals=None)
        if SURF_DEC is not None:
            model_surf = pp.poly_decimate(model_surf, SURF_DEC)
        center = .5 * np.asarray(model.shape, dtype=float)
        # Monomer centering
        model_surf = pp.poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = pp.poly_scale(model_surf, VOI_VSIZE)
        surf_diam = pp.poly_diam(model_surf)
        pol_l_generator = PGenHelixFiber()
        # Compute protein occupancy
        if protein.get_pmer_np() is not None:
            model_voxels = (model_mask == False).sum()
            voi_voxels = voi.sum()
            protein_occ = 100 * (protein.get_pmer_np() * model_voxels) / voi_voxels
        else:
            protein_occ = protein.get_pmer_occ()

        # Network generation
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, protein.get_pmer_l() * surf_diam, model_surf, protein.get_pmer_l_max(),
                             pol_l_generator, protein_occ, protein.get_pmer_over_tol(), poly=None,
                             svol=model < protein.get_iso())
        net_sawlc.build_network()
        voi = net_sawlc.get_voi()

        # Density tomogram updating
        net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min')
        net_sawlc.insert_density_svol(model, tomo_den, VOI_VSIZE, merge='max')
        hold_vtp = net_sawlc.get_vtp()
        pp.add_label_to_poly(hold_vtp, p_id, p_name=GTRUTH_VTP_LBLS)
        if poly_vtp is not None:
            poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
        else:
            poly_vtp = hold_vtp
        synth_tomo.add_network(net_sawlc, 'Protein', lbl=p_id, code=protein.get_mmer_id())

        # DEBUG
        """lio.write_mrc(voi.astype(np.float32), ROOT_PATH + '/hold_voi.mrc')
        lio.save_vtp(model_surf, ROOT_PATH + '/hold_model.vtp')
        lio.save_vtp(poly_vtp, ROOT_PATH + '/hold_poly.vtp')
        lio.write_mrc(tomo_den, ROOT_PATH + '/hold_den.mrc')"""

    # Storing simulated density results
    tomo_voi_out = TOMOS_DIR + '/tomo_voi_' + str(tomod_id) + '.mrc'
    lio.write_mrc(voi.astype(np.int16), tomo_voi_out, v_size=VOI_VSIZE)

    # Storing simulated density results
    tomo_den_out = TOMOS_DIR + '/tomo_den_' + str(tomod_id) + '.mrc'
    lio.write_mrc(tomo_den, tomo_den_out, v_size=VOI_VSIZE)
    synth_tomo.set_den(tomo_den_out)
    poly_den_out = TOMOS_DIR + '/poly_den_' + str(tomod_id) + '.vtp'
    lio.save_vtp(poly_vtp, poly_den_out)
    synth_tomo.set_poly(poly_den_out)

    # TEM for 3D reconstructions
    temic = tem.TEM(TEM_DIR)
    vol = lio.load_mrc(tomo_den_out)
    temic.gen_tilt_series_imod(vol, TILT_ANGS, ax='Y')
    temic.add_mics_misalignment(MALIGN_MN, MALIGN_MX, MALIGN_SG)
    if DETECTOR_SNR is not None:
        snr = round((DETECTOR_SNR[1] - DETECTOR_SNR[0])*random.random() + DETECTOR_SNR[0], 2)
        temic.add_detector_noise(snr)
    temic.invert_mics_den()
    temic.set_header(data='mics', p_size=(VOI_VSIZE, VOI_VSIZE, VOI_VSIZE))
    temic.recon3D_imod()
    temic.set_header(data='rec3d', p_size=(VOI_VSIZE, VOI_VSIZE, VOI_VSIZE), origin=(0, 0, 0))
    if DETECTOR_SNR is not None:
        out_mics, out_tomo_rec = TOMOS_DIR + '/tomo_mics_' + str(tomod_id) + '_snr' + str(snr) + '.mrc', TOMOS_DIR + '/tomo_rec_' \
                                 + str(tomod_id) + '_snr' + str(snr) + '.mrc'
    else:
        out_mics, out_tomo_rec = TOMOS_DIR + '/tomo_mics_' + str(tomod_id) + '.mrc', TOMOS_DIR + '/tomo_rec_' \
                                 + str(tomod_id) + '.mrc'
    shutil.copyfile(TEM_DIR + '/out_micrographs.mrc', out_mics)
    shutil.copyfile(TEM_DIR + '/out_rec3d.mrc', out_tomo_rec)
    synth_tomo.set_mics(out_mics)
    synth_tomo.set_tomo(out_tomo_rec)

    # Update the set
    set_stomos.add_tomos(synth_tomo)

# Storing tomograms CSV file
set_stomos.save_csv(OUT_DIR + '/tomos_motif_list.csv')