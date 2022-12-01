"""
Script for generating tomograms simulating all features available
monomers
    Input:
        - Number of tomograms to simulate
        - Tomogram dimensions parameter
        - Tomogram maximum occupancy
        - Features to simulate:
            + Membranes
            + Polymers:
                + Helicoidal fibers
                + Globular protein clusters
        - 3D reconstruction paramaters
    Output:
        - The simulated density maps
        - The 3D reconstructed tomograms
        - Micrograph stacks
        - Polydata files
        - STAR file mapping particle coordinates and orientations with tomograms
"""

__author__ = 'Antonio Martinez-Sanchez'

import sys
import time
import random

from polnet.utils import *
from polnet import lio
from polnet import tem
from polnet import poly as pp
from polnet.network import NetSAWLC, NetSAWLCInter, NetHelixFiber, NetHelixFiberB
from polnet.polymer import FiberUnitSDimer, MTUnit, MB_DOMAIN_FIELD_STR
from polnet.stomo import MmerFile, MbFile, SynthTomo, SetTomos, HelixFile, MTFile, ActinFile, MmerMbFile
from polnet.lrandom import EllipGen, SphGen, TorGen, PGenHelixFiberB, PGenHelixFiber, SGenUniform, SGenProp
from polnet.membrane import SetMembranes


##### Input parameters

# Common tomogram settings
ROOT_PATH = '/fs/pool/pool-lucic2/antonio/polnet/riboprot/synth_all' # '/home/antonio/workspace/synth_tomo/riboprot'
NTOMOS = 1 # 10 # 12
VOI_SHAPE = (400, 400, 236) # (1856, 1856, 464) # (1856, 1856, 236) # (1856, 1856, 464) # (400, 400, 464) # (924, 924, 300) # vx
VOI_OFFS =  ((4,396), (4,396), (4,232)) # ((4,1852), (4,1852), (32,432)) # ((4,1852), (4,1852), (4,232)) # vx
VOI_VSIZE = 2.2 # A/vx
GTRUTH_POINTS_RAD = 35 # nm

# Lists with the features to simulate
MEMBRANES_LIST = ['in/toroid.mbs']
PROTEINS_LIST = ['in/ribo_v2.pns', 'in/prot_v2.pns'] # ['in/ribo_v2.pns', 'in/prot_v2.pns', 'in/ribo_30S_v2.pns', 'in/ribo_50S_v2.pns', 'in/prot_sc_v2.pns',
                # 'in/prot_dc_v2.pns']
MB_PROTEINS_LIST = ['in/mb_5wek.pms']
HELIX_LIST = ['in/mt.hns', 'in/actin.hns']

# Proportions list, specifies the proportion for each protein, this proportion is tried to be achieved but no guaranteed
# The toal sum of this list must be 1
PROP_LIST = None # [.4, .6]
if PROP_LIST is not None:
    assert sum(PROP_LIST) == 1

DIST_OFF = 5 # A / vx
SURF_DEC = 0.9 # Target reduction factor for surface decimation (defatul None)

# Reconstruction tomograms
TILT_ANGS = range(-60, 61, 3) # np.arange(-60, 60, 3) # at MPI-B IMOD only works for ranges
DETECTOR_SNR = 0.2 # [.15, .25] # 0.2
MALIGN_MN = 1
MALIGN_MX = 1.5
MALIGN_SG = 0.2

# OUTPUT FILES
OUT_DIR = ROOT_PATH + '/out_rotations'
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
    poly_vtp, mbs_vtp = None, None
    entity_id = 1

    # Membranes loop
    for p_id, p_file in enumerate(MEMBRANES_LIST):

        print('PROCESSING FILE:', p_file)

        # Loading the membrane file
        memb = MbFile()
        memb.load_mb_file(ROOT_PATH + '/' + p_file)

        # Membrane random generation by type
        param_rg = (memb.get_min_rad(), math.sqrt(3) * max(VOI_SHAPE) * VOI_VSIZE, memb.get_max_ecc())
        if memb.get_type() == 'sphere':
            mb_sph_generator = SphGen(radius_rg=(.5 * param_rg[0], .5 * param_rg[1]))
            set_mbs = SetMembranes(voi, VOI_VSIZE, mb_sph_generator, param_rg, memb.get_thick_rg(),
                                   memb.get_layer_s_rg(), memb.get_occ(), memb.get_over_tol())
            set_mbs.build_set(verbosity=True)
        elif memb.get_type() == 'ellipse':
            mb_ellip_generator = EllipGen(radius_rg=param_rg[:2], max_ecc=param_rg[2])
            set_mbs = SetMembranes(voi, VOI_VSIZE, mb_ellip_generator, param_rg,  memb.get_thick_rg(),
                                   memb.get_layer_s_rg(), memb.get_occ(), memb.get_over_tol())
            set_mbs.build_set(verbosity=True)
        elif memb.get_type() == 'toroid':
            mb_tor_generator = TorGen(radius_rg=(.5 * param_rg[0], .5 * param_rg[1]))
            set_mbs = SetMembranes(voi, VOI_VSIZE, mb_tor_generator, param_rg, memb.get_thick_rg(), memb.get_layer_s_rg(),
                                   memb.get_occ(), memb.get_over_tol())
            set_mbs.build_set(verbosity=True)
        else:
            print('ERROR: Membrane type', memb.get_type(), 'not recognized!')
            sys.exit()

        # Density tomogram updating
        voi = set_mbs.get_voi()
        tomo_den = np.maximum(tomo_den, set_mbs.get_tomo())
        hold_vtp = set_mbs.get_vtp()
        if poly_vtp is None:
            poly_vtp = hold_vtp
        else:
            poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
        synth_tomo.add_set_mbs(set_mbs, 'Membrane', entity_id, 'Mb-' + memb.get_type())
        entity_id += 1

    # Get membranes poly
    if set_mbs is not None:
        mbs_vtp = vtk.vtkPolyData()
        mbs_vtp.DeepCopy(poly_vtp)

    # Loop for Helicoidal structures
    for p_id, p_file in enumerate(HELIX_LIST):

        print('PROCESSING FILE:', p_file)

        # Loading the helix file
        helix = HelixFile()
        helix.load_hx_file(ROOT_PATH + '/' + p_file)

        # Helicoida random generation by type
        if helix.get_type() == 'mt':

            helix = MTFile()
            helix.load_mt_file(ROOT_PATH + '/' + p_file)
            # Fiber unit generation
            funit = MTUnit(helix.get_mmer_rad(), helix.get_rad(), helix.get_nunits(), VOI_VSIZE)
            model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
            # Helix Fiber parameters model
            pol_generator = PGenHelixFiber()
            # Network generation
            net_helix = NetHelixFiber(voi, VOI_VSIZE, helix.get_hp_len() * helix.get_mmer_rad() * 2, model_surf,
                                      pol_generator, helix.get_occ(), helix.get_min_p_len(), helix.get_hp_len(),
                                      helix.get_mz_len(), helix.get_mz_len_f(), helix.get_over_tol())
            net_helix.build_network()
        elif helix.get_type() == 'actin':

            helix = ActinFile()
            helix.load_ac_file(ROOT_PATH + '/' + p_file)
            # Fiber unit generation
            funit = FiberUnitSDimer(helix.get_mmer_rad(), VOI_VSIZE)
            model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
            # Helix Fiber parameters model
            pol_generator = PGenHelixFiberB()
            # Network generation
            net_helix = NetHelixFiberB(voi, VOI_VSIZE, helix.get_hp_len() * helix.get_mmer_rad() * 2, model_surf,
                                       pol_generator, helix.get_occ(), helix.get_min_p_len(), helix.get_hp_len(),
                                       helix.get_mz_len(), helix.get_mz_len_f(), helix.get_bprop(),
                                       helix.get_p_branch(), helix.get_over_tol())
            net_helix.build_network()
        else:
            print('ERROR: Helicoidal type', helix.get_type(), 'not recognized!')
            sys.exit()

        # # DEBUG
        # lio.save_vtp(funit.get_vtp(), ROOT_PATH + '/hold_funit.vtp')
        # lio.save_vtp(net_helix.get_vtp(), ROOT_PATH + '/hold.vtp')

        # Density tomogram updating
        # voi = net_helix.get_voi()
        # tomo_den = np.maximum(tomo_den, net_helix.get_tomo())
        model_mask = model_svol < .25
        # off = .5 * np.asarray(model_svol.shape) - center
        net_helix.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min', off_svol=None)
        net_helix.insert_density_svol(model_svol, tomo_den, VOI_VSIZE, merge='max')
        hold_vtp = net_helix.get_vtp()
        if poly_vtp is None:
            poly_vtp = hold_vtp
        else:
            poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
        synth_tomo.add_network(net_helix, 'Helix', entity_id, code='Hx-' + helix.get_type())
        entity_id += 1

    # Loop for the list of input proteins loop
    model_surfs, models, model_masks, model_codes = list(), list(), list(), list()
    for p_id, p_file in enumerate(PROTEINS_LIST):

        print('PROCESSING FILE:', p_file)

        # Loading the protein
        protein = MmerFile(ROOT_PATH + '/' + p_file)

        # Genrate the SAWLC network associated to the input protein
        # Polymer parameters
        model = lio.load_mrc(protein.get_mmer_svol())
        model = lin_map(model, lb=0, ub=1)
        model = vol_cube(model)
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
        model_surfs.append(model_surf)
        surf_diam = pp.poly_diam(model_surf) * protein.get_pmer_l()
        models.append(model)
        model_masks.append(model_mask)
        model_codes.append(protein.get_mmer_id())

        # Network generation
        pol_l_generator = PGenHelixFiber()
        if PROP_LIST is None:
            pol_s_generator = SGenUniform()
        else:
            assert len(PROP_LIST) == len(PROTEINS_LIST)
            pol_s_generator = SGenProp(PROP_LIST)
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, protein.get_pmer_l() * surf_diam, model_surf, protein.get_pmer_l_max(),
                             pol_l_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(), poly=None,
                             svol=model < protein.get_iso())
        # net_sawlc = NetSAWLCInter(voi, VOI_VSIZE, surf_diams, model_surfs, protein.get_pmer_l_max(),
        #                           pol_l_generator, pol_s_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(),
        #                           poly=None, svols=model_masks, codes=model_codes, compaq=5.5)
        net_sawlc.build_network()

        # Density tomogram updating
        net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min')
        net_sawlc.insert_density_svol(model, tomo_den, VOI_VSIZE, merge='max')
        hold_vtp = net_sawlc.get_vtp()
        if poly_vtp is None:
            poly_vtp = hold_vtp
        else:
            poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
        synth_tomo.add_network(net_helix, 'Helix', entity_id, code='Hx-' + helix.get_type())
        entity_id += 1

    # Loop for the list of input proteins loop
    if mbs_vtp is None:
        print('WARNING: membrane proteins can not inserted because there is no membrane surfaces!')
    else:

        model_surfs, surf_diams, models, model_masks, model_codes = list(), list(), list(), list(), list()
        for p_id, p_file in enumerate(MB_PROTEINS_LIST):

            print('PROCESSING FILE:', p_file)

            # Loading the membrane protein
            protein = MmerMbFile(ROOT_PATH + '/' + p_file)

            # Insert membrane bound densities in a Polymer
            # Polymer parameters
            model = lio.load_mrc(protein.get_mmer_svol())
            model = lin_map(model, lb=0, ub=1)
            model_mask = model < protein.get_iso()
            model[model_mask] = 0
            model_surf = iso_surface(model, protein.get_iso(), closed=False, normals=None)
            center = np.asarray(protein.get_mmer_center())  # .5 * np.asarray(model.shape, dtype=float)
            off = .5 * np.asarray(model.shape) - center
            # Adding membrane domain to monomer surface
            mb_domain_mask = np.ones(shape=model.shape, dtype=bool)
            for z in range(protein.get_mb_z_height() + 1, model.shape[2]):
                mb_domain_mask[:, :, z] = 0
            pp.add_sfield_to_poly(model_surf, mb_domain_mask, MB_DOMAIN_FIELD_STR, dtype='float', interp='NN', mode='points')
            # Monomer centering
            model_surf = pp.poly_translate(model_surf, -center)
            # Voxel resolution scaling
            model_surf = pp.poly_scale(model_surf, VOI_VSIZE)
            surf_diam = pp.poly_diam(model_surf)
            pol_l_generator = PGenHelixFiber()
            # Network generation
            mb_poly = set_mbs.get_vtp()
            if protein.get_pmer_reverse_normals():
                mb_poly = pp.poly_reverse_normals(mb_poly)
            net_sawlc = NetSAWLC(voi, VOI_VSIZE, protein.get_pmer_l() * surf_diam, model_surf, protein.get_pmer_l_max(),
                                 pol_l_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(), poly=mb_poly,
                                 svol=model < protein.get_iso())
            # net_sawlc = NetSAWLCInter(voi, VOI_VSIZE, protein.get_pmer_l() * surf_diam, model_surf, protein.get_pmer_l_max(),
            #                      pol_l_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(), poly=mb_poly,
            #                      svol=model < protein.get_iso())
            net_sawlc.build_network()
            # voi = net_sawlc.get_voi()

            # lio.write_mrc(voi.astype(np.float32), ROOT_PATH + '/hold_voi.mrc')
            # lio.write_mrc(set_mbs.get_tomo().astype(np.float32), ROOT_PATH + '/hold_den.mrc')

            # Density tomogram updating
            net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min')
            net_sawlc.insert_density_svol(model, tomo_den, VOI_VSIZE, merge='max')
            hold_vtp = net_sawlc.get_vtp()
            if poly_vtp is None:
                poly_vtp = hold_vtp
            else:
                poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
            synth_tomo.add_network(net_sawlc, 'Helix', entity_id, code='Hx-' + helix.get_type())
            entity_id += 1

        # # Tomogram statistics
        # print('\t-TOMOGRAM', str(tomod_id), 'DENSITY STATISTICS:')
        # total_counts, total_vol = 0, 0
        # prot_counts = net_sawlc.count_proteins()
        # voi_vol = voi.sum() * VOI_VSIZE**3
        # for id in range(len(PROTEINS_LIST)):
        #     prot_vol = pp.poly_volume(model_surfs[id]) * prot_counts[id]
        #     print('\t\t+Protein', id, ':', prot_counts[id], '#, ', prot_vol, 'A**3, ', 100. * (prot_vol / voi_vol), '%')
        #     total_counts += prot_counts[id]
        #     total_vol += prot_vol
        # print('\t\t+Total:', total_counts, '#, ', voi_vol, 'A**3, ', 100. * (total_vol / voi_vol), '%')

        # DEBUG
        lio.write_mrc(voi.astype(np.float32), ROOT_PATH + '/hold_voi.mrc')
        lio.save_vtp(model_surf, ROOT_PATH + '/hold_model.vtp')
        lio.save_vtp(poly_vtp, ROOT_PATH + '/hold_poly.vtp')
        lio.write_mrc(tomo_den, ROOT_PATH + '/hold_den.mrc')

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
        if hasattr(DETECTOR_SNR, '__len__'):
            if len(DETECTOR_SNR) >= 2:
                snr = round((DETECTOR_SNR[1] - DETECTOR_SNR[0])*random.random() + DETECTOR_SNR[0], 2)
            else:
                snr = DETECTOR_SNR[0]
        else:
            snr = DETECTOR_SNR
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

print('Successfully terminated. (' + time.strftime("%c") + ')')