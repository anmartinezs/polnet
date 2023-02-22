"""
Script for evaluating time to reach maximum occupancy
monomers
    Input:
        - Range of occupancies to test
        - Tomogram dimensions parameter
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
import csv

from polnet.utils import *
import matplotlib.pyplot as plt
from polnet import lio
from polnet import poly as pp
from polnet.network import NetSAWLC, NetSAWLCInter, NetHelixFiber, NetHelixFiberB
from polnet.polymer import FiberUnitSDimer, MTUnit, MB_DOMAIN_FIELD_STR
from polnet.stomo import MmerFile, MbFile, SynthTomo, SetTomos, HelixFile, MTFile, ActinFile, MmerMbFile
from polnet.lrandom import EllipGen, SphGen, TorGen, PGenHelixFiberB, PGenHelixFiber, SGenUniform, SGenProp, OccGen
from polnet.membrane import SetMembranes


##### Input parameters

# Common tomogram settings
ROOT_PATH = '/fs/pool/pool-lucic2/antonio/polnet/riboprot/synth' #  '/fs/pool/pool-lucic2/antonio/polnet/riboprot/synth_all' # '/home/antonio/workspace/synth_tomo/riboprot'
RANGE_OCCUPANCIES = np.linspace(1, 35, 8) # 12
VOI_SHAPE = (500, 500, 200) # vx or a path to a mask (1-foreground, 0-background) tomogram
VOI_OFFS =  ((4,496), (4,496), (4,196)) # ((4,396), (4,396), (4,232)) # ((4,1852), (4,1852), (32,432)) # ((4,1852), (4,1852), (4,232)) # vx
VOI_VSIZE = 10 # 2.2 # A/vx
GTRUTH_POINTS_RAD = 35 # nm

# Lists with the features to simulate
MEMBRANES_LIST = [] # ['in_mbs/sphere.mbs', 'in_mbs/ellipse.mbs', 'in_mbs/toroid.mbs' ]

PROTEINS_LIST = ['in_10A/4v4r_10A.pns', 'in_10A/3j9i_10A.pns', 'in_10A/4v4r_50S_10A.pns', 'in_10A/4v4r_30S_10A.pns',
                 'in_10A/6utj_10A.pns', 'in_10A/5mrc_10A.pns', 'in_10A/4v94_10A.pns', 'in_10A/4cr2_10A.pns',
                 'in_10A/3qm1_10A.pns', 'in_10A/3h84_10A.pns', 'in_10A/3gl1_10A.pns', 'in_10A/3d2f_10A.pns',
                 'in_10A/3cf3_10A.pns', 'in_10A/2cg9_10A.pns', 'in_10A/1u6g_10A.pns', 'in_10A/1s3x_10A.pns',
                 'in_10A/1qvr_10A.pns', 'in_10A/1bxn_10A.pns']

MB_PROTEINS_LIST = [] # ['in_10A/mb_5wek_10A.pms', 'in_10A/mb_4pe5_10A.pms', 'in_10A/mb_5ide_10A.pms',
                    # 'in_10A/mb_5gjv_10A.pms', 'in_10A/mb_5kxi_10A.pms', 'in_10A/mb_5tj6_10A.pms',
                    # 'in_10A/mb_5tqq_10A.pms', 'in_10A/mb_5vai_10A.pms']

HELIX_LIST = [] # ['in_helix/mt.hns', 'in_helix/actin.hns']

# Proportions list, specifies the proportion for each protein, this proportion is tried to be achieved but no guaranteed
# The toal sum of this list must be 1
PROP_LIST = None # [.4, .6]
if PROP_LIST is not None:
    assert sum(PROP_LIST) == 1

DIST_OFF = 5 # A / vx
SURF_DEC = 0.9 # Target reduction factor for surface decimation (defatul None)

# Reconstruction tomograms
TILT_ANGS = range(-60, 61, 3) # np.arange(-60, 60, 3) # at MPI-B IMOD only works for ranges
DETECTOR_SNR = [1.0, 2.0] # 0.2 # [.15, .25] # 0.2
MALIGN_MN = 1
MALIGN_MX = 1.5
MALIGN_SG = 0.2

# OUTPUT FILES
OUT_DIR = ROOT_PATH + '/tests/occ_prot' # '/only_actin' # '/out_rotations'

# OUTPUT LABELS
LBL_MB = 1
LBL_AC = 2
LBL_MT = 3
LBL_CP = 4
LBL_MP = 5
LBL_BR = 6

##### Main procedure

vx_um3 = (VOI_VSIZE * 1e-4) ** 3
target_occs, meas_occs, meas_dens, meas_times = list(), list(), list(), list()

# Intalizing the CSV output file
with open(OUT_DIR + '/occupancies.csv', 'w') as file:
    writer = csv.DictWriter(file, fieldnames=['Target_Occupancy', 'Occupancy', 'Density', 'Time'])
    writer.writeheader()

# Loop for tomograms
for occ in RANGE_OCCUPANCIES:

    hold_time = time.time()

    print('GENERATING TOMOGRAM WITH TARGET OCCUPANCY:', occ)

    # Generate the VOI and tomogram density
    if isinstance(VOI_SHAPE, str):
        voi = lio.load_mrc(VOI_SHAPE) > 0
        voi_off = np.zeros(shape=voi.shape, dtype=bool)
        voi_off[VOI_OFFS[0][0]:VOI_OFFS[0][1], VOI_OFFS[1][0]:VOI_OFFS[1][1], VOI_OFFS[2][0]:VOI_OFFS[2][1]] = True
        voi = np.logical_and(voi, voi_off)
        del voi_off
    else:
        voi = np.zeros(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFFS[0][0]:VOI_OFFS[0][1], VOI_OFFS[1][0]:VOI_OFFS[1][1], VOI_OFFS[2][0]:VOI_OFFS[2][1]] = True
        voi_inital_invert = np.invert(voi)
    voi_voxels = voi.sum()
    tomo_lbls = np.zeros(shape=VOI_SHAPE, dtype=np.float32)
    tomo_den = np.zeros(shape=voi.shape, dtype=np.float32)
    synth_tomo = SynthTomo()
    poly_vtp, mbs_vtp, skel_vtp = None, None, None
    entity_id = 1
    mb_voxels, ac_voxels, mt_voxels, cp_voxels, mp_voxels = 0, 0, 0, 0, 0
    set_mbs = None

    # Divide occupancy by the input entities
    n_entities = len(MEMBRANES_LIST) + len(HELIX_LIST) + len(PROTEINS_LIST) + len(MB_PROTEINS_LIST)
    try:
        occ_entity = occ / n_entities
    except ZeroDivisionError:
        print('ERROR: no entity introduced!')
        sys.exit()

    # Membranes loop
    count_mbs, hold_den = 0, None
    for p_id, p_file in enumerate(MEMBRANES_LIST):

        print('PROCESSING FILE:', p_file)

        # Loading the membrane file
        memb = MbFile()
        memb.load_mb_file(ROOT_PATH + '/' + p_file)

        # Membrane random generation by type
        param_rg = (memb.get_min_rad(), math.sqrt(3) * max(VOI_SHAPE) * VOI_VSIZE, memb.get_max_ecc())
        if memb.get_type() == 'sphere':
            mb_sph_generator = SphGen(radius_rg=(param_rg[0], param_rg[1]))
            set_mbs = SetMembranes(voi, VOI_VSIZE, mb_sph_generator, param_rg, memb.get_thick_rg(),
                                   memb.get_layer_s_rg(), occ_entity, memb.get_over_tol())
            set_mbs.build_set(verbosity=True)
            hold_den = set_mbs.get_tomo()
            if memb.get_den_cf_rg() is not None:
                hold_den *= mb_sph_generator.gen_den_cf(memb.get_den_cf_rg()[0], memb.get_den_cf_rg()[1])
        elif memb.get_type() == 'ellipse':
            mb_ellip_generator = EllipGen(radius_rg=param_rg[:2], max_ecc=param_rg[2])
            set_mbs = SetMembranes(voi, VOI_VSIZE, mb_ellip_generator, param_rg,  memb.get_thick_rg(),
                                   memb.get_layer_s_rg(), occ_entity, memb.get_over_tol())
            set_mbs.build_set(verbosity=True)
            hold_den = set_mbs.get_tomo()
            if memb.get_den_cf_rg() is not None:
                hold_den *= mb_ellip_generator.gen_den_cf(memb.get_den_cf_rg()[0], memb.get_den_cf_rg()[1])
        elif memb.get_type() == 'toroid':
            mb_tor_generator = TorGen(radius_rg=(param_rg[0], param_rg[1]))
            set_mbs = SetMembranes(voi, VOI_VSIZE, mb_tor_generator, param_rg, memb.get_thick_rg(), memb.get_layer_s_rg(),
                                   occ_entity, memb.get_over_tol())
            set_mbs.build_set(verbosity=True)
            hold_den = set_mbs.get_tomo()
            if memb.get_den_cf_rg() is not None:
                hold_den *= mb_tor_generator.gen_den_cf(memb.get_den_cf_rg()[0], memb.get_den_cf_rg()[1])
        else:
            print('ERROR: Membrane type', memb.get_type(), 'not recognized!')
            sys.exit()

        # Density tomogram updating
        voi = set_mbs.get_voi()
        mb_mask = set_mbs.get_tomo() > 0
        mb_mask[voi_inital_invert] = False
        tomo_lbls[mb_mask] = entity_id
        count_mbs += set_mbs.get_num_mbs()
        mb_voxels += (tomo_lbls == entity_id).sum()
        tomo_den = np.maximum(tomo_den, hold_den)
        hold_vtp = set_mbs.get_vtp()
        pp.add_label_to_poly(hold_vtp, entity_id, 'Entity', mode='both')
        pp.add_label_to_poly(hold_vtp, LBL_MB, 'Type', mode='both')
        if poly_vtp is None:
            poly_vtp = hold_vtp
            skel_vtp = hold_vtp
        else:
            poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
            skel_vtp = pp.merge_polys(skel_vtp, hold_vtp)
        synth_tomo.add_set_mbs(set_mbs, 'Membrane', entity_id, memb.get_type())
        entity_id += 1

    # Get membranes poly
    if set_mbs is not None:
        mbs_vtp = vtk.vtkPolyData()
        mbs_vtp.DeepCopy(poly_vtp)

    # Loop for Helicoidal structures
    br_vtp, hold_br_vtp, hold_br_skel_vtp = None, None, None
    count_actins, count_mts = 0, 0
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
            net_helix = NetHelixFiber(voi, VOI_VSIZE, helix.get_l() * helix.get_mmer_rad() * 2, model_surf,
                                      pol_generator, occ_entity, helix.get_min_p_len(), helix.get_hp_len(),
                                      helix.get_mz_len(), helix.get_mz_len_f(), helix.get_over_tol())
            if helix.get_min_nmmer() is not None:
                net_helix.set_min_nmmer(helix.get_min_nmmer())
            net_helix.build_network()
        elif helix.get_type() == 'actin':

            br_vtp = pp.points_to_poly_spheres(points=[[0, 0, 0],], rad=helix.get_mmer_rad())
            helix = ActinFile()
            helix.load_ac_file(ROOT_PATH + '/' + p_file)
            # Fiber unit generation
            funit = FiberUnitSDimer(helix.get_mmer_rad(), VOI_VSIZE)
            model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
            # Helix Fiber parameters model
            pol_generator = PGenHelixFiberB()
            # Network generation
            net_helix = NetHelixFiberB(voi, VOI_VSIZE, helix.get_l() * helix.get_mmer_rad() * 2, model_surf,
                                       pol_generator, occ_entity, helix.get_min_p_len(), helix.get_hp_len(),
                                       helix.get_mz_len(), helix.get_mz_len_f(), helix.get_bprop(),
                                       helix.get_p_branch(), helix.get_over_tol())
            if helix.get_min_nmmer() is not None:
                net_helix.set_min_nmmer(helix.get_min_nmmer())
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
        model_mask = model_svol < .05
        # off = .5 * np.asarray(model_svol.shape) - center
        net_helix.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min', off_svol=None)
        if helix.get_den_cf_rg() is None:
            cte_val = 1
        else:
            cte_val = pol_generator.gen_den_cf(helix.get_den_cf_rg()[0], helix.get_den_cf_rg()[1])
        net_helix.insert_density_svol(model_svol * cte_val, tomo_den, VOI_VSIZE, merge='max')
        hold_lbls = np.zeros(shape=tomo_lbls.shape, dtype=np.float32)
        net_helix.insert_density_svol(np.invert(model_mask), hold_lbls, VOI_VSIZE, merge='max')
        tomo_lbls[hold_lbls > 0] = entity_id
        # lio.write_mrc(hold_lbls.astype(np.float32), '/fs/pool/pool-lucic2/antonio/polnet/riboprot/hold.mrc')
        hold_vtp = net_helix.get_vtp()
        hold_skel_vtp = net_helix.get_skel()
        pp.add_label_to_poly(hold_vtp, entity_id, 'Entity', mode='both')
        pp.add_label_to_poly(hold_skel_vtp, entity_id, 'Entity', mode='both')
        if helix.get_type() == 'mt':
            pp.add_label_to_poly(hold_vtp, LBL_MT, 'Type', mode='both')
            pp.add_label_to_poly(hold_skel_vtp, LBL_MT, 'Type', mode='both')
            count_mts += net_helix.get_num_pmers()
            mt_voxels += (tomo_lbls == entity_id).sum()
        elif helix.get_type() == 'actin':
            pp.add_label_to_poly(hold_vtp, LBL_AC, 'Type', mode='both')
            pp.add_label_to_poly(hold_skel_vtp, LBL_AC, 'Type', mode='both')
            count_actins += net_helix.get_num_pmers()
            ac_voxels += (tomo_lbls == entity_id).sum()
            # Adding branches to the vtps but as separated type
            hold_br_vtp, hold_br_skel_vtp = net_helix.get_branches_vtp(br_vtp), net_helix.get_branches_vtp()
            pp.add_label_to_poly(hold_br_vtp, entity_id, 'Entity', mode='both')
            pp.add_label_to_poly(hold_br_skel_vtp, entity_id, 'Entity', mode='both')
            pp.add_label_to_poly(hold_br_vtp, LBL_BR, 'Type', mode='both')
            pp.add_label_to_poly(hold_br_skel_vtp, LBL_BR, 'Type', mode='both')
        if poly_vtp is None:
            poly_vtp = hold_vtp
            skel_vtp = hold_skel_vtp
        else:
            poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
            skel_vtp = pp.merge_polys(skel_vtp, hold_skel_vtp)
            # if helix.get_type() == 'actin':
            #     poly_vtp = pp.merge_polys(poly_vtp, hold_br_vtp)
            #     skel_vtp = pp.merge_polys(skel_vtp, hold_br_skel_vtp)
        synth_tomo.add_network(net_helix, 'Helix', entity_id, code=helix.get_type())
        entity_id += 1

    # Loop for the list of input proteins loop
    count_prots = 0
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
                             pol_l_generator, occ_entity, protein.get_pmer_over_tol(), poly=None,
                             svol=model < protein.get_iso())
        # net_sawlc = NetSAWLCInter(voi, VOI_VSIZE, surf_diams, model_surfs, protein.get_pmer_l_max(),
        #                           pol_l_generator, pol_s_generator, protein.get_pmer_occ(), protein.get_pmer_over_tol(),
        #                           poly=None, svols=model_masks, codes=model_codes, compaq=5.5)
        net_sawlc.build_network()

        # Density tomogram updating
        net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min')
        net_sawlc.insert_density_svol(model, tomo_den, VOI_VSIZE, merge='max')
        hold_lbls = np.zeros(shape=tomo_lbls.shape, dtype=np.float32)
        net_sawlc.insert_density_svol(np.invert(model_mask), hold_lbls, VOI_VSIZE, merge='max')
        tomo_lbls[hold_lbls > 0] = entity_id
        count_prots += net_sawlc.get_num_mmers()
        cp_voxels += (tomo_lbls == entity_id).sum()
        hold_vtp = net_sawlc.get_vtp()
        hold_skel_vtp = net_sawlc.get_skel()
        pp.add_label_to_poly(hold_vtp, entity_id, 'Entity', mode='both')
        pp.add_label_to_poly(hold_skel_vtp, entity_id, 'Entity', mode='both')
        pp.add_label_to_poly(hold_vtp, LBL_CP, 'Type', mode='both')
        pp.add_label_to_poly(hold_skel_vtp, LBL_CP, 'Type', mode='both')
        if poly_vtp is None:
            poly_vtp = hold_vtp
            skel_vtp = hold_skel_vtp
        else:
            poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
            skel_vtp = pp.merge_polys(skel_vtp, hold_skel_vtp)
        synth_tomo.add_network(net_sawlc, 'SAWLC', entity_id, code=protein.get_mmer_id())
        entity_id += 1

    # Loop for the list of input proteins loop
    count_mb_prots = 0
    if mbs_vtp is None:
        if len(MB_PROTEINS_LIST) > 0:
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
            center = protein.get_mmer_center()  # .5 * np.asarray(model.shape, dtype=float)
            if center is None:
                center = .5 * (np.asarray(model.shape, dtype=np.float) - 1)
                off = np.asarray((.0, .0, .0))
            else:
                center = np.asarray(center)
                off = .5 * np.asarray(model.shape) - center
            # Adding membrane domain to monomer surface
            mb_domain_mask = np.ones(shape=model.shape, dtype=bool)
            hold_mb_z_height = protein.get_mb_z_height()
            if hold_mb_z_height is None:
                hold_mb_z_height = int(round(center[2] + 2.5 / VOI_VSIZE))
            for z in range(hold_mb_z_height + 1, model.shape[2]):
                mb_domain_mask[:, :, z] = 0
            pp.add_sfield_to_poly(model_surf, mb_domain_mask, MB_DOMAIN_FIELD_STR, dtype='float', interp='NN', mode='points')
            # Monomer centering
            model_surf = pp.poly_translate(model_surf, -center)
            # Voxel resolution scaling
            model_surf = pp.poly_scale(model_surf, VOI_VSIZE)
            surf_diam = pp.poly_diam(model_surf)
            pol_l_generator = PGenHelixFiber()
            # Network generation
            if protein.get_pmer_reverse_normals():
                mbs_vtp = pp.poly_reverse_normals(mbs_vtp)
            net_sawlc = NetSAWLC(voi, VOI_VSIZE, protein.get_pmer_l() * surf_diam, model_surf, protein.get_pmer_l_max(),
                                 pol_l_generator, occ_entity, protein.get_pmer_over_tol(), poly=mbs_vtp,
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
            hold_lbls = np.zeros(shape=tomo_lbls.shape, dtype=np.float32)
            net_sawlc.insert_density_svol(np.invert(model_mask), hold_lbls, VOI_VSIZE, merge='max')
            tomo_lbls[hold_lbls > 0] = entity_id
            count_mb_prots += net_sawlc.get_num_mmers()
            mp_voxels += (tomo_lbls == entity_id).sum()
            hold_vtp = net_sawlc.get_vtp()
            hold_skel_vtp = net_sawlc.get_skel()
            pp.add_label_to_poly(hold_vtp, entity_id, 'Entity', mode='both')
            pp.add_label_to_poly(hold_skel_vtp, entity_id, 'Entity', mode='both')
            pp.add_label_to_poly(hold_vtp, LBL_MP, 'Type', mode='both')
            pp.add_label_to_poly(hold_skel_vtp, LBL_MP, 'Type', mode='both')
            if poly_vtp is None:
                poly_vtp = hold_vtp
                skel_vtp = hold_skel_vtp
            else:
                poly_vtp = pp.merge_polys(poly_vtp, hold_vtp)
                skel_vtp = pp.merge_polys(skel_vtp, hold_skel_vtp)
            synth_tomo.add_network(net_sawlc, 'Mb-SAWLC', entity_id, code=protein.get_mmer_id())
            entity_id += 1

    # Tomogram statistics
    meas_time = time.time() - hold_time
    mb_occ = 100. * (mb_voxels / voi_voxels)
    ac_occ = 100. * (ac_voxels / voi_voxels)
    mt_occ = 100. * (mt_voxels / voi_voxels)
    pt_occ = 100. * (cp_voxels / voi_voxels)
    mp_occ = 100. * (mp_voxels / voi_voxels)
    meas_occ = mb_occ + ac_occ + mt_occ + pt_occ + mp_occ
    mb_den = mb_voxels * vx_um3
    ac_den = ac_voxels * vx_um3
    mt_den = mt_voxels * vx_um3
    pt_den = cp_voxels * vx_um3
    mp_den = mp_voxels * vx_um3
    meas_den = mb_den + ac_den + mt_den + pt_den + mp_den
    print('\t-TOMOGRAM WITH TARGET OCCUPANCY', occ, 'DENSITY STATISTICS:')
    print('\t\t+Membranes:', count_mbs, '#, ', mb_den, 'um**3, ', mb_occ, '%')
    print('\t\t+Actin:', count_actins, '#, ', ac_den, 'um**3, ', ac_occ, '%')
    print('\t\t+Microtublues:', count_mts, '#, ', mt_den, 'um**3, ', mt_occ, '%')
    print('\t\t+Proteins:', count_prots, '#, ', pt_den, 'um**3, ', pt_occ, '%')
    print('\t\t+Membrane proteins:', count_mb_prots, '#, ', mp_den, 'um**3, ', mp_occ, '%')
    counts_total = count_mbs + count_actins + count_mts + count_prots + count_mb_prots
    total_voxels = mb_voxels + ac_voxels + mt_voxels + cp_voxels + mp_voxels
    print('\t\t+Total:', counts_total, '#, ', total_voxels * vx_um3, 'um**3, ', 100. * (total_voxels / voi_voxels), '%')
    print('\t\t+Time for generation: ', meas_time, 'secs')
    target_occs.append(occ)
    meas_occs.append(meas_occ)
    meas_dens.append(meas_den)
    meas_times.append(meas_time)

    # Storing simulated density results
    tomo_den_out = OUT_DIR + '/tomo_den_occ_' + str(occ) + '.mrc'
    lio.write_mrc(tomo_den, tomo_den_out, v_size=VOI_VSIZE)
    synth_tomo.set_den(tomo_den_out)
    poly_den_out = OUT_DIR + '/poly_den_occ_' + str(occ) + '.vtp'
    lio.save_vtp(poly_vtp, poly_den_out)
    synth_tomo.set_poly(poly_den_out)
    poly_skel_out = OUT_DIR + '/poly_skel_occ_' + str(occ) + '.vtp'
    lio.save_vtp(skel_vtp, poly_skel_out)

    # Updating the CSV output file
    with open(OUT_DIR + '/occupancies.csv', 'a') as file:
        writer = csv.DictWriter(file, fieldnames=['Target_Occupancy', 'Occupancy', 'Density', 'Time'])
        writer.writerow({'Target_Occupancy': occ, 'Occupancy': meas_occ, 'Density': meas_den, 'Time': meas_time})


# Plotting occupancy results
plt.plot(target_occs, meas_occs)
plt.ylabel('Measured occupancy [%]')
plt.xlabel('Target occupancy [%]')
plt.savefig(OUT_DIR + '/Measured_occupancy.png')
plt.close()
plt.plot(target_occs, meas_dens)
plt.ylabel('Measured density [$um^3$]')
plt.xlabel('Target occupancy [%]')
plt.savefig(OUT_DIR + '/Measured_density.png')
plt.close()
plt.plot(target_occs, meas_times)
plt.ylabel('Time [secs]')
plt.xlabel('Target occupancy [%]')
plt.savefig(OUT_DIR + '/Measured_time.png')
plt.close()


print('Successfully terminated. (' + time.strftime("%c") + ')')