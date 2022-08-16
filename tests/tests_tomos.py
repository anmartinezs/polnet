"""
Test for checking a single simple tomogram with a categories of structures, currently:
    - Membranes with different local curvatures:
        + Sphere
        + Ellipsoid
        + Toroid
    - Independant fibers network: microtubule-like
    - Branched filament network: actin-like
    - Protein networks (cytosolic and membrane-bound):
        + Completely Spatial Randomness with Volume exclusion: Uncorrelated macromolecules
        + Self-Avoiding Worm-Like Chain (SAWLC): Polyribosomes
"""
import copy
from unittest import TestCase

from polnet.membrane import SetMembranes
from polnet.network import NetSAWLC, NetHelixFiber, NetHelixFiberB
from polnet.lrandom import EllipGen, SphGen, TorGen, PGenHelixFiber, PGenHelixFiberB
from polnet.polymer import FiberUnitSDimer, MTUnit
from polnet.polymer import MB_DOMAIN_FIELD_STR
from polnet.poly import *
from polnet.lio import *
from polnet.utils import *
from polnet.affine import *

# Common settings
VOI_SHAPE = (400, 400, 150) # (924, 924, 300) # vx
VOI_OFF = 4 # vx
VOI_VSIZE = 13.68 # A/vx
GTRUTH_VTP_LBLS = 'gt_labels'
GTRUTH_POINTS_RAD = 35 # nm

# NetSAWLC settings
MMER_MRC = './in/emd_5225_13.68A_30.mrc'
RIBO_MMER_ISO = .543
RIBO_MMER_MRC = './in/emd_5225_13.68A_30.mrc'
RIBO_PMER_L = 1.2 # Number of times wrt MMER diameters
RIBO_PMER_OCC = 1 # 5 # %
PMER_L_MAX = 3000
PMER_OVER_TOL = 0.001

# NetHelixFiber settings
A_MMER_ISO = 0.9
A_MMER_RAD = 25 # A
A_PMER_OCC = 0.1 # 5 # %
A_MIN_P_LEN = 17.7e4 # 5200e4 # Actin filament 17.7um according to 10.1083/jcb.130.4.909
A_HP_LEN = 720 # A
A_MZ_LEN = 50 # A
A_MZ_LEN_F = 0.2
A_BPROP = 1 / 2
A_MAX_P_BRANCH = 2

# NetHelixFiber for MT settings
MT_MMER_ISO = 0.9
MT_MMER_RAD = 40 # A according to http://web.physics.ucsb.edu/~jennyr/research1.html
MT_RAD = 100.5 # A according to http://web.physics.ucsb.edu/~jennyr/research1.html
MT_NUNITS = 13
MT_PMER_OCC = 0.1 # 5 # %
MT_MIN_P_LEN = 5200e4 # MT 5200um according to 10.1083/jcb.130.4.909
MT_HP_LEN = 216000 # A
MT_MZ_LEN = 50 # 80 # A
MT_MZ_LEN_F = 0.4

# Generic settings for membranes
MB_OCC = 0.002 # %
MB_THICK_RG = (25, 35) # A
MB_LAYER_S_RG = (5, 10) # A
MB_MAX_ECC = .75
MB_OVER_TOL = 0 # 1e-20
MB_MIN_RAD = 75 # A
MB_VOI_GROWTH = 2 # voxels

# Membrane bound protein
MMER_MRC = './in/pdb_5wek_13.68A_uncentered.mrc'
MB_MMER_ISO = .046
MMER_CENTER = [13, 14.5, 11] # [57, 76, 36] # voxels
MB_Z_HEIGHT = 16 # voxels
PMER_L = 1.2 # Number of times wrt MMER diameters
PMER_OCC = 0.03 # 5 # %
PMER_L_MAX = 3000
PMER_OVER_TOL = 0.01
PMER_REVERSE_NORMALS = True

# Paths to tomograms outputs
STOMO_TOMO_OUT = './out/stomo_density.mrc'
STOMO_VOI_OUT = './out/stomo_voi.mrc'
STOMO_VTP_OUT = './out/stomo.vtp'
STOMO_GTRUTH_VTP_OUT = './out/stomo_gtruth.vtp'


class TestTomos(TestCase):

    def test_build_tomo(self):

        # return

        #### MEMBRANES Y MEMBRANE-BOUND PROTEINS

        # Generate the VOI and tomogram density
        voi = np.ones(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFF:VOI_SHAPE[0]-VOI_OFF, VOI_OFF:VOI_SHAPE[1]-VOI_OFF, VOI_OFF:VOI_SHAPE[2]-VOI_OFF] = True
        tomo = np.zeros(shape=voi.shape, dtype=np.float32)

        # Membrane random generation parametrization
        param_rg = (MB_MIN_RAD, math.sqrt(3) * max(VOI_SHAPE) * VOI_VSIZE, MB_MAX_ECC)
        mb_ellip_generator = EllipGen(radius_rg=param_rg[:2], max_ecc=param_rg[2])
        mb_sph_generator = SphGen(radius_rg=(.5*param_rg[0], .5*param_rg[1]))
        mb_tor_generator = TorGen(radius_rg=(.5*param_rg[0], .5*param_rg[1]))

        # Network generation for Ellipsoids
        print('Generating Ellipsoids...')
        set_mbs = SetMembranes(voi, VOI_VSIZE, mb_ellip_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, 3 * MB_OCC,
                               MB_OVER_TOL, MB_VOI_GROWTH)
        set_mbs.build_set(verbosity=True)
        voi = set_mbs.get_voi()
        tomo = np.maximum(tomo, set_mbs.get_tomo())
        gtruth_vtp = set_mbs.get_vtp()
        add_label_to_poly(gtruth_vtp, 1, p_name=GTRUTH_VTP_LBLS)
        poly_vtp = vtk.vtkPolyData()
        poly_vtp.DeepCopy(gtruth_vtp)

        # Insert membrane bound densities in a Polymer
        # Polymer parameters
        model = load_mrc(MMER_MRC)
        model_mask = model < MB_MMER_ISO
        model[model_mask] = 0
        model_surf = iso_surface(model, MB_MMER_ISO, closed=False, normals=None)
        center = np.asarray(MMER_CENTER) # .5 * np.asarray(model.shape, dtype=float)
        off = .5 * np.asarray(model.shape) - center
        # Adding membrane domain to monomer surface
        mb_domain_mask = np.ones(shape=model.shape, dtype=bool)
        for z in range(MB_Z_HEIGHT+1, model.shape[2]):
            mb_domain_mask[:, :, z] = 0
        add_sfield_to_poly(model_surf, mb_domain_mask, MB_DOMAIN_FIELD_STR, dtype='float', interp='NN', mode='points')
        # Monomer centering
        model_surf = poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = poly_scale(model_surf, VOI_VSIZE)
        surf_diam = poly_diam(model_surf)
        pol_l_generator = PGenHelixFiber()
        # Network generation
        mb_poly = set_mbs.get_vtp()
        if PMER_REVERSE_NORMALS:
            mb_poly = poly_reverse_normals(mb_poly)
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, PMER_L * surf_diam, model_surf, PMER_L_MAX, pol_l_generator, PMER_OCC,
                             PMER_OVER_TOL, poly=mb_poly, svol=model < MB_MMER_ISO)
        net_sawlc.build_network()
        voi = net_sawlc.get_voi()

        # Adding densities
        net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min', off_svol=off)
        net_sawlc.insert_density_svol(model, tomo, VOI_VSIZE, merge='max', off_svol=off)
        hold_gtruth_vtp = net_sawlc.get_skel(add_lines=False, verts_rad=GTRUTH_POINTS_RAD)
        add_label_to_poly(hold_gtruth_vtp, 7, p_name=GTRUTH_VTP_LBLS)
        gtruth_vtp = merge_polys(gtruth_vtp, hold_gtruth_vtp)
        hold_vtp = net_sawlc.get_vtp()
        add_label_to_poly(hold_vtp, 7, p_name=GTRUTH_VTP_LBLS)
        poly_vtp = merge_polys(poly_vtp, hold_vtp)

        # Network generation for Spheres
        print('Generating Spheres...')
        set_mbs = SetMembranes(voi, VOI_VSIZE, mb_sph_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
                               MB_OVER_TOL)
        set_mbs.build_set(verbosity=True)
        voi = set_mbs.get_voi()
        tomo = np.maximum(tomo, set_mbs.get_tomo())
        hold_gtruth_vtp = set_mbs.get_vtp()
        add_label_to_poly(hold_gtruth_vtp, 2, p_name=GTRUTH_VTP_LBLS)
        gtruth_vtp = merge_polys(gtruth_vtp, hold_gtruth_vtp)
        poly_vtp = merge_polys(poly_vtp, hold_gtruth_vtp)

        # Network generation for Toroids
        print('Generating Toroids...')
        set_mbs = SetMembranes(voi, VOI_VSIZE, mb_tor_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
                               MB_OVER_TOL)
        set_mbs.build_set(verbosity=True)
        voi = set_mbs.get_voi()
        tomo = np.maximum(tomo, set_mbs.get_tomo())
        hold_gtruth_vtp = set_mbs.get_vtp()
        add_label_to_poly(hold_gtruth_vtp, 3, p_name=GTRUTH_VTP_LBLS)
        gtruth_vtp = merge_polys(gtruth_vtp, hold_gtruth_vtp)
        poly_vtp = merge_polys(poly_vtp, hold_gtruth_vtp)

        #### BRANCHED ACTIN-LIKE

        # Fiber unit generation
        funit = FiberUnitSDimer(A_MMER_RAD, VOI_VSIZE)
        model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
        model_mask = model_mask < A_MMER_ISO

        # Helix Fiber parameters model
        pol_generator = PGenHelixFiberB()

        # Network generation
        net_helix = NetHelixFiberB(voi, VOI_VSIZE, PMER_L * A_MMER_RAD * 2, model_surf, pol_generator, A_PMER_OCC,
                                   A_MIN_P_LEN, A_HP_LEN, A_MZ_LEN, A_MZ_LEN_F, A_BPROP, A_MAX_P_BRANCH, PMER_OVER_TOL)
        net_helix.build_network()
        voi = net_helix.get_voi()

        # Density tomogram generation
        net_helix.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min', off_svol=off)
        net_helix.insert_density_svol(model_svol, tomo, VOI_VSIZE, merge='max')
        hold_gtruth_vtp = net_helix.get_skel()
        add_label_to_poly(hold_gtruth_vtp, 4, p_name=GTRUTH_VTP_LBLS)
        gtruth_vtp = merge_polys(gtruth_vtp, hold_gtruth_vtp)
        hold_vtp = net_helix.get_vtp()
        add_label_to_poly(hold_vtp, 4, p_name=GTRUTH_VTP_LBLS)
        poly_vtp = merge_polys(poly_vtp, hold_vtp)

        #### MICROTUBULE-LIKE

        # Fiber unit generation
        funit = MTUnit(MT_MMER_RAD, MT_RAD, MT_NUNITS, VOI_VSIZE)
        model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
        model_mask = model_mask < MT_MMER_ISO

        # Helix Fiber parameters model
        pol_generator = PGenHelixFiber()

        # Network generation
        net_helix = NetHelixFiber(voi, VOI_VSIZE, PMER_L * MT_MMER_RAD * 2, model_surf, pol_generator, MT_PMER_OCC,
                                  MT_MIN_P_LEN, MT_HP_LEN, MT_MZ_LEN, MT_MZ_LEN_F, PMER_OVER_TOL)
        net_helix.build_network()
        voi = net_helix.get_voi()

        # Density tomogram generation
        net_helix.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min', off_svol=off)
        net_helix.insert_density_svol(model_svol, tomo, VOI_VSIZE, merge='max')
        hold_gtruth_vtp = net_helix.get_skel()
        add_label_to_poly(hold_gtruth_vtp, 5, p_name=GTRUTH_VTP_LBLS)
        gtruth_vtp = merge_polys(gtruth_vtp, hold_gtruth_vtp)
        hold_vtp = net_helix.get_vtp()
        add_label_to_poly(hold_vtp, 5, p_name=GTRUTH_VTP_LBLS)
        poly_vtp = merge_polys(poly_vtp, hold_vtp)

        #### POLYRIBOSOME

        # Polymer parameters
        model = load_mrc(RIBO_MMER_MRC)
        model_mask = model < RIBO_MMER_ISO
        model[model_mask] = 0
        model_surf = iso_surface(model, RIBO_MMER_ISO, closed=False, normals=None)
        center = .5 * np.asarray(model.shape, dtype=float)
        # Monomer centering
        model_surf = poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = poly_scale(model_surf, VOI_VSIZE)
        surf_diam = poly_diam(model_surf)
        pol_l_generator = PGenHelixFiber()

        # Network generation
        # net_sawlc = NetSAWLC(voi, VOI_VSIZE, PMER_L * surf_diam, model_surf, PMER_L_MAX, pol_l_generator, RIBO_PMER_OCC,
        #                      PMER_OVER_TOL)
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, PMER_L * surf_diam, model_surf, PMER_L_MAX, pol_l_generator, RIBO_PMER_OCC,
                             PMER_OVER_TOL, poly=None, svol=model < RIBO_MMER_ISO)
        net_sawlc.build_network()
        voi = net_sawlc.get_voi()

        # Density tomogram generation
        net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min')
        net_sawlc.insert_density_svol(model, tomo, VOI_VSIZE, merge='max')
        hold_gtruth_vtp = net_sawlc.get_skel(add_lines=False, verts_rad=GTRUTH_POINTS_RAD)
        add_label_to_poly(hold_gtruth_vtp, 6, p_name=GTRUTH_VTP_LBLS)
        gtruth_vtp = merge_polys(gtruth_vtp, hold_gtruth_vtp)
        hold_vtp = net_sawlc.get_vtp()
        add_label_to_poly(hold_vtp, 6, p_name=GTRUTH_VTP_LBLS)
        poly_vtp = merge_polys(poly_vtp, hold_vtp)

        # Storing the final results
        write_mrc(tomo, STOMO_TOMO_OUT, v_size=VOI_VSIZE)
        write_mrc(net_sawlc.get_voi().astype(np.float32), STOMO_VOI_OUT, v_size=VOI_VSIZE)
        save_vtp(gtruth_vtp, STOMO_GTRUTH_VTP_OUT)
        save_vtp(poly_vtp, STOMO_VTP_OUT)