from unittest import TestCase

import numpy as np

from polnet.network import NetSAWLC, NetHelixFiber, NetHelixFiberB
from polnet.lrandom import PGenHelixFiber, PGenHelixFiberB
from polnet.polymer import FiberUnitSDimer, MTUnit
from polnet.lio import *
from polnet.utils import *
from polnet.affine import *

# Common settings
VOI_SHAPE = (924, 924, 300) # vx
VOI_OFF = 4 # vx
VOI_VSIZE = 13.68 # A/vx

# NetSAWLC settings
MMER_MRC = './in/emd_5225_13.68A_30.mrc'
MMER_ISO = .543
PMER_L = 1.2 # Number of times wrt MMER diameters
PMER_OCC = 1 # 5 # %
PMER_L_MAX = 3000
PMER_OVER_TOL = 0.01
MMER_OUT = './out/sawlc_mmer.vtp'
NET_OUT = './out/sawlc_net.vtp'
NET_SKEL_OUT = './out/sawlc_net_skel.vtp'
NET_TOMO_OUT = './out/sawlc_net_tomo.mrc'

# NetHelixFiber settings
A_MMER_RAD = 25 # A
A_PMER_OCC = 0.1 # 5 # %
A_MIN_P_LEN = 17.7e4 # 5200e4 # Actin filament 17.7um according to 10.1083/jcb.130.4.909
A_HP_LEN = 720 # A
A_MZ_LEN = 50 # A
A_MZ_LEN_F = 0.2
A_BPROP = 1 / 2
A_MAX_P_BRANCH = 2
A_UNIT_OUT = './out/helix_unit.vtp'
A_UNIT_TOMO_OUT = './out/helix_unit.mrc'
A_NET_OUT = './out/helix_net.vtp'
A_NET_SKEL_OUT = './out/helix_net_skel.vtp'
A_NET_TOMO_OUT = './out/helix_net_tomo.mrc'

# NetHelixFiber for MT settings
MT_MMER_RAD = 40 # A according to http://web.physics.ucsb.edu/~jennyr/research1.html
MT_RAD = 100.5 # A according to http://web.physics.ucsb.edu/~jennyr/research1.html
MT_NUNITS = 13
MT_PMER_OCC = 0.1 # 5 # %
MT_MIN_P_LEN = 5200e4 # MT 5200um according to 10.1083/jcb.130.4.909
MT_HP_LEN = 216000 # A
MT_MZ_LEN = 50 # 80 # A
MT_MZ_LEN_F = 0.4
MT_UNIT_OUT = './out/mt_unit.vtp'
MT_UNIT_TOMO_OUT = './out/mt_unit.mrc'
MT_NET_OUT = './out/mt_net.vtp'
MT_NET_SKEL_OUT = './out/mt_net_skel.vtp'
MT_NET_TOMO_OUT = './out/mt_net_tomo.mrc'


class TestNetSAWLC(TestCase):

    def test_build_network(self):

        # return

        # Generate the VOI
        voi = np.ones(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFF:VOI_SHAPE[0]-VOI_OFF, VOI_OFF:VOI_SHAPE[1]-VOI_OFF, VOI_OFF:VOI_SHAPE[2]-VOI_OFF] = True

        # Polymer parameters
        model = load_mrc(MMER_MRC)
        model_surf = iso_surface(model, MMER_ISO, closed=False, normals=None)
        center = .5 * np.asarray(model.shape, dtype=float)
        # Monomer centering
        model_surf = poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = poly_scale(model_surf, VOI_VSIZE)
        surf_diam = poly_max_distance(model_surf)
        save_vtp(model_surf, MMER_OUT)
        pol_l_generator = PGenHelixFiber()

        # Network generation
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, PMER_L*surf_diam, model_surf, PMER_L_MAX, pol_l_generator, PMER_OCC, PMER_OVER_TOL)
        net_sawlc.build_network()

        # Density tomogram generation
        tomo = np.ones(shape=net_sawlc.get_voi().shape, dtype=np.float32)
        net_sawlc.insert_density_svol(model, tomo, VOI_VSIZE, merge='min')

        # Save the results
        save_vtp(net_sawlc.get_vtp(), NET_OUT)
        save_vtp(net_sawlc.get_skel(), NET_SKEL_OUT)
        write_mrc(tomo, NET_TOMO_OUT, v_size=VOI_VSIZE)


class TestNetHelixFiber(TestCase):

    def test_build_network(self):

        # Generate the VOI
        voi = np.ones(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFF:VOI_SHAPE[0]-VOI_OFF, VOI_OFF:VOI_SHAPE[1]-VOI_OFF, VOI_OFF:VOI_SHAPE[2]-VOI_OFF] = True

        # Fiber unit generation
        funit = FiberUnitSDimer(A_MMER_RAD, VOI_VSIZE)
        model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
        save_vtp(model_surf, A_UNIT_OUT)
        write_mrc(model_svol, A_UNIT_TOMO_OUT)

        # Helix Fiber parameters model
        pol_generator = PGenHelixFiberB()

        # Network generation
        net_helix = NetHelixFiberB(voi, VOI_VSIZE, PMER_L*A_MMER_RAD*2, model_surf, pol_generator, A_PMER_OCC,
                                  A_MIN_P_LEN, A_HP_LEN, A_MZ_LEN, A_MZ_LEN_F, A_BPROP, A_MAX_P_BRANCH, PMER_OVER_TOL)
        net_helix.build_network()

        # Density tomogram generation
        tomo = np.ones(shape=net_helix.get_voi().shape, dtype=np.float32)
        net_helix.insert_density_svol(model_svol, tomo, VOI_VSIZE, merge='min')

        # Save the results
        save_vtp(net_helix.get_vtp(), A_NET_OUT)
        save_vtp(net_helix.get_skel(), A_NET_SKEL_OUT)
        write_mrc(tomo, A_NET_TOMO_OUT, v_size=VOI_VSIZE)

    def test_build_network_mt(self):

        # return

        # Generate the VOI
        voi = np.ones(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFF:VOI_SHAPE[0]-VOI_OFF, VOI_OFF:VOI_SHAPE[1]-VOI_OFF, VOI_OFF:VOI_SHAPE[2]-VOI_OFF] = True

        # Fiber unit generation
        funit = MTUnit(MT_MMER_RAD, MT_RAD, MT_NUNITS, VOI_VSIZE)
        model_svol, model_surf = funit.get_tomo(), funit.get_vtp()
        save_vtp(model_surf, MT_UNIT_OUT)
        write_mrc(model_svol, MT_UNIT_TOMO_OUT)

        # Helix Fiber parameters model
        pol_generator = PGenHelixFiber()

        # Network generation
        net_helix = NetHelixFiber(voi, VOI_VSIZE, PMER_L*MT_MMER_RAD*2, model_surf, pol_generator, MT_PMER_OCC,
                                  MT_MIN_P_LEN, MT_HP_LEN, MT_MZ_LEN, MT_MZ_LEN_F, PMER_OVER_TOL)
        net_helix.build_network()

        # Density tomogram generation
        tomo = np.ones(shape=net_helix.get_voi().shape, dtype=np.float32)
        net_helix.insert_density_svol(model_svol, tomo, VOI_VSIZE, merge='min')

        # Save the results
        save_vtp(net_helix.get_vtp(), MT_NET_OUT)
        save_vtp(net_helix.get_skel(), MT_NET_SKEL_OUT)
        write_mrc(tomo, MT_NET_TOMO_OUT, v_size=VOI_VSIZE)