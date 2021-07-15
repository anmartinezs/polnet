from unittest import TestCase

from polnet.network import NetSAWLC
from polnet.lrandom import PGenUniformInRange
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
PMER_OCC = .1 # 1 # 5 # %
PMER_L_RANGE = (1000, 3000)
PMER_OVER_TOL = 0.01
MMER_OUT = './out/sawlc_mmer.vtp'
NET_OUT = './out/sawlc_net.vtp'
NET_SKEL_OUT = './out/sawlc_net_skel.vtp'
NET_TOMO_OUT = './out/sawlc_net_tomo.mrc'


class TestNetSAWLC(TestCase):

    def test_build_network(self):

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
        pol_l_generator = PGenUniformInRange(PMER_L_RANGE[0], PMER_L_RANGE[1])

        # Network generation
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, PMER_L*surf_diam, model_surf, pol_l_generator, PMER_OCC, PMER_OVER_TOL)
        net_sawlc.build_network()

        # Density tomogram generation
        tomo = np.zeros(shape=net_sawlc.get_voi().shape, dtype=np.float32)
        net_sawlc.insert_density_svol(model, tomo, VOI_VSIZE, merge='max')

        # Save the results
        save_vtp(net_sawlc.get_vtp(), NET_OUT)
        save_vtp(net_sawlc.get_skel(), NET_SKEL_OUT)
        write_mrc(tomo, NET_TOMO_OUT, v_size=VOI_VSIZE)
