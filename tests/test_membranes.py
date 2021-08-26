from unittest import TestCase

from polnet.membrane import SetEllipMembranes
from polnet.lrandom import EllipGen
from polnet.lio import *
from polnet.utils import *
from polnet.affine import *

# Common settings
VOI_SHAPE = (924, 924, 300) # vx
VOI_OFF = 4 # vx
VOI_VSIZE = 13.68 # A/vx

# Set Ellipsoids membrane settings
MB_OCC = 1 # %
MB_PARAM_RG = (200, 20e4) # A
MB_THICK_RG = (50, 1000) # A
MB_LAYER_S_RG = (5, 15) # A
MB_OVER_TOL = 0.01
MB_OUT = './out/ellip_mbs.vtp'
MB_TOMO_OUT = './out/ellip_mbs.mrc'
MB_GTRUTH_OUT = './out/ellip_mbs_gt.mrc'


class TestEllipMembranes(TestCase):

    def test_build_network(self):

        # return

        # Generate the VOI
        voi = np.ones(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFF:VOI_SHAPE[0]-VOI_OFF, VOI_OFF:VOI_SHAPE[1]-VOI_OFF, VOI_OFF:VOI_SHAPE[2]-VOI_OFF] = True

        # Membrane random generation parametrization
        mb_ellip_generator = EllipGen()

        # Network generation
        set_mbs = SetEllipMembranes(voi, VOI_VSIZE, mb_ellip_generator, MB_PARAM_RG, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC, MB_OVER_TOL)
        set_mbs.build_set()

        # Save the results
        save_vtp(set_mbs.get_vtp(), MB_OUT)
        write_mrc(set_mbs.get_tomo(), MB_TOMO_OUT, v_size=VOI_VSIZE)
        write_mrc(set_mbs.get_gtruth(), MB_GTRUTH_OUT, v_size=VOI_VSIZE)


