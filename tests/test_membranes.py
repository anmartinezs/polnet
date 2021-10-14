import math
from unittest import TestCase

import numpy as np

from polnet.membrane import SetMembranes
from polnet.lrandom import EllipGen, SphGen, TorGen
from polnet.lio import *
from polnet.utils import *
from polnet.affine import *

# Common settings
VOI_SHAPE = (400, 400, 150) # (924, 924, 300) # vx
VOI_OFF = 4 # vx
VOI_VSIZE = 13.68 # A/vx

# Generic settings for membranes
MB_OCC = 0.05 # %
MB_THICK_RG = (25, 35) # A
MB_LAYER_S_RG = (1, 1.5) # A
MB_MAX_ECC = .75
MB_OVER_TOL = 1e-8
MB_MIN_RAD = 75 # A

# Set Ellipsoids membrane settings
MB_OUT = './out/ellip_mbs.vtp'
MB_TOMO_OUT = './out/ellip_mbs.mrc'
MB_GTRUTH_OUT = './out/ellip_mbs_gt.mrc'

# Set Spherical membrane settings
MB_SPH_OUT = './out/sph_mbs.vtp'
MB_SPH_TOMO_OUT = './out/sph_mbs.mrc'
MB_SPH_GTRUTH_OUT = './out/sph_mbs_gt.mrc'

# Set Toroidal membrane settings
MB_TOR_OUT = './out/torus_mbs.vtp'
MB_TOR_TOMO_OUT = './out/torus_mbs.mrc'
MB_TOR_GTRUTH_OUT = './out/torus_mbs_gt.mrc'


class TestEllipMembranes(TestCase):

    def test_build_network(self):

        # return

        # Generate the VOI
        voi = np.ones(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFF:VOI_SHAPE[0]-VOI_OFF, VOI_OFF:VOI_SHAPE[1]-VOI_OFF, VOI_OFF:VOI_SHAPE[2]-VOI_OFF] = True

        # Membrane random generation parametrization
        param_rg = (MB_MIN_RAD, math.sqrt(3) * max(VOI_SHAPE) * VOI_VSIZE, MB_MAX_ECC)
        mb_ellip_generator = EllipGen(radius_rg=param_rg[:2], max_ecc=param_rg[2])
        mb_sph_generator = SphGen(radius_rg=(.5*param_rg[0], .5*param_rg[1]))
        mb_tor_generator = TorGen(radius_rg=(.5*param_rg[0], .5*param_rg[1]))

        # Network generation for Ellipsoids
        # print('Generating Ellipsoids...')
        # set_mbs = SetMembranes(voi, VOI_VSIZE, mb_ellip_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
        #                        MB_OVER_TOL)
        # set_mbs.build_set()
        # save_vtp(set_mbs.get_vtp(), MB_OUT)
        # write_mrc(set_mbs.get_tomo(), MB_TOMO_OUT, v_size=VOI_VSIZE, dtype=np.float32)
        # write_mrc(set_mbs.get_gtruth().astype(np.int16), MB_GTRUTH_OUT, v_size=VOI_VSIZE, dtype=np.float32)

        # # Network generation for Spheres
        # print('Generating Spheres...')
        # set_mbs = SetMembranes(voi, VOI_VSIZE, mb_sph_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
        #                        MB_OVER_TOL)
        # set_mbs.build_set()
        # save_vtp(set_mbs.get_vtp(), MB_SPH_OUT)
        # write_mrc(set_mbs.get_tomo(), MB_SPH_TOMO_OUT, v_size=VOI_VSIZE, dtype=np.float32)
        # write_mrc(set_mbs.get_gtruth().astype(np.int16), MB_SPH_GTRUTH_OUT, v_size=VOI_VSIZE, dtype=np.float32)

        # Network generation for Toroids
        print('Generating Toroids...')
        set_mbs = SetMembranes(voi, VOI_VSIZE, mb_tor_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
                               MB_OVER_TOL)
        set_mbs.build_set()
        save_vtp(set_mbs.get_vtp(), MB_TOR_OUT)
        write_mrc(set_mbs.get_tomo(), MB_TOR_TOMO_OUT, v_size=VOI_VSIZE, dtype=np.float32)
        write_mrc(set_mbs.get_gtruth().astype(np.int16), MB_TOR_GTRUTH_OUT, v_size=VOI_VSIZE, dtype=np.float32)

    def test_ellip_surf(self):

        ellipsoid = vtk.vtkParametricEllipsoid()
        ellipsoid.SetXRadius(1)
        ellipsoid.SetYRadius(0.25)
        ellipsoid.SetZRadius(0.75)
        ellipsoidSource = vtk.vtkParametricFunctionSource()
        ellipsoidSource.SetParametricFunction(ellipsoid)
        ellipsoidSource.SetScalarModeToZ()
        ellipsoidSource.Update()

        hold = ellipsoidSource.GetOutput()
        save_vtp(hold, "./out/ellip_surf.vtp")

        distance = vtk.vtkSignedDistance()
        distance.SetInputData(hold)
        distance.SetDimensions(400, 400, 150)
        bounds = hold.GetBounds()
        range = np.asarray((5, 5, 5))
        distance.SetBounds(bounds[0] - range[0] * 0.1, bounds[1] + range[0] * 0.1,
                            bounds[2] - range[1] * 0.1, bounds[3] + range[1] * 0.1,
                            bounds[4] - range[2] * 0.1, bounds[5] + range[2] * 0.1)
        distance.Update()
        hold2 = distance.GetOutput()
        save_vti(hold2, "./out/ellip_dist.vti")

        # ellipsoidMapper = vtk.vtkPolyDataMapper()
        # ellipsoidMapper.SetInputConnection(ellipsoidSource.GetOutputPort())
        # ellipsoidMapper.SetScalarRange(-0.5, 0.5)


