import math
from unittest import TestCase

import numpy as np

from polnet.membrane import SetMembranes
from polnet.network import NetSAWLC, PGenHelixFiber
from polnet.lrandom import EllipGen, SphGen, TorGen
from polnet.polymer import MB_DOMAIN_FIELD_STR
from polnet.poly import poly_reverse_normals, add_sfield_to_poly
from polnet.lio import *
from polnet.utils import *
from polnet.affine import *

# Common settings
VOI_SHAPE = (400, 400, 150) # (924, 924, 300) # vx
VOI_OFF = 4 # vx
VOI_VSIZE = 13.68 # A/vx

# Generic settings for membranes
MB_OCC = 0.002 # %
MB_THICK_RG = (25, 35) # A
MB_LAYER_S_RG = (5, 10) # A
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

# Membrane bound protein
MMER_MRC = './in/pdb_5wek_13.68A_uncentered.mrc'
MMER_ISO = .046
MMER_CENTER = [13, 14.5, 11] # [57, 76, 36] # voxels
MB_Z_HEIGHT = 16 # voxels
MMER_OUT_VTK = './out/mb_mmer.vtp'
MMER_OUT_VTK_2 = './out/mb_mmer_2.vtp'
MMER_OUT_MRC = './out/mb_mmer.mrc'
PMER_L = 1.2 # Number of times wrt MMER diameters
PMER_OCC = 0.05 # 5 # %
PMER_L_MAX = 3000
PMER_OVER_TOL = 0.01
PMER_REVERSE_NORMALS = True
MMER_OUT = './out/mb_sawlc_mmer.vtp'
NET_OUT = './out/mb_sawlc_net.vtp'
NET_SKEL_OUT = './out/mb_sawlc_net_skel.vtp'
NET_TOMO_OUT = './out/mb_sawlc_net_tomo.mrc'
MB_POLY_OUT = './out/mb_poly.vtp'
NET_VOI_OUT = './out/mb_sawlc_voi.mrc'


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
        print('Generating Ellipsoids...')
        set_mbs = SetMembranes(voi, VOI_VSIZE, mb_ellip_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
                               MB_OVER_TOL)
        set_mbs.build_set(verbosity=True)
        save_vtp(set_mbs.get_vtp(), MB_OUT)
        write_mrc(set_mbs.get_tomo(), MB_TOMO_OUT, v_size=VOI_VSIZE, dtype=np.float32)
        write_mrc(set_mbs.get_gtruth().astype(np.int16), MB_GTRUTH_OUT, v_size=VOI_VSIZE, dtype=np.float32)

        # Insert membrane bound densities in a Polymer
        # Polymer parameters
        model = load_mrc(MMER_MRC)
        model_mask = model < MMER_ISO
        model[model_mask] = 0
        model_surf = iso_surface(model, MMER_ISO, closed=False, normals=None)
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
        surf_diam = poly_max_distance(model_surf)
        save_vtp(model_surf, MMER_OUT)
        pol_l_generator = PGenHelixFiber()
        # Network generation
        mb_poly = set_mbs.get_vtp()
        if PMER_REVERSE_NORMALS:
            mb_poly = poly_reverse_normals(mb_poly)
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, PMER_L * surf_diam, model_surf, PMER_L_MAX, pol_l_generator, PMER_OCC,
                             PMER_OVER_TOL, poly=mb_poly, svol=model < MMER_ISO)
        net_sawlc.build_network()
        voi = net_sawlc.get_voi()
        net_sawlc.insert_density_svol(model_mask, voi, VOI_VSIZE, merge='min', off_svol=off)
        net_sawlc.set_voi(voi)
        # Density tomogram generation
        tomo = np.zeros(shape=net_sawlc.get_voi().shape, dtype=np.float32)
        net_sawlc.insert_density_svol(model, tomo, VOI_VSIZE, merge='max', off_svol=off)
        # Save the results
        save_vtp(net_sawlc.get_vtp(), NET_OUT)
        save_vtp(net_sawlc.get_skel(), NET_SKEL_OUT)
        write_mrc(tomo, NET_TOMO_OUT, v_size=VOI_VSIZE)
        write_mrc(net_sawlc.get_voi().astype(np.float32), NET_VOI_OUT, v_size=VOI_VSIZE)

        # Network generation for Spheres
        print('Generating Spheres...')
        set_mbs = SetMembranes(voi, VOI_VSIZE, mb_sph_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
                               MB_OVER_TOL)
        set_mbs.build_set(verbosity=True)
        save_vtp(set_mbs.get_vtp(), MB_SPH_OUT)
        save_vtp(mb_poly, MB_POLY_OUT)
        write_mrc(set_mbs.get_tomo(), MB_SPH_TOMO_OUT, v_size=VOI_VSIZE, dtype=np.float32)
        write_mrc(set_mbs.get_gtruth().astype(np.int16), MB_SPH_GTRUTH_OUT, v_size=VOI_VSIZE, dtype=np.float32)

        # Network generation for Toroids
        print('Generating Toroids...')
        set_mbs = SetMembranes(voi, VOI_VSIZE, mb_tor_generator, param_rg, MB_THICK_RG, MB_LAYER_S_RG, MB_OCC,
                               MB_OVER_TOL)
        set_mbs.build_set(verbosity=True)
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


