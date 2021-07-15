from unittest import TestCase

import vtk
from polnet.polymer import Monomer
from polnet.lio import *
from polnet.utils import *
from polnet.affine import *

# MMER SETTINGS
VOI_VSIZE = 13.68 # A/vx
MMER_MRC = './in/emd_5225_13.68A_30.mrc'
MMER_ISO = .543
MMER_OUT_VTK = './out/mmer.vtp'
MMER_OUT_VTK_2 = './out/mmer_2.vtp'
MMER_OUT_MRC = './out/mmer.mrc'


class TestMonomer(TestCase):

    def test_insert_density_svol(self):

        # Generate the surface
        model = load_mrc(MMER_MRC)
        model_surf = iso_surface(model, MMER_ISO, closed=False, normals=None)
        center = .5 * np.asarray(model.shape, dtype=float)
        # Monomer centering
        model_surf = poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = poly_scale(model_surf, VOI_VSIZE)
        m_diam = poly_max_distance(model_surf)

        # Create the mmer
        mmer = Monomer(model_surf, m_diam)

        # Apply rigid transformations
        q = gen_rand_unit_quaternion()
        # q = np.asarray((0, 0, 0, 1.))
        q /= vector_module(q)
        mmer.rotate_q(q)
        mmer.translate(center * VOI_VSIZE)
        # m_mrc = tomo_rotate_sitk(model, q, active=True)
        m_mrc = tomo_rotate(model, q)

        # Store the output
        save_vtp(mmer.get_vtp(), MMER_OUT_VTK)
        write_mrc(m_mrc, MMER_OUT_MRC, v_size=VOI_VSIZE)
        model_surf_2 = iso_surface(m_mrc, MMER_ISO, closed=False, normals=None)
        model_surf_2 = poly_scale(model_surf_2, VOI_VSIZE)
        save_vtp(model_surf_2, MMER_OUT_VTK_2)