from unittest import TestCase

from polnet.lrandom import PGenHelixFiber
from polnet.polymer import Monomer, FiberUnitSDimer, HelixFiber
from polnet.network import NetSAWLC
from polnet.lio import *
from polnet.utils import *
from polnet.affine import *
from polnet.poly import poly_diam

# MMER SETTINGS
VOI_VSIZE = 10 # A/vx
VOI_SHAPE = (125, 125, 125)
VOI_OFFS =  ((4,121), (4,121), (4,121))
VOI_OUT = './out/voi.mrc'
MMER_MRC = './in/5mrc.mrc'
MMER_ISO = .35
MMER_OUT_VTK = './out/mmer.vtp'
MMER_OUT_VTK_2 = './out/mmer_2.vtp'
MMER_OUT_MRC = './out/mmer.mrc'

# Helix polymer settings
HLIX_MMER_RAD = 25
HLIX_PMER_L = 1.2
HLIX_P_LEN = 1e3
HLIX_HP_LEN = 720
HLIX_MZ_LEN = 50
HLIX_MZ_LEN_F = 0.2
HLIX_MAX_LEN = 7500
HLIX_OUT_VTK = './out/helix_pmmer.vtp'
HLIX_OUT_SKEL_VTK = './out/helix_pmmer_skel.vtp'

# SALWC polymer setting
SAWLC_PMER_L = 1.2
SAWLC_PMER_L_MAX = 1000
SAWLC_OCC = 3
SAWLC_OVER = 0
SAWLC_OUT_VTK = './out/sawlc_pmmer.vtp'
SAWLC_OUT_SKEL_VTK = './out/sawlc_pmmer_skel.vtp'


class TestMonomer(TestCase):

    def test_insert_density_svol(self):

        return

        # Generate the surface
        model = load_mrc(MMER_MRC)
        model_surf = iso_surface(model, MMER_ISO, closed=False, normals=None)
        center = .5 * np.asarray(model.shape, dtype=float)
        # Monomer centering
        model_surf = poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = poly_scale(model_surf, VOI_VSIZE)
        m_diam = poly_diam(model_surf)

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

    def test_helix_polymer(self):

        return

        # Fiber unit generation
        funit = FiberUnitSDimer(HLIX_MMER_RAD, VOI_VSIZE)

        # Network generation
        polymer = HelixFiber(HLIX_PMER_L, funit.get_vtp(), HLIX_P_LEN, HLIX_HP_LEN, HLIX_MZ_LEN,
                             HLIX_MZ_LEN_F, p0=(0, 0, 0), vz=(0, 0, 1), rot_rand=False)
        monomer_data = polymer.gen_new_monomer(100, voi=None, v_size=VOI_VSIZE)
        while polymer.get_total_len() < HLIX_MAX_LEN:
            polymer.add_monomer(monomer_data[0], monomer_data[1], monomer_data[2], monomer_data[3])
            monomer_data = polymer.gen_new_monomer(monomer_data, voi=None, v_size=VOI_VSIZE)

        # Save the results
        save_vtp(polymer.get_vtp(), HLIX_OUT_VTK)
        save_vtp(polymer.get_skel(), HLIX_OUT_SKEL_VTK)

    def test_SAWLC_polymer(self):

        # return

        # VOI
        voi = np.zeros(shape=VOI_SHAPE, dtype=bool)
        voi[VOI_OFFS[0][0]:VOI_OFFS[0][1], VOI_OFFS[1][0]:VOI_OFFS[1][1], VOI_OFFS[2][0]:VOI_OFFS[2][1]] = True

        # Polymer parameters
        model = load_mrc(MMER_MRC)
        model = lin_map(model, lb=0, ub=1)
        model = vol_cube(model)
        model_mask = model < MMER_ISO
        model[model_mask] = 0
        model_surf = iso_surface(model, MMER_ISO, closed=False, normals=None)
        center = .5 * np.asarray(model.shape, dtype=float)
        # Monomer centering
        model_surf = poly_translate(model_surf, -center)
        # Voxel resolution scaling
        model_surf = poly_scale(model_surf, VOI_VSIZE)
        surf_diam = poly_diam(model_surf) * SAWLC_PMER_L

        # Network generation (Clustered)
        pol_l_generator = PGenHelixFiber()
        net_sawlc = NetSAWLC(voi, VOI_VSIZE, SAWLC_PMER_L * surf_diam, model_surf, SAWLC_PMER_L_MAX,
                             pol_l_generator, SAWLC_OCC, SAWLC_OVER, poly=None, svol= model_mask)
        net_sawlc.build_network()

        # Save the results
        save_vtp(net_sawlc.get_vtp(), SAWLC_OUT_VTK)
        save_vtp(net_sawlc.get_skel(), SAWLC_OUT_SKEL_VTK)
        write_mrc(voi.astype(np.float32), VOI_OUT, v_size=VOI_VSIZE)

