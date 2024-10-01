import unittest
import numpy as np
from polnet.membrane import SetMembranes
from polnet.lrandom import EllipGen
from polnet.lio import write_mrc, save_vtp
import random
from gui.core.all_features2 import all_features2
import os
import shutil
import hashlib

class TestRepeatability(unittest.TestCase):

    def setUp(self):
        self.voi_shape = (100, 100, 50)
        self.voi_vsize = 10
        self.mb_occ = 0.002
        self.mb_thick_rg = (25, 35)
        self.mb_layer_s_rg = (5, 10)
        self.mb_max_ecc = 0.75
        self.mb_over_tol = 0
        self.mb_min_rad = 75
        self.random_seed = 42

    def generate_membranes(self):
        random.seed(self.random_seed)
        np.random.seed(self.random_seed)

        voi = np.ones(shape=self.voi_shape, dtype=bool)
        param_rg = (self.mb_min_rad, np.sqrt(3) * max(self.voi_shape) * self.voi_vsize, self.mb_max_ecc)
        mb_ellip_generator = EllipGen(radius_rg=param_rg[:2], max_ecc=param_rg[2])

        set_mbs = SetMembranes(voi, self.voi_vsize, mb_ellip_generator, param_rg, self.mb_thick_rg, 
                               self.mb_layer_s_rg, self.mb_occ, self.mb_over_tol)
        
        set_mbs.build_set(verbosity=False)

        tomo = set_mbs.get_tomo()
        vtp = set_mbs.get_vtp()

        return tomo, vtp

    def test_repeatability(self):
        tomo1, vtp1 = self.generate_membranes()
        tomo2, vtp2 = self.generate_membranes()

        np.testing.assert_array_equal(tomo1, tomo2, "Tomograms are not identical")
        
        self.assertEqual(vtp1.GetNumberOfPoints(), vtp2.GetNumberOfPoints(), 
                         "VTPs have different number of points")
        self.assertEqual(vtp1.GetNumberOfCells(), vtp2.GetNumberOfCells(), 
                         "VTPs have different number of cells")

        for i in range(vtp1.GetNumberOfPoints()):
            point1 = vtp1.GetPoint(i)
            point2 = vtp2.GetPoint(i)
            np.testing.assert_array_almost_equal(point1, point2, 
                                                 err_msg=f"Points at index {i} are not equal")

class TestAllFeatures2Repeatability(unittest.TestCase):
    def setUp(self):
        self.out_dir = "./test_output_light"
        if not os.path.exists(self.out_dir):
            os.makedirs(self.out_dir)
        
        # Get the absolute path to the data directory
        self.data_dir = os.path.abspath(os.path.join(os.path.dirname(__file__), '..', 'data'))
        
        self.light_params = {
            "NTOMOS": 1,
            "VOI_SHAPE": (100, 100, 50),
            "VOI_OFFS": ((2, 98), (2, 98), (2, 48)),
            "VOI_VSIZE": 10,
            "MMER_TRIES": 10,
            "PMER_TRIES": 50,
            "MEMBRANES_LIST": [os.path.join(self.data_dir, "in_mbs", "sphere.mbs")],
            "HELIX_LIST": [os.path.join(self.data_dir, "in_helix", "actin.hns")],
            "PROTEINS_LIST": [os.path.join(self.data_dir, "in_10A", "4v7r_10A.pns")],
            "MB_PROTEINS_LIST": [],
            "SURF_DEC": 0.9,
            "TILT_ANGS": (-30, 31, 5),
            "DETECTOR_SNR": 1,
            "MALIGN_MN": 1,
            "MALIGN_MX": 1.5,
            "MALIGN_SG": 0.2
        }

    def tearDown(self):
        if os.path.exists(self.out_dir):
            shutil.rmtree(self.out_dir)

    def hash_file(self, filename):
        with open(filename, "rb") as f:
            return hashlib.md5(f.read()).hexdigest()

    def test_repeatability(self):
        # Run all_features2 twice with the same random seed
        seed = 42
        
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/run1", random_seed=seed)
        
        # Check only specific output files instead of entire directory
        key_files = [
            "tomos/tomo_den_0.mrc",
            "tomos/tomo_lbls_0.mrc",
            "tomos/poly_den_0.vtp",
            "tomos/poly_skel_0.vtp",
            "tomos_motif_list.csv"
        ]
        
        hashes1 = {file: self.hash_file(os.path.join(self.out_dir, "run1", file)) for file in key_files}
        
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/run2", random_seed=seed)
        hashes2 = {file: self.hash_file(os.path.join(self.out_dir, "run2", file)) for file in key_files}
        
        self.assertEqual(hashes1, hashes2, "Key output files are not identical between runs")

    def test_seed_sensitivity(self):
        # Run with different seeds and ensure key outputs are different
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/seed1", random_seed=1)
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/seed2", random_seed=2)
        
        key_file = "tomos/tomo_den_0.mrc"
        hash1 = self.hash_file(os.path.join(self.out_dir, "seed1", key_file))
        hash2 = self.hash_file(os.path.join(self.out_dir, "seed2", key_file))
        
        self.assertNotEqual(hash1, hash2, "Output should be different with different seeds")


if __name__ == '__main__':
    unittest.main()