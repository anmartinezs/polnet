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
import mrcfile
import numpy as np
import vtk
from vtkmodules.util.numpy_support import vtk_to_numpy
import pandas as pd

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

    def load_mrc(self, filename):
        with mrcfile.open(filename) as mrc:
            return mrc.data

    def load_vtp(self, filename):
        reader = vtk.vtkXMLPolyDataReader()
        reader.SetFileName(filename)
        reader.Update()
        polydata = reader.GetOutput()
        if polydata.GetNumberOfPoints() == 0:
            print(f"Warning: VTP file {filename} contains no points.")
            return np.array([])
        points = vtk_to_numpy(polydata.GetPoints().GetData())
        return points

    def compare_arrays(self, arr1, arr2):
        return np.allclose(arr1, arr2, rtol=1e-5, atol=1e-8)

    def compare_csv(self, file1, file2):
        df1 = pd.read_csv(file1, delimiter='\t')
        df2 = pd.read_csv(file2, delimiter='\t')
        
        # Columns to ignore in the comparison
        ignore_columns = ['Density', 'Micrographs', 'PolyData', 'Tomo3D']
        
        columns_to_compare = [col for col in df1.columns if col not in ignore_columns]
        
        if not df1[columns_to_compare].equals(df2[columns_to_compare]):
            print(f"Differences in CSV files:")
            for column in columns_to_compare:
                if not df1[column].equals(df2[column]):
                    print(f"  Column '{column}' differs:")
                    diff = ~(df1[column] == df2[column])
                    print(f"    File 1:\n{df1.loc[diff, column]}")
                    print(f"    File 2:\n{df2.loc[diff, column]}")
            return False
        return True

    def test_repeatability(self):
        seed = 17
        
        print("Running first simulation...")
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/run1", random_seed=seed)
        print("First simulation completed.")
        
        print("Running second simulation...")
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/run2", random_seed=seed)
        print("Second simulation completed.")
        
        key_files = [
            "tomos/tomo_den_0.mrc",
            "tomos/tomo_lbls_0.mrc",
            "tomos/poly_den_0.vtp",
            "tomos/poly_skel_0.vtp",
            "tomos_motif_list.csv"
        ]
        
        results = {}
        
        for file in key_files:
            file1 = os.path.join(self.out_dir, "run1", file)
            file2 = os.path.join(self.out_dir, "run2", file)
            
            print(f"Comparing {file}...")
            
            if file.endswith('.mrc'):
                arr1 = self.load_mrc(file1)
                arr2 = self.load_mrc(file2)
                results[file] = self.compare_arrays(arr1, arr2)
            elif file.endswith('.vtp'):
                points1 = self.load_vtp(file1)
                points2 = self.load_vtp(file2)
                if len(points1) == 0 and len(points2) == 0:
                    print(f"Both VTP files {file} are empty. Skipping comparison.")
                    results[file] = "Skipped (Empty)"
                else:
                    results[file] = self.compare_arrays(points1, points2)
            elif file.endswith('.csv'):
                results[file] = self.compare_csv(file1, file2)
            
            print(f"Comparison of {file} completed.")
        
        print("\nComparison Results:")
        for file, result in results.items():
            print(f"{file}: {'Passed' if result == True else 'Failed' if result == False else result}")
        
        self.assertTrue(all(result == True for result in results.values() if result != "Skipped (Empty)"), 
                        "Not all files were identical between runs")

    def test_seed_sensitivity(self):
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/seed1", random_seed=1)
        all_features2(**self.light_params, OUT_DIR=f"{self.out_dir}/seed2", random_seed=2)
        
        key_file = "tomos/tomo_den_0.mrc"
        arr1 = self.load_mrc(os.path.join(self.out_dir, "seed1", key_file))
        arr2 = self.load_mrc(os.path.join(self.out_dir, "seed2", key_file))
        
        self.assertFalse(self.compare_arrays(arr1, arr2), "Output should be different with different seeds")

if __name__ == '__main__':
    unittest.main()