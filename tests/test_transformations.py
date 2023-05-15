"""
Test for checking if the transformations of a CSV motif list are correct and agrees with the density and label maps
"""
import os.path

import numpy as np
import pandas as pd

from unittest import TestCase

from polnet import lio, affine, utils

# Common settings
IN_CSV = '/media/martinez/Sistema/Users/Antonio/workspace/data/polnet/out_all_tomos_5-6/tomos_motif_list.csv'
IN_DEN_TOMO = '/media/martinez/Sistema/Users/Antonio/workspace/data/polnet/out_all_tomos_5-6/tomos/tomo_den_0.mrc'
IN_LABEL_CODE = 11
IN_STRUCT_MODEL = '/media/martinez/Sistema/Users/Antonio/workspace/data/polnet/templates/mrcs_10A/5mrc.mrc'
VOI_VSIZE = 10 # A/vx

# Paths to tomograms outputs
STOMO_TOMO_OUT = './out/density_postion_rotations.mrc'


class TestTransformations(TestCase):

    def test_positions_rotations(self):

        tomo_model = lio.load_mrc(IN_STRUCT_MODEL)
        den_tomo = lio.load_mrc(IN_DEN_TOMO, mmap=True)
        tomo = np.zeros(shape=den_tomo.shape, dtype=tomo_model.dtype)
        den_file_name = os.path.split(IN_DEN_TOMO)[1]

        # Load the motif list and filter the structure of the specific model
        df = pd.read_csv(IN_CSV, delimiter='\t')
        df_filt = df[df['Density'].str.contains(den_file_name) & (df['Label'] == IN_LABEL_CODE)]

        # Loop for processing any entry of the model
        for row in df_filt.iterrows():

            # Get positions and rotations
            center = np.asarray((row[1]['X'], row[1]['Y'], row[1]['Z']))
            q_rots = np.asarray((row[1]['Q1'], row[1]['Q2'], row[1]['Q3'], row[1]['Q4']))

            # Applying the transformations
            tomo_rot = affine.tomo_rotate(tomo_model, q_rots)
            utils.insert_svol_tomo(tomo_rot, tomo, center / VOI_VSIZE, merge='max')

        # Storing the output tomogram
        lio.write_mrc(tomo, STOMO_TOMO_OUT, v_size=VOI_VSIZE)

