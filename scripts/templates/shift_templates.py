"""
This script shifts all input templates (shifting done in real space by voxels)
"""

__author__ = 'Antonio Martinez-Sanchez'

import os
import time
import math

import numpy as np
from polnet import lio, affine


ROOT_DIR = '/fs/pool/pool-lucic2/antonio/polnet/riboprot/templates/mb_mrcs_10A'

# Input list of MRCs to process
in_mrc_l = ['5wek_mbz_align.mrc', '4pe5_mbz_align.mrc', '5ide_mbz_align.mrc', '5gjv_mbz_align.mrc',
            '5kxi_mbz_align.mrc', '5tj6_mbz_align.mrc', '5tqq_mbz_align.mrc', '5vai_mbz_align.mrc']
shift_x = 0 # vx
shift_y = 0 # vx
shift_z = 2 # vx
vsize = 10 # A/vx


# Output directory
out_dir = '../mb_mrcs_10A'


# Loop for PDBs
for in_mrc in in_mrc_l:

    print('Processing MRC:', ROOT_DIR + '/' + in_mrc)
    tomo = lio.load_mrc(ROOT_DIR + '/' + in_mrc)

    # Shifting
    if (shift_x is not None) and (shift_x != 0):
        tomo = np.roll(tomo, shift_x, axis=0)
    if (shift_y is not None) and (shift_y != 0):
        tomo = np.roll(tomo, shift_y, axis=1)
    if (shift_z is not None) and (shift_z != 0):
        tomo = np.roll(tomo, shift_z, axis=2)

    # First command: guess minimal dimension
    hold_out_mrc = ROOT_DIR + '/' + out_dir + '/' + os.path.splitext(os.path.split(in_mrc)[1])[0] + '_shift.mrc'
    print('\t-Saving the output file in:', hold_out_mrc)
    lio.write_mrc(tomo, hold_out_mrc, v_size=vsize)


print('Successfully terminated. (' + time.strftime("%c") + ')')