"""
This script has the purpose of generate cubic density maps to work as references to align membrane bound proteins
"""

__author__ = 'Antonio Martinez-Sanchez'

import os
import time
import math

import numpy as np
from polnet import lio, affine


ROOT_DIR = '/fs/pool/pool-lucic2/antonio/polnet/riboprot/templates/mb_mrcs_10A'

# Input list of MRCs to process
in_mrc_l = ['5wek.mrc', '4pe5.mrc', '5ide.mrc', '5gjv.mrc', '5kxi.mrc', '5tj6.mrc', '5tqq.mrc', '5vai.mrc']
mb_size = 50 # A
zaxis_rad = 5 # A
vsize = 10 # 2.2 # A/vx


# Output directory
out_dir = '../mb_mrcs_10A' # '../mb_mrcs_2.2A'


# Loop for PDBs
mb_size_vx, zaxis_rad_vx = mb_size / vsize, zaxis_rad / vsize
for in_mrc in in_mrc_l:

    print('Processing MRC:', ROOT_DIR + '/' + in_mrc)
    tomo = lio.load_mrc(ROOT_DIR + '/' + in_mrc)

    # Constructing the reference subvolume
    center = .5 * (np.asarray(tomo.shape, dtype=np.float32) - 1)
    X, Y, Z = np.meshgrid(np.arange(tomo.shape[0]), np.arange(tomo.shape[1]), np.arange(tomo.shape[2]), indexing='ij')
    X, Y, Z = (X - center[0]).astype(np.float32), (Y - center[1]).astype(np.float32), (Z - center[2]).astype(np.float32)
    mask_mb = (Z <= 0) * (Z >= -mb_size_vx)
    zaxis_mask = X*X + Y*Y <= zaxis_rad
    ref_svol = np.zeros(shape=tomo.shape, dtype=np.int16)
    ref_svol[mask_mb] = 1
    ref_svol[zaxis_mask] = 1

    # First command: guess minimal dimension
    hold_out_mrc = ROOT_DIR + '/' + out_dir + '/' + os.path.splitext(os.path.split(in_mrc)[1])[0] + '_mb_ref.mrc'
    print('\t-Saving the output file in:', hold_out_mrc)
    lio.write_mrc(ref_svol, hold_out_mrc, v_size=vsize)


print('Successfully terminated. (' + time.strftime("%c") + ')')