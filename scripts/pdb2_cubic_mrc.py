"""
This script has the purpose of generate cubic density maps from pdb atomic models
It is based on e2pdbmrc.py script, so EMAN is required
"""

__author__ = 'Antonio Martinez-Sanchez'

import os
import time
import subprocess

import numpy as np
from polnet import lio


ROOT_DIR = '/fs/pool/pool-lucic2/antonio/polnet/riboprot/templates/pdbs' # '/fs/pool/pool-lucic2/antonio/polnet/riboprot/templates/mb_pdbs'
CMD_PDB2MRC = 'e2pdb2mrc.py'

# Input list of PDBs to process
in_pdb_l = ['3ipm.pdb', '1bxn.pdb', '1qvr.pdb', '1s3x.pdb', '1u6g.pdb', '2cg9.pdb', '3cf3.pdb', '3d2f.pdb', '3gl1.pdb',
            '3h84.pdb', '3qm1.pdb', '4cr2.pdb', '4v94.pdb', '5mrc.pdb', '6utj.pdb', '4v4r_30S.pdb', '4v4r_50S.pdb',
            '3j9i.pdb', '4v4r.pdb']
           # ['5wek.pdb', '4pe5.pdb', '5ide.pdb', '5gjv.pdb', '5kxi.pdb', '5tj6.pdb', '5tqq.pdb', '5vai.pdb']

# Output directory
out_dir = '../mrcs_10A' # '../mb_mrcs_10A'

# Input parameters
v_size = 10 # 2.2 # A/voxel
res = 30 # 8.8 # resolution in A (low-pass Gaussian 1/res)
offset = 20 # 60 # voxels


# Loop for PDBs
for in_pdb in in_pdb_l:

    print('Processing PDB:', in_pdb)

    # First command: guess minimal dimension
    cmd_1 = [CMD_PDB2MRC,]
    cmd_1 += ['--apix', str(v_size)]
    cmd_1 += [ROOT_DIR + '/' + in_pdb]
    hold_out_mrc = ROOT_DIR + '/' + out_dir + '/' + os.path.splitext(os.path.split(in_pdb)[1])[0] + '_hold.mrc'
    cmd_1 += [hold_out_mrc]
    print('\t-Running commad:', ''.join(cmd_1))
    subprocess.run(cmd_1)

    # Computing the cubic bounding box
    hold_mrc = lio.load_mrc(hold_out_mrc)
    max_dim = str(np.asarray(hold_mrc.shape).max() + 2*offset)
    cubic_shape = [max_dim, max_dim, max_dim]

    # Second command: generate the cubic model
    cmd_2 = [CMD_PDB2MRC, ]
    cmd_2 += ['--box', max_dim + ',' + max_dim + ',' + max_dim]
    cmd_2 += ['--apix', str(v_size)]
    cmd_2 += ['--res', str(res)]
    cmd_2 += [ROOT_DIR + '/' + in_pdb]
    out_mrc = ROOT_DIR + '/' + out_dir + '/' + os.path.splitext(os.path.split(in_pdb)[1])[0] + '.mrc'
    cmd_2 += [out_mrc]
    print('\t-Running commad:', cmd_2)
    subprocess.run(cmd_2)


print('Successfully terminated. (' + time.strftime("%c") + ')')