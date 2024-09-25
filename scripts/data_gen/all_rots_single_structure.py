"""
This script generates subvolumes sampling uniformly the whole space of rotation of single structure
"""

import random
import shutil

import numpy as np

from polnet import lio, affine, tem

ROOT_DIR = r'C:\Users\amart\workspace\data'
IN_TEMPL = ROOT_DIR + r'\synth_data\templates\mrcs_10A\3j9i.mrc'
N_ROTS = 150
TILT_ANGS = range(-60, 61, 3) # np.arange(-60, 60, 3) # at MPI-B IMOD only works for ranges
DETECTOR_SNR = [0.1, 0.2] # [.15, .25]
MALIGN_MN = 1
MALIGN_MX = 1.5
MALIGN_SG = 0.2
VOI_VSIZE = 10 # 2.2 # A/vx
TEM_DIR = ROOT_DIR + r'\tda\tem'
TOMOS_DIR = ROOT_DIR + r'\tda\snr_1-2_tilt_60\3j9i'

templ = lio.load_mrc(IN_TEMPL)
if N_ROTS <= 0:
    rots = np.asarray([[.0, .0, .0, .0]])
else:
    rots = affine.uniform_sampling_so3(N_ROTS)

for i, rot in enumerate(rots):

    print('\t-Processing particle ' + str(i) + '/' + str(N_ROTS))

    # Rotations
    if N_ROTS <= 0:
        templ_r = templ
    else:
        templ_r = affine.tomo_rotate(templ, rot)

    # TEM
    temic = tem.TEM(TEM_DIR)
    temic.gen_tilt_series_imod(templ_r, TILT_ANGS, ax='Y')
    temic.add_mics_misalignment(MALIGN_MN, MALIGN_MX, MALIGN_SG)
    if DETECTOR_SNR is not None:
        if hasattr(DETECTOR_SNR, '__len__'):
            if len(DETECTOR_SNR) >= 2:
                snr = round((DETECTOR_SNR[1] - DETECTOR_SNR[0]) * random.random() + DETECTOR_SNR[0], 2)
            else:
                snr = DETECTOR_SNR[0]
        else:
            snr = DETECTOR_SNR
        temic.add_detector_noise(snr)
    temic.invert_mics_den()
    temic.set_header(data='mics', p_size=(VOI_VSIZE, VOI_VSIZE, VOI_VSIZE))
    temic.recon3D_imod()
    temic.set_header(data='rec3d', p_size=(VOI_VSIZE, VOI_VSIZE, VOI_VSIZE), origin=(0, 0, 0))
    out_mics, out_tomo_rec = TOMOS_DIR + '/tomo_mics_' + str(i) + '.mrc', TOMOS_DIR + '/tomo_rec_' + str(i) + '.mrc'
    out_tomo_vol = TOMOS_DIR + '/tomo_vol_' + str(i) + '.mrc'
    shutil.copyfile(TEM_DIR + '/out_micrographs.mrc', out_mics)
    shutil.copyfile(TEM_DIR + '/out_rec3d.mrc', out_tomo_rec)
    shutil.copyfile(TEM_DIR + '/in_vol.mrc', out_tomo_vol)
