from unittest import TestCase

import numpy as np
from polnet.tem import TEM
from polnet import lio

TEM_DIR = "./out/tem_dir"
IN_VOL = "./out/stomo_density.mrc"  # './out/mt_net_tomo.mrc'
TILT_ANGS = np.arange(-60, 60, 2)
DETECTOR_SNR = 0.3


class TestTEM(TestCase):

    def test_tem(self):

        # Simulates a 3D TEM process (currently no MTF or CTF determination are applied)
        tem = TEM(TEM_DIR)
        vol = lio.load_mrc(IN_VOL)
        tem.gen_tilt_series_imod(vol, TILT_ANGS)
        tem.add_detector_noise(DETECTOR_SNR)
        tem.invert_mics_den()
        tem.recon3D_imod()
