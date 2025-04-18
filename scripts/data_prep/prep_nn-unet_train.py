"""
This script prepares a dataset generated by PolNet to serve as input data to training nn-UNet for segmentation
    Input:
        - A directory with the output of a PolNet dataset (see script data_gen/all_features.py)
        - Labels dictionary to match each structural with its label, labels not included in this dictionary are
          considered background
    Ouput:
        - A new directory with the structure required by nn-UNet for training a model for semantic segmentation
"""

__author__ = "Antonio Martinez-Sanchez"

import os
import shutil
import nrrd
import json
import time

import numpy as np
import pandas

from pathlib import Path
from polnet import lio


ROOT_DIR = Path(__file__)
ROOT_DIR = ROOT_DIR.joinpath(Path("../../../data")).resolve()
in_csv = ROOT_DIR.joinpath(
    Path("data_generated/all_v11/tomos_motif_list.csv")
).resolve()
out_dir = os.getenv("nnUNet_raw") or os.path.join(ROOT_DIR + "/data_prepared")
dataset_id = "099"
dataset_suffix = "cetnet"
fg_labels = {
    "membrane": (1, 2, 3),
    "microtuble": (4,),
    "actin": (5,),
    "ribo": (6, 11, 12),
    "cprots": tuple(np.arange(7, 11).tolist() + np.arange(13, 26).tolist()),
    "mb_prot": tuple(range(26, 35)),
}

# Parsing tomograms filenames from the CSV file
df = pandas.read_csv(in_csv, delimiter="\t")
tomos = set(df["Tomo3D"].tolist())
segs = dict()
for tomo in tomos:
    tomo_path, tomo_fname = os.path.split(tomo)
    segs[tomo] = tomo_path + "/tomo_lbls_" + tomo_fname.split("_")[2] + ".mrc"
assert len(tomos) == len(segs.keys())

# Create the dataset in nn-UNet format
out_dataset = out_dir + "/Dataset" + dataset_id + "_" + dataset_suffix
if os.path.exists(out_dataset):
    shutil.rmtree(out_dataset)
os.mkdir(out_dataset)
imgs_tr_dir, lbls_ts_dir = out_dataset + "/imagesTr", out_dataset + "/labelsTr"
os.mkdir(imgs_tr_dir)
os.mkdir(lbls_ts_dir)
out_labels = {"background": 0}
for tomo_id, tomo_in in enumerate(tomos):
    print("Processing tomogram:", tomo_in)
    tomo = lio.load_mrc(tomo_in)
    seg = lio.load_mrc(segs[tomo_in])
    seg_post = np.zeros(shape=seg.shape, dtype=np.uint8)
    for i, key in enumerate(fg_labels.keys()):
        print("\tProcessing label:", key)
        for lbl in fg_labels[key]:
            seg_post[seg == lbl] = i + 1
        out_labels[key] = i + 1
        out_labels[key] = i + 1
    nrrd.write(
        imgs_tr_dir + "/tomo_" + str(tomo_id).zfill(3) + "_0000.nrrd", tomo
    )
    nrrd.write(
        lbls_ts_dir + "/tomo_" + str(tomo_id).zfill(3) + ".nrrd", seg_post
    )

# Json configuration file
dict_json = {
    "channel_names": {"0": "rescale_to_0_1"},
    "labels": out_labels,
    "numTraining": len(tomos),
    "file_ending": ".nrrd",
}
with open(out_dataset + "/dataset.json", "w") as outfile:
    outfile.write(json.dumps(dict_json, indent=4))

print("Successfully terminated. (" + time.strftime("%c") + ")")
