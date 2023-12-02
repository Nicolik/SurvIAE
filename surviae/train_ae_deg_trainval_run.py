import os
import sys
import json
import subprocess

sys.path.append(os.path.dirname(os.path.abspath(os.path.join(__file__, '..'))))

from surviae.data_deg import Target, FeatureSet, Preprocessing
from surviae.model import AutoencoderSize
from definitions import ROOT_DIR

script_path = os.path.join(ROOT_DIR, "surviae", "train_ae_deg_trainval_rna.py")
target = Target.OS
feature_set = FeatureSet.DEG
model = AutoencoderSize.XS
preprocessing = Preprocessing.NONE

train_ae = 0
train_clf = 0
input_clinical = 0
do_xai = 1
xai_on_val = 1

epochs_ae = 2000
epochs_clf = 200

for xai_on_val in [0, 1, 2]:
    for preprocessing in [Preprocessing.NONE, Preprocessing.MIN_MAX_SCALER, Preprocessing.STANDARD_SCALER]:
        for target in [Target.OS, Target.PFS]:
            for fu in [12, 36, 60]:
                for feature_set in [FeatureSet.DEG]:
                    for model in [AutoencoderSize.BASELINE,
                                  AutoencoderSize.XS, AutoencoderSize.S, AutoencoderSize.M, AutoencoderSize.L]:
                        subprocess.run(["python", script_path, str(int(target)), str(int(feature_set)), str(int(model)),
                                        "--train-ae", str(train_ae),
                                        "--epochs-ae", str(epochs_ae),
                                        "--train-clf", str(train_clf),
                                        "--epochs-clf", str(epochs_clf),
                                        "--do-xai", str(do_xai),
                                        "--xai-on-val", str(xai_on_val),
                                        "--preprocessing", str(int(preprocessing)),
                                        "--follow-up", str(fu)
                                        ])
