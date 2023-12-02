import os
from collections import defaultdict
import pandas as pd

from surviae.data_deg import (DICT_PREPROCESSING_DIRNAME, DICT_TARGET, DICT_FEATURESET,
                              Preprocessing, Target, FeatureSet)
from surviae.model import AutoencoderSize, AE_DIM_DICT
from definitions import ROOT_DIR

figure_base_dir = os.path.join(ROOT_DIR, 'figures')

data_dict = {}
data_dict_val = defaultdict(list)
data_dict_test = defaultdict(list)

model = AutoencoderSize.BASELINE
dimension = AE_DIM_DICT[model]

for name in ["GradientBoostingClassifier", "AdaBoostClassifier", "RandomForestClassifier",
             "ExtraTreesClassifier", "XGBClassifier"]:
    preprocessing_model = []

    for target in [Target.OS, Target.PFS]:
        for fu in [12, 36, 60]:
            for fs in [FeatureSet.DEG]:
                for preprocessing in [Preprocessing.NONE, Preprocessing.MIN_MAX_SCALER,
                                      Preprocessing.STANDARD_SCALER]:
                    preprocessing_dir = DICT_PREPROCESSING_DIRNAME[preprocessing]
                    target_dir = DICT_TARGET[target]
                    feature_set_dir = DICT_FEATURESET[fs]
                    figure_path = os.path.join(figure_base_dir, preprocessing_dir, f"{target_dir}{fu}-{feature_set_dir}")

                    model_metrics_csv = os.path.join(figure_path, f'results_METRICS_{dimension.__name__}_{name}_fs_{fs}_target_{target}.csv')

                    if os.path.exists(model_metrics_csv):
                        model_metrics_df = pd.read_csv(model_metrics_csv, index_col=0)
                        auc_roc = model_metrics_df.iloc[0, :]
                        auc_pr = model_metrics_df.iloc[1, :]
                        auc_roc_rounded = [round(a, 3) for a in auc_roc[1:]]
                        auc_pr_rounded = [round(a, 3) for a in auc_pr[1:]]
                        preprocessing_model.extend(list(auc_roc_rounded))
                        preprocessing_model.extend(list(auc_pr_rounded))

                        auc_roc_r_val = auc_roc_rounded[1]
                        auc_roc_r_test = auc_roc_rounded[2]
                        auc_pr_r_val = auc_pr_rounded[1]
                        auc_pr_r_test = auc_pr_rounded[2]

                        data_dict_val[f'{target_dir}{fu}-ROC-AUC'].append(auc_roc_r_val)
                        data_dict_val[f'{target_dir}{fu}-PR-AUC'].append(auc_pr_r_val)

                        data_dict_test[f'{target_dir}{fu}-ROC-AUC'].append(auc_roc_r_test)
                        data_dict_test[f'{target_dir}{fu}-PR-AUC'].append(auc_pr_r_test)

                    else:
                        preprocessing_model.extend([-1] * 6)
        ae_name = f'{preprocessing_dir}-{dimension.__name__}'
        data_dict[ae_name] = preprocessing_model


data_val_df = pd.DataFrame(data=data_dict_val)
data_val_df.to_csv(os.path.join(figure_base_dir, 'metrics_baseline_val.csv'), sep=';')

data_test_df = pd.DataFrame(data=data_dict_test)
data_test_df.to_csv(os.path.join(figure_base_dir, 'metrics_baseline_test.csv'), sep=';')
