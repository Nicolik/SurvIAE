import os
from enum import IntEnum
import pandas as pd
from definitions import ROOT_DIR


class FeatureSet(IntEnum):
    CLINIC = 0
    RNA = 1
    CLINIC_RNA = 2
    DEG = 3
    CLINIC_DEG = 4
    DEG_XS = 5
    CLINIC_DEG2 = 6


class Target(IntEnum):
    OS = 0
    PFS = 1


class Preprocessing(IntEnum):
    NONE = 0
    MIN_MAX_SCALER = 1
    STANDARD_SCALER = 2


DICT_PREPROCESSING_DIRNAME = {
    Preprocessing.NONE: "genes-clf-nonorm",
    Preprocessing.MIN_MAX_SCALER: "genes-clf-minmax",
    Preprocessing.STANDARD_SCALER: "genes-clf-standard"
}


DICT_FEATURESET = {
    FeatureSet.CLINIC: "Clinical",
    FeatureSet.RNA: "RNA",
    FeatureSet.CLINIC_RNA: "Clinical+RNA",
    FeatureSet.DEG: "DEG",
    FeatureSet.CLINIC_DEG: "Clinical+DEG"
}


DICT_TARGET = {
    Target.OS: "OS",
    Target.PFS: "PFS"
}


def get_os(clinic_df, get_event_time=False):
    clinic_df['OS'] = clinic_df['OS'].astype(int)
    clinic_df['OS_time'] = clinic_df['OS_time'].astype(float)
    if get_event_time:
        return clinic_df[['OS', 'OS_time']]
    else:
        return clinic_df['OS']


def get_pfs(clinic_df, get_event_time=False):
    clinic_df['PFS'] = clinic_df['PFS'].astype(int)
    clinic_df['PFS_time'] = clinic_df['PFS_time'].astype(float)
    if get_event_time:
        return clinic_df[['PFS', 'PFS_time']]
    else:
        return clinic_df['PFS']


def get_dataset_file(preprocessing, fu,  target, feature_set):
    out_dir = os.path.join(ROOT_DIR, 'dataset')
    preproc_dir = DICT_PREPROCESSING_DIRNAME[preprocessing]
    target_dir = DICT_TARGET[target]
    feature_set_dir = DICT_FEATURESET[feature_set]
    data_dir = os.path.join(out_dir, preproc_dir, f"{target_dir}{fu}-{feature_set_dir}")

    X_train_read = pd.read_csv(os.path.join(data_dir, 'X_train.csv'), index_col=0)
    y_train_read = pd.read_csv(os.path.join(data_dir, 'y_train.csv'), index_col=0).iloc[:, 0]
    t_train_read = pd.read_csv(os.path.join(data_dir, 't_train.csv'), index_col=0).iloc[:, 0]

    X_val_read = pd.read_csv(os.path.join(data_dir, 'X_val.csv'), index_col=0)
    y_val_read = pd.read_csv(os.path.join(data_dir, 'y_val.csv'), index_col=0).iloc[:, 0]
    t_val_read = pd.read_csv(os.path.join(data_dir, 't_val.csv'), index_col=0).iloc[:, 0]

    X_test_read = pd.read_csv(os.path.join(data_dir, 'X_test.csv'), index_col=0)
    y_test_read = pd.read_csv(os.path.join(data_dir, 'y_test.csv'), index_col=0).iloc[:, 0]
    t_test_read = pd.read_csv(os.path.join(data_dir, 't_test.csv'), index_col=0).iloc[:, 0]

    X_schmitz_read = pd.read_csv(os.path.join(data_dir, 'X_schmitz.csv'), index_col=0)
    y_schmitz_read = pd.read_csv(os.path.join(data_dir, 'y_schmitz.csv'), index_col=0).iloc[:, 0]
    t_schmitz_read = pd.read_csv(os.path.join(data_dir, 't_schmitz.csv'), index_col=0).iloc[:, 0]

    return ((X_train_read, y_train_read, t_train_read), (X_val_read, y_val_read, t_val_read),
            (X_test_read, y_test_read, t_test_read), (X_schmitz_read, y_schmitz_read, t_schmitz_read))
