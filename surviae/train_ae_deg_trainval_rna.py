import datetime
import filecmp
import os
import sys
import warnings
import numpy as np
import pandas as pd
import matplotlib.pyplot as plt
import logging
from sklearn.metrics import (roc_auc_score, average_precision_score, roc_curve, precision_recall_curve,
                             precision_score, recall_score, accuracy_score)
from sklearn.ensemble import GradientBoostingClassifier, AdaBoostClassifier, RandomForestClassifier, ExtraTreesClassifier
import xgboost as xgb

warnings.filterwarnings('ignore')

import shap
from keras.metrics import AUC, MeanAbsoluteError, MeanSquaredError
from keras.callbacks import EarlyStopping
from keras.saving.save import load_model

sys.path.append(os.path.dirname(os.path.abspath(os.path.join(__file__, '..'))))

from definitions import ROOT_DIR
from surviae.data_deg import FeatureSet, DICT_FEATURESET, DICT_TARGET, DICT_PREPROCESSING_DIRNAME, get_dataset_file
from surviae.model import (build_autoencoder, AE_DIM_DICT, extract_encoder, extract_decoder,
                           extract_hidden, build_classifier, combine_models, extract_encoder_hidden)
from surviae.plot import tsne_plot, umap_plot, find_best_th_min, plot_km, plot_coxph, plot_trend
from surviae.surv import multivariate_surv_analysis, univariate_surv_analysis, multivariate_surv_analysis_lifelines

logging.basicConfig(
    filename=os.path.join(ROOT_DIR, 'logs', 'train_ae_deg_trainval_rna.log'),
    level=logging.INFO
)

if __name__ == '__main__':
    import argparse
    parser = argparse.ArgumentParser(
        prog='SurvIAE',
        description='SurvIAE: Survival prediction with Interpretable Autoencoders from '
                    'Diffuse Large B-Cells Lymphoma gene expression data',
        epilog='Nicola Altini (2023)'
    )
    parser.add_argument('target', type=int, help='0=OS, 1=PFS')
    parser.add_argument('fs', type=int, help='1=RNA, 3=DEG, 5=DEG_XS')
    parser.add_argument('model', type=int, help='-1=BASELINE, 0=XS, 1=S, 2=M, 3=L, 4=XL, 5=XXL')
    parser.add_argument('--epochs-ae', type=int, default=2000, help='Training epochs for autoencoder')
    parser.add_argument('--epochs-clf', type=int, default=200, help='Training epochs for classifier')
    parser.add_argument('--patience-clf', type=int, default=100, help='Training patience for classifier')
    parser.add_argument('--train-ae', type=int, default=1, help='0=load, 1=train')
    parser.add_argument('--train-clf', type=int, default=1, help='0=load, 1=train')
    parser.add_argument('--do-xai', type=int, default=1, help='0=no, 1=yes')
    parser.add_argument('--xai-on-val', type=int, default=1, help='0=sha, 1=chapuy, 2=schmitz')
    parser.add_argument('--follow-up', type=int, default=None, help='months of follow-up to consider')
    parser.add_argument('--preprocessing', type=int, default=0, help='0=no, 1=min-max-scaling, 2=std-scaling')
    args = parser.parse_args()

    target = args.target
    fs = args.fs
    dim_int = args.model
    epochs_ae = args.epochs_ae
    epochs_clf = args.epochs_clf
    patience_clf = args.patience_clf
    train_ae = args.train_ae
    train_clf = args.train_clf
    do_xai = args.do_xai
    xai_on_val = args.xai_on_val
    input_clinical = 0
    follow_up = args.follow_up
    preprocessing = args.preprocessing

    preprocessing_dir = DICT_PREPROCESSING_DIRNAME[preprocessing]
    target_dir = DICT_TARGET[target]
    follow_up_dir = str(follow_up)
    feature_set_dir = DICT_FEATURESET[fs]

    logging.info(f"==============================\n\n"
                 f"Starting run at {datetime.datetime.now()}.\n"
                 f"Config: Preprocessing={preprocessing_dir}; Target={target_dir}-{follow_up}; "
                 f"FeatureSet={feature_set_dir}; Architecture={dim_int}\n"
                 f"Train-AE={train_ae}; Train-CLF={train_clf}")

    figure_path = os.path.join(ROOT_DIR, 'figures', preprocessing_dir, f"{target_dir}{follow_up_dir}-{feature_set_dir}")
    logs_path = os.path.join(ROOT_DIR, 'logs', preprocessing_dir, target_dir, follow_up_dir, feature_set_dir)
    os.makedirs(figure_path, exist_ok=True)
    os.makedirs(logs_path, exist_ok=True)

    do_ae = dim_int >= 0
    dimension = AE_DIM_DICT[dim_int]

    filename_autoencoder = os.path.join(logs_path, f'{dimension.__name__}_target_{target}_fs_{fs}.h5')
    filename_classifier = os.path.join(logs_path, f'clf_{dimension.__name__}_target_{target}_fs_{fs}.h5')

    data_sha_train, data_sha_val, data_chapuy, data_schmitz = get_dataset_file(preprocessing, follow_up, target, fs)

    X_deg_sha_train, y_sha_train, time_sha_train = data_sha_train
    X_deg_sha_val, y_sha_val, time_sha_val = data_sha_val
    X_deg_chapuy, y_chapuy, time_chapuy = data_chapuy
    X_deg_schmitz, y_schmitz, time_schmitz = data_schmitz

    X_deg_sha = pd.concat([X_deg_sha_train, X_deg_sha_val], axis=0).astype('float64')
    y_sha = pd.concat([y_sha_train, y_sha_val], axis=0).astype('bool')
    time_sha = pd.concat([time_sha_train, time_sha_val], axis=0).astype('float64')

    print(f"X (Sha) [{X_deg_sha.dtypes}]: {X_deg_sha.shape}")
    print(f"y (Sha) [{y_sha.dtypes}]: {y_sha.shape}")
    print(f"t (Sha) [{time_sha.dtypes}]: {time_sha.shape}")

    print(f"X (Chapuy) [{X_deg_chapuy.dtypes}]: {X_deg_chapuy.shape}")
    print(f"y (Chapuy) [{y_chapuy.dtypes}]: {y_chapuy.shape}")
    print(f"t (Chapuy) [{time_chapuy.dtypes}]: {time_chapuy.shape}")

    if follow_up:
        y_sha_train_fu = np.logical_and(y_sha_train, time_sha_train <= follow_up)
        y_sha_val_fu = np.logical_and(y_sha_val, time_sha_val <= follow_up)
        y_sha_fu = np.logical_and(y_sha, time_sha <= follow_up)
        y_chapuy_fu = np.logical_and(y_chapuy, time_chapuy <= follow_up)
        y_schmitz_fu = np.logical_and(y_schmitz, time_schmitz <= follow_up)

        print(f"y (Sha) Train [{y_sha_train_fu.dtypes}]: {y_sha_train_fu.shape}")
        print(f"y (Sha) Val   [{y_sha_val_fu.dtypes}]: {y_sha_val_fu.shape}")
        print(f"y (Sha)       [{y_sha_fu.dtypes}]: {y_sha_fu.shape}")
        print(f"y (Chapuy)    [{y_chapuy_fu.dtypes}]: {y_chapuy_fu.shape}")
        print(f"y (Schmitz)   [{y_schmitz_fu.dtypes}]: {y_schmitz_fu.shape}")
    else:
        y_sha_train_fu = y_sha_train
        y_sha_val_fu = y_sha_val
        y_sha_fu = y_sha
        y_chapuy_fu = y_chapuy
        y_schmitz_fu = y_schmitz

    print(f"Dataset (Sha)     ratio: {y_sha.sum() / len(y_sha):.3f}, FU-{follow_up} ratio: {y_sha_fu.sum() / len(y_sha_fu):.3f}")
    print(f"Dataset (Chapuy)  ratio: {y_chapuy.sum() / len(y_chapuy):.3f}, FU-{follow_up} ratio: {y_chapuy_fu.sum() / len(y_chapuy_fu):.3f}")
    print(f"Dataset (Schmitz) ratio: {y_schmitz.sum() / len(y_schmitz):.3f}, FU-{follow_up} ratio: {y_schmitz_fu.sum() / len(y_schmitz_fu):.3f}")

    logging.info(f"Val Data Info: {X_deg_sha_val.min().min()}, {X_deg_sha_val.min().mean()}, {X_deg_sha_val.mean().mean()}, "
                 f"{X_deg_sha_val.max().mean()}, {X_deg_sha_val.max().max()}")

    X_deg = np.concatenate([X_deg_sha, X_deg_chapuy, X_deg_schmitz])
    label = np.concatenate([y_sha_fu, y_chapuy_fu, y_schmitz_fu])
    label_str = np.where(label == 1, f"{DICT_TARGET[target]}{follow_up}-1", f"{DICT_TARGET[target]}{follow_up}-0")
    subset = ["Train"] * len(y_sha) + ["Val"] * len(y_chapuy) + ["Test"] * len(y_schmitz)
    input_rna = X_deg_sha_train.shape[1]

    if not do_ae:
        x_train = X_deg_sha
        X_valid = X_deg_chapuy
        X_test = X_deg_schmitz

        gbc = GradientBoostingClassifier(n_estimators=100, learning_rate=1.0, max_depth=1, random_state=0)
        abc = AdaBoostClassifier(n_estimators=100, random_state=0)
        rfc = RandomForestClassifier(n_estimators=10, max_depth=None, min_samples_split=2, random_state=0)
        xtc = ExtraTreesClassifier(n_estimators=10, max_depth=None, min_samples_split=2, random_state=0)
        xgc = xgb.XGBClassifier()

        gbc.fit(x_train, y_sha_fu)
        abc.fit(x_train, y_sha_fu)
        rfc.fit(x_train, y_sha_fu)
        xtc.fit(x_train, y_sha_fu)
        xgc.fit(x_train, y_sha_fu)

        for clf, name in zip([gbc, abc, rfc, xtc, xgc],
                             ["GradientBoostingClassifier",
                              "AdaBoostClassifier",
                              "RandomForestClassifier",
                              "ExtraTreesClassifier",
                              "XGBClassifier"]):
            p_train = clf.predict_proba(x_train)[:, 1]
            p_valid = clf.predict_proba(X_valid)[:, 1]
            p_test = clf.predict_proba(X_test)[:, 1]

            auc_roc_train = roc_auc_score(y_sha_fu, p_train)
            auc_roc_valid = roc_auc_score(y_chapuy_fu, p_valid)
            auc_roc_test = roc_auc_score(y_schmitz_fu, p_test)

            auc_pr_train = average_precision_score(y_sha_fu, p_train)
            auc_pr_valid = average_precision_score(y_chapuy_fu, p_valid)
            auc_pr_test = average_precision_score(y_schmitz_fu, p_test)

            best_threshold_train = find_best_th_min(y_sha_fu, p_train)
            best_threshold_valid = find_best_th_min(y_chapuy_fu, p_valid)
            best_threshold_test = find_best_th_min(y_schmitz_fu, p_test)

            y_hat_train = (p_train > best_threshold_train).astype(int)
            y_hat_valid = (p_valid > best_threshold_valid).astype(int)
            y_hat_test = (p_test > best_threshold_test).astype(int)

            acc_train = accuracy_score(y_sha_fu, y_hat_train)
            acc_valid = accuracy_score(y_chapuy_fu, y_hat_valid)
            acc_test = accuracy_score(y_schmitz_fu, y_hat_test)

            prec_train = precision_score(y_sha_fu, y_hat_train)
            prec_valid = precision_score(y_chapuy_fu, y_hat_valid)
            prec_test = precision_score(y_schmitz_fu, y_hat_test)

            reca_train = recall_score(y_sha_fu, y_hat_train)
            reca_valid = recall_score(y_chapuy_fu, y_hat_valid)
            reca_test = recall_score(y_schmitz_fu, y_hat_test)

            logging.info(f"[AUROC] Train: {auc_roc_train:.3f}, Val: {auc_roc_valid:.3f}, Test: {auc_roc_test:.3f}")
            logging.info(f"[AUPRC] Train: {auc_pr_train:.3f}, Val: {auc_pr_valid:.3f}, Test: {auc_pr_test:.3f}")

            dict_metrics = {
                'Metric': ['AUC-ROC', 'AUC-PR', 'Accuracy', 'Precision', 'Recall', 'Threshold'],
                'Sha': [auc_roc_train, auc_pr_train, acc_train, prec_train, reca_train, best_threshold_train],
                'Chapuy': [auc_roc_valid, auc_pr_valid, acc_valid, prec_valid, reca_valid, best_threshold_valid],
                'Schmitz': [auc_roc_test, auc_pr_test, acc_test, prec_test, reca_test, best_threshold_test],
            }

            df_metrics = pd.DataFrame(data=dict_metrics)
            df_metrics.to_csv(os.path.join(figure_path, f'results_METRICS_{dimension.__name__}_{name}_fs_{fs}_target_{target}.csv'))

    else:
        tsne_plot({'Features': X_deg, 'Label': label_str, 'Subset': subset}, n_components=2,
                  fig_path=os.path.join(figure_path,
                                        f'fig_TSNE_original_{dimension.__name__}_fs_{fs}_target_{target}.png'))
        umap_plot({'Features': X_deg, 'Label': label_str, 'Subset': subset}, n_components=2,
                  fig_path=os.path.join(figure_path,
                                        f'fig_UMAP_original_{dimension.__name__}_fs_{fs}_target_{target}.png'))

        if train_ae:
            print("Training the autoencoder!")
            autoencoder, encoder, decoder = build_autoencoder(dimension=dimension, input_rna=input_rna)
            autoencoder.compile(optimizer='adam', loss='mse', metrics=[MeanSquaredError(), MeanAbsoluteError()])
            # Train the autoencoder
            autoencoder.fit(X_deg_sha, X_deg_sha,
                            epochs=epochs_ae,
                            batch_size=64,
                            shuffle=True,
                            validation_data=(X_deg_chapuy, X_deg_chapuy))
            logging.info(f"Saving Autoencoder to {filename_autoencoder}")
            autoencoder.save_weights(filename_autoencoder)
            autoencoder.summary()

        else:
            logging.info(f"Loading Autoencoder from {filename_autoencoder}")
            autoencoder, _, _ = build_autoencoder(dimension=dimension, input_rna=input_rna)
            autoencoder.compile(optimizer='adam', loss='mse', metrics=[MeanSquaredError(), MeanAbsoluteError()])
            autoencoder.load_weights(filename_autoencoder)

            autoencoder.summary()
            encoder = extract_encoder(autoencoder)
            decoder = extract_decoder(autoencoder)

        # Evaluate the autoencoder
        mse_train = autoencoder.evaluate(X_deg_sha, X_deg_sha, verbose=0)
        mse_valid = autoencoder.evaluate(X_deg_chapuy, X_deg_chapuy, verbose=0)
        mse_test = autoencoder.evaluate(X_deg_schmitz, X_deg_schmitz, verbose=0)

        logging.info(f'Train  (Sha)     {len(mse_train)} MSE: {mse_train[0]:.3f}, MAE: {mse_train[2]:.3f}')
        logging.info(f'Valid  (Chapuy)  {len(mse_valid)} MSE: {mse_valid[0]:.3f}, MAE: {mse_valid[2]:.3f}')
        logging.info(f'Test   (Schmitz) {len(mse_test)} MSE: {mse_test[0]:.3f}, MAE: {mse_test[2]:.3f}')

        X_deg_train_pred = autoencoder.predict(X_deg_sha)
        X_deg_train_enc = encoder.predict(X_deg_sha)
        X_deg_train_pred_v2 = decoder.predict(X_deg_train_enc)

        X_deg_valid_pred = autoencoder.predict(X_deg_chapuy)
        X_deg_valid_enc = encoder.predict(X_deg_chapuy)
        X_deg_valid_pred_v2 = decoder.predict(X_deg_valid_enc)

        X_deg_test_pred = autoencoder.predict(X_deg_schmitz)
        X_deg_test_enc = encoder.predict(X_deg_schmitz)
        X_deg_test_pred_v2 = decoder.predict(X_deg_test_enc)

        logging.debug(f"(Train) Encoded shape: {X_deg_train_enc.shape}")
        logging.debug(f"(Valid) Encoded shape: {X_deg_valid_enc.shape}")
        logging.debug(f"(Test)  Encoded shape: {X_deg_test_enc.shape}")

        base_clf = build_classifier(dimension=dimension, input_encoded=X_deg_train_enc.shape[1], use_hidden_1=False)
        clf = combine_models(encoder, base_clf)

        X_deg = np.concatenate([X_deg_train_enc, X_deg_valid_enc, X_deg_test_enc])

        tsne_plot({'Features': X_deg, 'Label': label_str, 'Subset': subset}, n_components=2, fig_path=os.path.join(figure_path, f'fig_TSNE_complete_{dimension.__name__}_fs_{fs}_target_{target}.png'))
        umap_plot({'Features': X_deg, 'Label': label_str, 'Subset': subset}, n_components=2, fig_path=os.path.join(figure_path, f'fig_UMAP_complete_{dimension.__name__}_fs_{fs}_target_{target}.png'))

        assert np.allclose(X_deg_train_pred, X_deg_train_pred_v2)
        assert np.allclose(X_deg_valid_pred, X_deg_valid_pred_v2)
        assert np.allclose(X_deg_test_pred, X_deg_test_pred_v2)

        # Train the classifier
        clf.compile(optimizer='adam', loss='binary_crossentropy', metrics=['accuracy', AUC(curve="ROC"), AUC(curve="PR")])

        x_train = X_deg_sha
        X_valid = X_deg_chapuy
        X_test = X_deg_schmitz

        if train_clf:
            callback = EarlyStopping(monitor='val_loss', patience=patience_clf, restore_best_weights=True)
            clf.fit(x_train, y_sha_fu,
                    epochs=epochs_clf,
                    batch_size=64,
                    shuffle=True,
                    callbacks=[callback],
                    validation_data=(X_valid, y_chapuy_fu))

            print(f"Saving Classifier to {filename_classifier}...")
            clf.save(filename_classifier)
            clf.summary()
            hidden = extract_encoder_hidden(clf)

        else:
            print("Loading the Classifier!")
            clf = load_model(filename_classifier)
            clf.summary()
            hidden = extract_encoder_hidden(clf)

        X_train_hid = hidden.predict(x_train)
        X_valid_hid = hidden.predict(X_valid)
        X_test_hid = hidden.predict(X_test)

        logging.debug(f"(Train) Hidden shape: {X_train_hid.shape}")
        logging.debug(f"(Valid) Hidden shape: {X_valid_hid.shape}")
        logging.debug(f"(Test)  Hidden shape: {X_test_hid.shape}")

        X_hid = np.concatenate([X_train_hid, X_valid_hid, X_test_hid])

        tsne_plot({'Features': X_hid, 'Label': label_str, 'Subset': subset}, n_components=2, fig_path=os.path.join(figure_path, f'fig_hidden_TSNE_complete_{dimension.__name__}_fs_{fs}_target_{target}.png'))
        umap_plot({'Features': X_hid, 'Label': label_str, 'Subset': subset}, n_components=2, fig_path=os.path.join(figure_path, f'fig_hidden_UMAP_complete_{dimension.__name__}_fs_{fs}_target_{target}.png'))

        # Evaluate the classifier
        ce_train = clf.evaluate(x_train, y_sha_fu, verbose=0)
        ce_valid = clf.evaluate(X_valid, y_chapuy_fu, verbose=0)
        ce_test = clf.evaluate(X_test, y_schmitz_fu, verbose=0)

        logging.info(f'Dimension: {dimension}')
        logging.info(f'Train (Sha)     CE: {ce_train[0]:.3f}, Acc: {ce_train[1]:.3f}')
        logging.info(f'Valid (Chapuy)  CE: {ce_valid[0]:.3f}, Acc: {ce_valid[1]:.3f}')
        logging.info(f'Test  (Schmitz) CE: {ce_test[0]:.3f}, Acc: {ce_test[1]:.3f}')

        p_train = clf.predict(x_train)[:, 0]
        p_valid = clf.predict(X_valid)[:, 0]
        p_test = clf.predict(X_test)[:, 0]

        y_train_thrs = []
        y_val_thrs = []
        y_test_thrs = []
        thrs_time = []

        aucs_train = []
        aucs_valid = []
        aucs_test = []

        for thr_time in range(12, 210, 12):
            thrs_time.append(thr_time)

            y_train_thr = np.logical_and(y_sha, time_sha <= thr_time)
            y_valid_thr = np.logical_and(y_chapuy, time_chapuy <= thr_time)
            y_test_thr = np.logical_and(y_schmitz, time_schmitz <= thr_time)

            y_train_thrs.append(y_train_thr)
            y_val_thrs.append(y_valid_thr)
            y_test_thrs.append(y_test_thr)

            auc_roc_train_thr = roc_auc_score(y_train_thr, p_train)
            auc_roc_valid_thr = roc_auc_score(y_valid_thr, p_valid)
            auc_roc_test_thr = roc_auc_score(y_test_thr, p_test)

            aucs_train.append(auc_roc_train_thr)
            aucs_valid.append(auc_roc_valid_thr)
            aucs_test.append(auc_roc_test_thr)

        plot_trend(thrs_time, [aucs_train, aucs_valid, aucs_test], ["Sha", "Chapuy", "Schmitz"],
                   filename=os.path.join(figure_path, f'fig_FU_AUC_{dimension.__name__}_fs_{fs}_target_{target}.png'))

        fpr_train, tpr_train, thresholds_train = roc_curve(y_sha_fu, p_train)
        fpr_valid, tpr_valid, thresholds_valid = roc_curve(y_chapuy_fu, p_valid)
        fpr_test, tpr_test, thresholds_test = roc_curve(y_schmitz_fu, p_test)

        precision_train, recall_train, _ = precision_recall_curve(y_sha_fu, p_train)
        precision_valid, recall_valid, _ = precision_recall_curve(y_chapuy_fu, p_valid)
        precision_test, recall_test, _ = precision_recall_curve(y_schmitz_fu, p_test)

        auc_roc_train = roc_auc_score(y_sha_fu, p_train)
        auc_roc_valid = roc_auc_score(y_chapuy_fu, p_valid)
        auc_roc_test = roc_auc_score(y_schmitz_fu, p_test)

        auc_pr_train = average_precision_score(y_sha_fu, p_train)
        auc_pr_valid = average_precision_score(y_chapuy_fu, p_valid)
        auc_pr_test = average_precision_score(y_schmitz_fu, p_test)

        best_threshold_train = find_best_th_min(y_sha_fu, p_train)
        best_threshold_valid = find_best_th_min(y_chapuy_fu, p_valid)
        best_threshold_test = find_best_th_min(y_schmitz_fu, p_test)

        logging.info(f"Best Threshold, (Theory) 0.5, (Sha) {best_threshold_train:.3f}, "
                     f"(Chapuy) {best_threshold_valid:.3f}, (Schmitz) {best_threshold_test:.3f}")

        logging.info(f'AUC ROC: {auc_roc_train:.3f}, {auc_roc_valid:.3f}, {auc_roc_test:.3f}')
        logging.info(f'AUC PR : {auc_pr_train:.3f}, {auc_pr_valid:.3f}, {auc_pr_test:.3f}')

        y_hat_train = (p_train > best_threshold_train).astype(int)
        y_hat_valid = (p_valid > best_threshold_valid).astype(int)
        y_hat_test = (p_test > best_threshold_test).astype(int)

        plot_km(y_sha, time_sha, y_hat_train, os.path.join(figure_path, f'fig_KM_sha_{dimension.__name__}_fs_{fs}_target_{target}.png'), show_confidence=False)
        plot_km(y_chapuy, time_chapuy, y_hat_valid, os.path.join(figure_path, f'fig_KM_chapuy_{dimension.__name__}_fs_{fs}_target_{target}.png'), show_confidence=False)
        plot_km(y_schmitz, time_schmitz, y_hat_test, os.path.join(figure_path, f'fig_KM_schmitz_{dimension.__name__}_fs_{fs}_target_{target}.png'), show_confidence=False)

        acc_train = accuracy_score(y_sha_fu, y_hat_train)
        acc_valid = accuracy_score(y_chapuy_fu, y_hat_valid)
        acc_test = accuracy_score(y_schmitz_fu, y_hat_test)

        prec_train = precision_score(y_sha_fu, y_hat_train)
        prec_valid = precision_score(y_chapuy_fu, y_hat_valid)
        prec_test = precision_score(y_schmitz_fu, y_hat_test)

        reca_train = recall_score(y_sha_fu, y_hat_train)
        reca_valid = recall_score(y_chapuy_fu, y_hat_valid)
        reca_test = recall_score(y_schmitz_fu, y_hat_test)

        logging.info(f'Accuracy : {acc_train:.3f}, {acc_valid:.3f}, {acc_test:.3f}')
        logging.info(f'Precision: {prec_train:.3f}, {prec_valid:.3f}, {prec_test:.3f}')
        logging.info(f'Recall   : {reca_train:.3f}, {reca_valid:.3f}, {reca_test:.3f}')

        dict_metrics = {
            'Metric': ['AUC-ROC', 'AUC-PR', 'Accuracy', 'Precision', 'Recall', 'Threshold'],
            'Sha': [auc_roc_train, auc_pr_train, acc_train, prec_train, reca_train, best_threshold_train],
            'Chapuy': [auc_roc_valid, auc_pr_valid, acc_valid, prec_valid, reca_valid, best_threshold_valid],
            'Schmitz': [auc_roc_test, auc_pr_test, acc_test, prec_test, reca_test, best_threshold_test],
        }

        dict_metrics_ae = {
            'Metric': ['MSE', 'MAE'],
            'Sha': [mse_train[0], mse_train[2]],
            'Chapuy': [mse_valid[0], mse_valid[2]],
            'Schmitz': [mse_test[0], mse_test[2]],
        }

        dict_pred_train = {
            'Patient': x_train.index,
            'Prob': p_train,
            'Risk': y_hat_train
        }

        dict_pred_valid = {
            'Patient': X_valid.index,
            'Prob': p_valid,
            'Risk': y_hat_valid
        }

        dict_pred_test = {
            'Patient': X_test.index,
            'Prob': p_test,
            'Risk': y_hat_test
        }

        df_pred_train = pd.DataFrame(data=dict_pred_train)
        df_pred_train.to_csv(os.path.join(figure_path, f'results_RISK_sha_{dimension.__name__}_fs_{fs}_target_{target}.csv'))

        df_pred_valid = pd.DataFrame(data=dict_pred_valid)
        df_pred_valid.to_csv(os.path.join(figure_path, f'results_RISK_chapuy_{dimension.__name__}_fs_{fs}_target_{target}.csv'))

        df_pred_test = pd.DataFrame(data=dict_pred_test)
        df_pred_test.to_csv(os.path.join(figure_path, f'results_RISK_schmitz_{dimension.__name__}_fs_{fs}_target_{target}.csv'))

        df_metrics = pd.DataFrame(data=dict_metrics)
        df_metrics.to_csv(os.path.join(figure_path, f'results_METRICS_{dimension.__name__}_fs_{fs}_target_{target}.csv'))

        df_metrics_ae = pd.DataFrame(data=dict_metrics_ae)
        df_metrics_ae.to_csv(os.path.join(figure_path, f'results_METRICS_AE_{dimension.__name__}_fs_{fs}_target_{target}.csv'))

        xy_label_fontsize = 24
        legend_fontsize = 20
        linewidth = 3.0

        # Plot the ROC curve
        plt.figure(figsize=(10, 10))
        plt.style.use('classic')
        plt.plot(fpr_train, tpr_train, label=f'Train AUROC = {round(auc_roc_train, 3):.3f}', color='g', linewidth=linewidth)
        plt.plot(fpr_valid, tpr_valid, label=f'Val AUROC = {round(auc_roc_valid, 3):.3f}', color='b', linewidth=linewidth)
        plt.plot(fpr_test, tpr_test, label=f'Test AUROC = {round(auc_roc_test, 3):.3f}', color='c', linewidth=linewidth)
        plt.plot([0, 1], [0, 1], '--', color='y', label='Random Guess', linewidth=linewidth)
        plt.xlim([-0.02, 1.02])
        plt.ylim([-0.02, 1.02])
        plt.xticks(fontsize=legend_fontsize)
        plt.yticks(fontsize=legend_fontsize)
        plt.xlabel('False Positive Rate', fontsize=xy_label_fontsize)
        plt.ylabel('True Positive Rate', fontsize=xy_label_fontsize)
        plt.legend(loc='lower right', fontsize=legend_fontsize)
        plt.grid(True)
        plt.savefig(os.path.join(figure_path, f'fig_ROC_{dimension.__name__}_fs_{fs}_target_{target}.png'))
        plt.close()

        plt.figure(figsize=(10, 10))
        plt.style.use('classic')
        plt.plot(recall_train, precision_train, label=f'Train AUPRC = {round(auc_pr_train, 3):.3f}', color='g', linewidth=linewidth)
        plt.plot(recall_valid, precision_valid, label=f'Val AUPRC = {round(auc_pr_valid, 3):.3f}', color='b', linewidth=linewidth)
        plt.plot(recall_test, precision_test, label=f'Test AUPRC = {round(auc_pr_test, 3):.3f}', color='c', linewidth=linewidth)
        plt.xlim([-0.02, 1.02])
        plt.ylim([-0.02, 1.02])
        plt.xticks(fontsize=legend_fontsize)
        plt.yticks(fontsize=legend_fontsize)
        plt.xlabel('Recall', fontsize=xy_label_fontsize)
        plt.ylabel('Precision', fontsize=xy_label_fontsize)
        plt.legend(loc='upper right', fontsize=legend_fontsize)
        plt.grid(True)
        plt.savefig(os.path.join(figure_path, f'fig_PR_{dimension.__name__}_fs_{fs}_target_{target}.png'))
        plt.close()

        if do_xai:
            print("Doing XAI...")

            X_deg_train_np = X_deg_sha.values
            X_deg_val_np = X_deg_chapuy.values
            X_deg_test_np = X_deg_schmitz.values

            X_train_np = X_deg_train_np
            X_val_np = X_deg_val_np
            X_test_np = X_deg_test_np

            clf_layer = 'clf_layer'
            explainer = shap.DeepExplainer(clf, X_train_np)
            shap.explainers._deep.deep_tf.op_handlers["AddV2"] = shap.explainers._deep.deep_tf.passthrough

            if xai_on_val == 1:
                shap_values = explainer.shap_values(X_val_np, check_additivity=False)
                X_deg_shap = X_deg_val_np
            elif xai_on_val == 2:
                shap_values = explainer.shap_values(X_test_np, check_additivity=False)
                X_deg_shap = X_deg_test_np
            else:
                shap_values = explainer.shap_values(X_train_np, check_additivity=False)
                X_deg_shap = X_deg_train_np

            shap_values_deg = shap_values[0]

            feature_names_deg = X_deg_sha_val.columns.to_numpy()

            dict_deg_shap = {
                'Gene': feature_names_deg,
                'SHAP Value': np.mean(np.absolute(shap_values_deg), axis=0)
            }

            df_deg_shap = pd.DataFrame(data=dict_deg_shap)
            df_deg_shap.to_csv(os.path.join(figure_path, f'results_SHAP_{dimension.__name__}_fs_{fs}_target_{target}_val_{xai_on_val}.csv'))

            # Summarize the SHAP values
            plt.figure()
            shap.summary_plot(shap_values_deg, X_deg_shap, feature_names=feature_names_deg, show=False, plot_type="bar")
            plt.savefig(os.path.join(figure_path, f'fig_SHAP_{dimension.__name__}_fs_{fs}_target_{target}_bar_val_{xai_on_val}.png'))
            plt.tight_layout()
            plt.close()

            plt.figure()
            shap.summary_plot(shap_values_deg, X_deg_shap, feature_names=feature_names_deg, show=False)
            plt.savefig(os.path.join(figure_path, f'fig_SHAP_{dimension.__name__}_fs_{fs}_target_{target}_val_{xai_on_val}.png'))
            plt.tight_layout()
            plt.close()
