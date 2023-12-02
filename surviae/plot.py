import numpy as np
import matplotlib.pyplot as plt
import seaborn as sns
from sklearn.manifold import TSNE
from sklearn.metrics import confusion_matrix, roc_curve, silhouette_score, silhouette_samples
from umap import UMAP
from sksurv.nonparametric import kaplan_meier_estimator

from surviae.surv import multivariate_surv_analysis_lifelines


def tsne_plot(data, n_components=2, fig_path=None):
    tsne = TSNE(n_components=n_components)
    X = tsne.fit_transform(data["Features"])
    scatter_plot(X, data, fig_path=fig_path, title="TSNE")


def umap_plot(data, n_components=2, fig_path=None):
    umap = UMAP(n_components=n_components)
    X = umap.fit_transform(data["Features"])
    scatter_plot(X, data, fig_path=fig_path, title="UMAP", silhouette=True)


def scatter_plot(X, data, fig_path=None, title=None, silhouette=False):
    plt.figure(figsize=(10, 10))
    style = data["Subset"] if "Subset" in data else None
    hue = data["Label"] if "Label" in data else None
    sns.set(font_scale=2)
    sns.set_theme(style='white')
    sns.scatterplot(x=X[:, 0], y=X[:, 1], hue=hue, style=style, legend="full", s=100)
    xl, xr = plt.xlim()
    yl, yr = plt.ylim()

    if silhouette:
        silhouette = silhouette_score(X, data["Subset"])
        sv = 1
        s = sv * 50
        if silhouette > 0.2:
            c = (sv, 0, 0)
            marker = "v"
        elif silhouette > 0:
            sv = .5
            c = (sv, sv, sv)
            marker = "s"
        else:
            c = (0, sv, 0)
            marker = "^"
        # plt.autoscale(False)
        plt.scatter(-100, -100, s, [c], marker, label=f'Sil={silhouette:.3f}')

    plt.legend(fontsize=20, markerscale=2.)
    if title:
        # plt.title(f"{title} Plot")
        plt.xlim(xl, xr)
        plt.ylim(yl, yr)
        plt.xlabel(f"{title}-1", fontsize=24)
        plt.ylabel(f"{title}-2", fontsize=24)
        plt.xticks(fontsize=20)
        plt.yticks(fontsize=20)
    plt.tight_layout()
    if fig_path:
        plt.savefig(fig_path)
        plt.close()


def heatmap_plot(values, title, filename, xticks):
    plt.figure()
    plt.title(title)
    plt.ylabel("METHOD")
    plt.xlabel("GENES NUMBER")
    hm = sns.heatmap(values, annot=True)
    hm.set_yticklabels(["DEG", "LFC", "HR", "PCOX", "AE-L", "AE-M", "AE-S", "AE-XS"])
    hm.set_xticklabels(xticks)
    plt.tight_layout()
    plt.savefig(filename)
    plt.close()


def find_best_th_roc_simple(y, y_hat):
    fpr, tpr, thresholds = roc_curve(y, y_hat)
    optimal_idx = np.argmax(tpr - fpr)
    optimal_threshold = thresholds[optimal_idx]
    return optimal_threshold


def find_best_th_min(y, y_hat):
    ths = np.linspace(0.001, 0.999, 198)

    recalls    = np.zeros(ths.shape)
    precisions = np.zeros(ths.shape)

    min_to_beat = 0
    th_to_beat  = 0
    pre_to_beat = 0
    rec_to_beat = 0
    best_index  = 0
    for i, th in enumerate(ths):
        y_hat_l = y_hat >= th
        CM = confusion_matrix(y, y_hat_l)
        (tn, fp, fn, tp) = CM.ravel()
        if tp + fp > 0:
            precision = tp / (tp + fp)
        else:
            precision = 0
        recall = tp / (tp + fn)
        min_pr = min(precision, recall)
        if min_pr > min_to_beat:
            min_to_beat = min_pr
            th_to_beat = th
            pre_to_beat = precision
            rec_to_beat = recall
            best_index = i
    print(f"Best Threshold = {th_to_beat:.3f}")
    print(f"[Precision = {pre_to_beat:.3f}, Recall = {rec_to_beat:.3f}]")

    diff_pre_rec = abs(pre_to_beat - rec_to_beat)
    if diff_pre_rec > 0.1:
        print("Recalls    = ", recalls[best_index-10:best_index+10])
        print("Precisions = ", recalls[best_index-10:best_index+10])
        print("Thresholds = ", ths[best_index-10:best_index+10])
    return th_to_beat


def plot_km(event, time, event_pred, filename=None, show_confidence=True):
    plt.figure(figsize=(10, 6))
    for risk_type in [0, 1]:
        mask_risk = event_pred == risk_type
        event_risk = event[mask_risk]
        time_risk = time[mask_risk]
        print(f"mask_risk: {mask_risk.shape} -- {mask_risk.dtype}, "
              f"event_risk: {event_risk.shape} -- {event_risk.dtype}, "
              f"time_risk: {time_risk.shape} -- {time_risk.dtype} -- {time_risk.min():.1f} -- {time_risk.max():.1f}")
        if len(event_risk) == 0:
            event_risk_km = event
            time_risk_km = time
        else:
            event_risk_km = event_risk
            time_risk_km = time_risk

        if show_confidence:
            time_treatment, survival_prob_treatment, conf_int = kaplan_meier_estimator(
                event_risk_km.astype(bool),
                time_risk_km,
                conf_type="log-log"
            )

            plt.step(time_treatment, survival_prob_treatment, where="post", label=f"Risk = {risk_type}")
            plt.fill_between(time_treatment, conf_int[0], conf_int[1], alpha=0.25, step="post")
        else:
            time_treatment, survival_prob_treatment = kaplan_meier_estimator(
                event_risk_km.astype(bool),
                time_risk_km
            )
            plt.step(time_treatment, survival_prob_treatment, where="post", label=f"Risk = {risk_type}")

    plt.yticks([0, 0.25, 0.50, 0.75, 1.00])
    plt.ylim(0, 1)
    plt.title("Kaplan-Meier Curves")
    plt.ylabel("est. probability of survival $\hat{S}(t)$")
    plt.xlabel("time $t$")
    plt.legend(loc="best")
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()


def plot_coxph(event, time, predictors, filename=None):
    cph = multivariate_surv_analysis_lifelines(event, time, predictors)
    plt.figure(figsize=(10, 10))
    cph.plot()
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()


def plot_trend(x, trends, names, filename=None):
    plt.figure(figsize=(10, 10))
    for trend, name in zip(trends, names):
        plt.plot(x, trend, marker='o', linestyle='dashed', label=name)
    plt.title("AUC varying years of follow-up")
    plt.xlabel("FU [months]")
    plt.ylabel("AUC")
    plt.yticks(np.linspace(0.35, 1.0, 14))
    plt.legend()
    plt.tight_layout()
    if filename:
        plt.savefig(filename)
        plt.close()
