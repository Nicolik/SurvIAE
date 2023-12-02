import numpy as np
import pandas as pd
from sksurv.linear_model import CoxPHSurvivalAnalysis
from lifelines import CoxPHFitter


def univariate_surv_analysis(event, time, predictors):
    structured_event_time = array2structured(event, time)
    n_features = predictors.shape[1]
    scores = np.empty(n_features)
    m = CoxPHSurvivalAnalysis()
    for j, feat_name in enumerate(predictors.columns):
        Xj = np.expand_dims(predictors[feat_name], 1)
        m.fit(Xj, structured_event_time)
        scores[j] = m.score(Xj, structured_event_time)
    print(pd.Series(scores, index=predictors.columns))
    return scores


def multivariate_surv_analysis(event, time, predictors):
    structured_event_time = array2structured(event, time)
    estimator = CoxPHSurvivalAnalysis()
    estimator.fit(predictors, structured_event_time)
    # print(estimator.get_params())
    print(pd.Series(estimator.coef_, index=predictors.columns))
    return estimator


def multivariate_surv_analysis_lifelines(event, time, predictors):
    event_df = pd.DataFrame(data={'Event': event})
    event_df.index = predictors.index
    time_df = pd.DataFrame(data={'Time': time})
    time_df.index = predictors.index
    data = pd.concat([event_df, time_df, predictors], axis=1)

    cph = CoxPHFitter()
    cph.fit(data, duration_col='Time', event_col='Event')
    # print(cph.print_summary())
    return cph


def array2structured(event, time):
    aux = [(e, t) for e, t in zip(event, time)]
    return np.array(aux, dtype=[('Status', '?'), ('Survival_in_days', '<f8')])
