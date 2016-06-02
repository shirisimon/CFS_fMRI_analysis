# leave subject out

import numpy as np
import preprocessing as prp
import configuration as cfg
import pickle
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import accuracy_score


def sub_analysis(subs, sub_data_generator):
    for sub in subs:
        print("process sub: %0.0f" % sub)
        X, y, _, _, _ = sub_data_generator.next()
        X, y = avg_sub_trials(X, y)
        yield X, y

def avg_sub_trials(X, y):
    Xavg = []
    for c in list(np.unique(y)):
        Xavg.append(np.mean(X[np.where(y==c)[0],:], axis=0))
    return np.vstack(Xavg), np.unique(y)

def voi_analysis(vois, subs):
    for voi in vois:
        print("process voi: %0.0f" % voi)
        sub_data_generator = prep.get_sub_data(voi)
        yield zip(*sub_analysis(subs, sub_data_generator))

def concat_subs(subs, vois, data):
    for voi in vois:
        X, y = [], []
        for sub in subs:
            X.append(data[0][voi][sub])
            y.append(data[1][voi][sub])
        yield np.vstack(X),np.hstack(y)

def run_ml(X, y):
    # do machine learning
    for voi in range(params.voi_num-1):
        te_scores = []
        for cv in range(params.cv_test_num):
            print("process voi: %0.0f | cv %0.0f" % (voi, cv))
            X_tr, X_te, y_tr, y_te = prep.chunks_split(X[voi], y[voi], 14, cv)
            pipel = prep.skl_pipeline(X_tr, y_tr)
            algo, algo_params, score = prep.params.model
            clf = GridSearchCV(algo, algo_params, cv=3, scoring=score, n_jobs=-1, verbose=0)
            clf.fit(prep .convert(pipel, X_tr), y_tr)
            y_true, y_pred = y_te, clf.predict(prep.convert(pipel, X_te))
            te_scores.append(accuracy_score(y_true, y_pred))
        yield np.mean(np.asarray(te_scores)), np.std(np.asarray(te_scores))


if __name__ == '__main__':
    FILE_NAME = 'data2skl_high.mat'
    prep = prp.Preprocessing(FILE_NAME)
    params = cfg.Params()
    vois = xrange(prep.voidata.chunks[0]-1)
    subs = xrange(len(cfg.Params().sublist))

    # get data
    if params.generate_data:
        X, y  = zip(*voi_analysis(vois, subs))
        pickle.dump([X, y], open("lso_data.p", "wb"))

    data = pickle.load(open("lso_data.p", "rb"))
    X, y = zip(*concat_subs(subs, vois, data))

    # do machine learning
    voi_scores_mean, voi_scores_std = zip(*run_ml(X, y))
    pickle.dump([voi_scores_mean, voi_scores_std], open("lso_results.p", "wb"))

