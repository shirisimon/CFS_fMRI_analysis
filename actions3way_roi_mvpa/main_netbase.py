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
        yield np.mean(X, axis=1), y

def voi_analysis(vois, subs):
    for voi in vois:
        print("process voi: %0.0f" % voi)
        sub_data_generator = prep.get_sub_data(voi)
        yield zip(*sub_analysis(subs, sub_data_generator))

def concat_vois(subs, vois, data):
    for sub in subs:
        X, y, pipe = [], [], []
        for voi in vois:
            X.append(data[0][voi][sub])
        y.append(data[1][0][sub])
        yield np.asarray(X).T, np.asarray(y).T

def run_ml(X, y):
    # do machine learning
    for sub in range(params.sub_num):
        te_scores = []
        for cv in range(params.cv_test_num):
            print("process sub: %0.0f | cv %0.0f" % (sub, cv))
            X_tr, X_te, y_tr, y_te = prep.stratified_split(X[sub], y[sub], cv)
            y_tr = np.resize(y_tr, (y_tr.shape[0]))
            y_te = np.resize(y_te, (y_te.shape[0]))
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
        pickle.dump([X, y], open("netbase_data.p", "wb"))

    data = pickle.load(open("netbase_data.p", "rb"))
    X, y = zip(*concat_vois(subs, vois, data))

    # do machine learning
    sub_scores_mean, sub_scores_std = zip(*run_ml(X, y))
    pickle.dump([sub_scores_mean, sub_scores_std] , open("netbase_results.p", "wb"))


