import numpy as np
import preprocessing as prp
import configuration as cfg
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import accuracy_score
import csv
import cProfile, pstats, StringIO
pr = cProfile.Profile()
pr.enable()


def sub_analysis(subs, sub_data_generator):
        for sub in subs:
            print("process sub: %0.0f" % sub)
            # X_tr, y_tr, X_te, y_te, pipe = sub_data_generator.next()
            X, y, _, _, pipe = sub_data_generator.next()
            cv_scores  = []
            te1_scores = []
            te2_scores = []
            for state in range(3): #
                X_tr, X_te, y_tr, y_te = prep.stratified_split(X, y, state)
                # train on nmsk:
                algo, algo_params, score = prep.params.model
                clf = GridSearchCV(algo, algo_params, cv=3, scoring=score, n_jobs=-1, verbose=0)
                clf.fit(prep.convert(pipe, X_tr), y_tr)
                # test on msk:
                if prep.params.collapse_opacities:
                    X_te_high, y_te_high, X_te_low, y_te_low = prep.split_high_low_sets(X_te, y_te)
                    y_te_high = prep.collapse_lbls(y_te_high, (3,4,5), (0,1,2))
                    y_true_high, y_pred_high = y_te_high, clf.predict(prep.convert(pipe, X_te_high))
                    y_true_low,  y_pred_low  = y_te_low,  clf.predict(prep.convert(pipe, X_te_low))
                    cv_scores.append(clf.best_score_)
                    te1_scores.append(accuracy_score(y_true_high, y_pred_high))
                    te2_scores.append(accuracy_score(y_true_low, y_pred_low))
                else:
                    y_true, y_pred = y_te, clf.predict(prep.convert(pipe, X_te))
                    cv_scores.append(clf.best_score_)
                    te1_scores.append(accuracy_score(y_true, y_pred))
                    te2_scores = None
            yield cv_scores, te1_scores, te2_scores

def voi_analysis(vois, subs):
    for voi in vois:
        print("process voi: %0.0f" % voi)
        sub_data_generator = prep.get_sub_data(voi)
        yield zip(*sub_analysis(subs, sub_data_generator))

def save_as_csv(filename, var):
    with open(filename, 'w+') as result:
        writer = csv.writer(result, dialect='excel')
        writer.writerows(var)


if __name__ == '__main__':
    FILE_NAME = 'data2skl_high.mat'
    prep = prp.Preprocessing(FILE_NAME)
    vois = xrange(prep.voidata.chunks[0]-1)
    subs = xrange(len(cfg.Params().sublist))
    if prep.params.collapse_opacities:
        nmsk, msk_high, msk_low = zip(*voi_analysis(vois, subs))
        save_as_csv('results_msk_high.csv', msk_high)
        save_as_csv('results_msk_low.csv', msk_low)
    else:
        nmh_cv, nmh_test, _ = zip(*voi_analysis(vois, subs))
        print np.mean(nmh_cv), np.mean(nmh_test)
        print np.median(nmh_cv), np.median(nmh_test)
        # save_as_csv('results_NMH.csv', msk)
        # save_as_csv('results_NMH_cv.csv', nmsk)


pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()


# TRIED:
# with smoothed, with snr, with lasso feature selection
# pca optimization - (best for v1) is with 0.75 variance
# only zscore and zscore + epoch based (% signal change)
# classifers - svm rbf, random forest, lda (none is best)
# separate and collapsed high and low