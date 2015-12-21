import preprocessing as prp
import configuration as cfg
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import accuracy_score
import csv
import cProfile, pstats, StringIO
# pr = cProfile.Profile()
# pr.enable()


def sub_analysis(subs, sub_data_generator):
        for sub in subs:
            X_tr, y_tr, X_te, y_te, pipe = sub_data_generator.next()
            # train on nmsk:
            algo, algo_params, score = prep.params.model
            clf = GridSearchCV(algo, algo_params, cv=9, scoring=score, n_jobs=-1, verbose=1)
            clf.fit(pipe.transform(X_tr), y_tr)
            # test on msk:
            X_te_high, y_te_high, X_te_low, y_te_low = prep.split_high_low_sets(X_te, y_te)
            y_te_high = prep.collapse_lbls(y_te_high, (3,4,5), (0,1,2))
            y_true_high, y_pred_high = y_te_high, clf.predict(pipe.transform(X_te_high))
            y_true_low,  y_pred_low  = y_te_low,  clf.predict(pipe.transform(X_te_low))
            yield clf.best_score_, accuracy_score(y_true_high, y_pred_high), accuracy_score(y_true_low, y_pred_low)

def voi_analysis(vois, subs):
    for voi in vois:
        sub_data_generator = prep.get_sub_data(voi)
        yield zip(*sub_analysis(subs, sub_data_generator))

def save_as_csv(filename, var):
    with open(filename, 'w+') as result:
        writer = csv.writer(result, dialect='excel')
        writer.writerows(var)


if __name__ == '__main__':
    FILE_NAME = 'data2skl.mat'
    prep = prp.Preprocessing(FILE_NAME)
    vois = xrange(prep.voidata.chunks[0])
    subs = xrange(len(cfg.Params().sublist))
    nmsk, msk_high, msk_low = zip(*voi_analysis(vois, subs))
    # save results to csv:
    save_as_csv('nmsk_cv.csv', nmsk)
    save_as_csv('msk_high.csv', msk_high)
    save_as_csv('msk_low.csv', msk_low)





# pr.disable()
# s = StringIO.StringIO()
# sortby = 'cumulative'
# ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
# ps.print_stats()


# TODO: try classifiers