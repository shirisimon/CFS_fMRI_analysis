import numpy as np
import preprocessing as prp
import configuration as cfg
from sklearn.grid_search import GridSearchCV
from sklearn.metrics import accuracy_score




def voi_analysis(vois, subs):
    for voi in vois:
        print("process voi: %0.0f" % voi)
        sub_data_generator = prep.get_sub_data(voi)
        yield zip(*sub_analysis(subs, sub_data_generator))

def sub_analysis(subs, sub_data_generator):
    for sub in subs:
        X, y, _, _, pipe = sub_data_generator.next()
        X_tr, X_te, y_tr, y_te = prep.stratified_split(X, y, 0)
        yield np.mean(X_tr, axis=1), np.mean(X_te, axis=1), y_tr, y_te


def concat_vois(subs, vois):
    for sub in subs:
        X_tr_new = []
        for voi in vois:
            X_tr_new.append(X_tr[voi][sub])
        yield np.asarray(X_tr_new).T



if __name__ == '__main__':
    FILE_NAME = 'data2skl_high.mat'
    prep = prp.Preprocessing(FILE_NAME)
    vois = xrange(prep.voidata.chunks[0]-1)
    subs = xrange(len(cfg.Params().sublist))

    X_tr, X_te, y_tr, y_te = zip(*voi_analysis(vois, subs))
    a = zip(*concat_vois(subs, vois))

