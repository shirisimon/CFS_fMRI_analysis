import numpy as np
from sklearn import svm
from sklearn import ensemble

class Params(object):
    def __init__(self):
        self.collapse_opacities = False
        self.sublist         = ['379', '381', '382', '383', '384',
                                '385', '386', '388', '389', '390',
                                '391', '392', '393', '394', '395']  # TODO: get it from hdf file
        self.tr_pattern      = '_nmsk'
        self.te_pattern      = '_msk'
        self.conds_num       = 3                                    # conditions number in one functional run
        self.do_pca          = True
        self.voi_num         = 15                                   # TODO: get it from hdf file
        self.sub_num         = 15                                   # TODO: get it from hdf file
        self.tpoi            = 2                                    # 'Time Point Of Interest' where the signal reach to its peak relative to onset
        self.baseline        = -1
        self.baseline_method = 'epoch'                              # epoch based (divide by baseline) or file based (just zscore)
        self.model           = Classifer().svc()



class Classifer():

    def svc(self):
        algo = svm.SVC()
        kernel = ['rbf']
        C_range = np.logspace(-3, 3, 7)
        gamma_range = np.logspace(-3, 3, 7)
        # degree_range = [2]
        gs_params = dict(kernel=kernel, C=C_range, gamma=gamma_range) # degree=degree_range)
        scores = 'accuracy'
        return algo, gs_params, scores

    def random_forest(self):
        algo = ensemble.RandomForestClassifier()
        gs_params = [{'n_estimators': [1000], 'criterion': ['gini', 'entropy']}]
        scores = 'accuracy'
        return algo, gs_params, scores
