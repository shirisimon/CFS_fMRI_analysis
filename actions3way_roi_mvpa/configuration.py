import numpy as np
from sklearn import linear_model
from sklearn import svm

from sklearn import ensemble
# from sklearn import lda


class Params(object):
    def __init__(self):
        self.collapse_opacities   = False
        self.sublist              = ['379', '381', '382', '383', '384',
                                     '385', '386', '388', '389', '390',
                                     '391', '392', '393', '394', '395']  # TODO: get it from hdf file
        self.tr_pattern           = '_nmsk'
        self.te_pattern           = '_msk'
        self.conds_num            = 3                                    # conditions number in one functional run
        self.do_standard_scaling  = True
        self.do_pca               = False
        self.pca_variance         = 0.75
        self.do_feature_selection = True
        self.voi_num              = 15                                   # TODO: get it from hdf file
        self.sub_num              = 15                                   # TODO: get it from hdf file
        self.fmri_data_type       = 'betas'                              # 'raw' (% signal change) OR 'betas'
        self.tpoi                 = 2                                    # 'Time Point Of Interest' where the signal reach to its peak relative to onset
        self.baseline             = 0
        self.baseline_method      = 'file'                               # epoch based (divide by baseline) ORr file based (just zscore)
        self.snr_factor           = 3                                    # 1 (no snr), 3 (average across 3 trials)
        self.model                = Classifer().svc()



class Classifer():

    def svc(self):
        algo = svm.SVC()
        kernel = ['rbf', 'poly', 'sigmoid']
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

    def lda(self):
        algo = lda.LDA()
        gs_params = [{'solver': ['svd', 'lsqr', 'eigen']}]
        scores = 'accuracy'
        return algo, gs_params, scores

class FeatureSelector():

    def logistic_regression(self):
        return linear_model.LogisticRegression(penalty='l1')



