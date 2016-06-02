import numpy as np
from sklearn import linear_model
from sklearn import svm
from sklearn import mixture

from sklearn import ensemble
from sklearn import decomposition
from sklearn.feature_selection import SelectFromModel


class Params(object):
    def __init__(self):
        self.generate_data          = True
        self.collapse_opacities     = False
        self.voilist                = ['SMA', 'lSMC','lM1', 'rpMdsup',
                                       'rpMdinf', 'rTPJ', 'rIPS', 'rpMv',
                                       'rIFG', 'rpreSMA', 'lpreSMA',
                                       'lTPJ', 'rV1', 'lMT', 'rMT']
        self.sublist                = ['379', '381', '382', '383', '384',
                                       '385', '386', '388', '389', '390',
                                       '391', '392', '393', '394', '395']  # TODO: get it from hdf file
        self.tr_pattern             = '_nmsk'
        self.te_pattern             = '_msk'
        self.conds_num              = 3                                    # conditions number in one functional run
        self.voi_num                = 15                                   # TODO: get it from hdf file
        self.sub_num                = 15                                   # TODO: get it from hdf file
        self.fmri_data_type         = 'raw'                                # 'raw' (% signal change) OR 'betas'

        self.do_standard_scaling    = True
        self.do_pca                 = False
        self.pca_variance           = 50                                  # try: 1, 0.95, 0.75, 50, 25
        self.pca_type               = decomposition.KernelPCA(n_components=self.pca_variance)                           # KernelPCA(kernel='rbf'), FactorAnalysis(), VarianceThreshold()
        self.do_feature_selection   = False
        self.feature_selection_type = SelectFromModel(FeatureSelector().logistic_regression())
        self.tpoi                   = 2                                    # 'Time Point Of Interest' where the signal reach to its peak relative to onset
        self.baseline               = 0
        self.baseline_method        = 'epoch'                               # epoch based (divide by baseline) OR file based (just zscore)
        self.snr_factor             = 1                                    # 1 (no snr), 3 (average across 3 trials)
        self.model                  = Classifer().svc()

        self.cv_test_num            = 15



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

    def gmm(self):
        algo = mixture.GMM()
        gs_params = [{'n_components': [25,50,100], 'covariance_type':['spherical', 'tied', 'diag', 'full']}]
        scores = 'accuracy'
        return algo, gs_params, scores

    # def lda(self):
    #     algo = lda.LDA()
    #     gs_params = [{'solver': ['svd', 'lsqr', 'eigen']}]
    #     scores = 'accuracy'
    #     return algo, gs_params, scores

class FeatureSelector():

    def logistic_regression(self):
        return linear_model.LogisticRegression(penalty='l1')



