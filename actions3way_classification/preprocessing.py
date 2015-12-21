
import configuration
import h5py
from scipy import stats
import numpy as np
from sklearn.pipeline import Pipeline
from sklearn.preprocessing import StandardScaler
from sklearn import decomposition




class DataLoader(object):
    def __init__(self, path):
        """
        This class deals with whatever related to access to HDF5 file parts
        :param path:
        :return:
        """
        self.params = configuration.Params()
        self.path = path
        self.file = h5py.File(self.path)
        self.voidata  = self.file['voidata']
        self.vtcfiles = self.file['vtcfiles']
        self.prtdata = self.file['prtdata']

    def get_vtc_data(self, vtc_idx, voi_idx):
        return self.file[self.voidata[voi_idx, vtc_idx]][:].transpose()

    def get_prt_data(self, prt_idx):
        """
        get the onsets of each condition from prt
        :param prt_idx:
        :return:
        """
        prt = self.file[self.prtdata[prt_idx,0]]
        conds = []
        for cond in range(self.params.conds_num):
            tmps = self.file[prt[cond, 0]][:][0]
            tmps = [int(i) for i in tmps]
            conds.append(tmps)
        return conds

    def get_vtc_filename(self, vtc_idx):
        s = []
        [s.append(u''.join(unichr(c))) for c in self.file[self.vtcfiles[0,vtc_idx]]]
        return ''.join(s)


class Preprocessing(DataLoader):

    def get_sub_data(self, voi):
        for sub in self.params.sublist:
            X_tr, y_tr = self.get_concat_vols_and_labels(sub+self.params.tr_pattern, voi)
            X_te, y_te = self.get_concat_vols_and_labels(sub+self.params.te_pattern, voi)
            pipe = self.skl_pipeline(X_tr)
            y_tr = self.collapse_lbls(y_tr, (3,4,5), (0,1,2))
            yield X_tr, y_tr, X_te, y_te, pipe

    def get_concat_vols_and_labels(self, vtc_pattern, voi_idx):
        """
        concat relevant volumes and labels from one subject' vtc
        :param vtc_pattern: with subject name and run name
        :param voi_idx:
        :return: X, y
        """
        vtcs_idx = self.find_vtcs_idx(vtc_pattern)
        cvols = []
        clbls = []
        for idx in vtcs_idx:
            vtc = self.scale_vtc(idx, voi_idx)                              # zscore vtc
            vols, lbls = self.extract_vtc_vols_and_labels(idx, vtc)         # extract relevant volumes and labels
            cvols.append(vols)
            clbls.append(lbls)
        cvols = [c for c in np.array(cvols)]
        clbls = [l for l in np.array(clbls)]
        return np.concatenate((cvols)), np.concatenate((clbls))

    def extract_vtc_vols_and_labels(self, vtc_idx, vtc_data):
        """
        extract relevant volumes and labels from each vtc
        :param vtc_idx:
        :param vtc_data:
        :return:
        """
        prt = self.get_prt_data(vtc_idx)
        cvols = []
        clbls = []  # labels
        for cond in xrange(len(prt)):
            vols = self.extract_cond_vols(vtc_data, prt[cond], self.params.tpoi)
            lbls = vols.shape[0] * [cond]
            cvols.append(vols)                                             # concat conditions volumes of one vtc file
            clbls.append(lbls)                                             # concat conditions labels of one vtc file
        cvols = [c for c in np.array(cvols)]
        clbls = [l for l in np.array(clbls)]
        return np.concatenate((cvols)), np.concatenate((clbls))

    def extract_cond_vols(self, vtc_data, onsets, tpoi):
        # TODO (if not working): add function of % signal change
        """
        convert vols by time points
        :param vtc_data:
        :param onsets: onset of events (list)
        :param tpoi: time point of interest
        :return:
        """
        cond_tpoi = [x+tpoi for x in onsets]
        if cond_tpoi[-1]>vtc_data.shape[0]-1:
            cond_tpoi = cond_tpoi[:-1]
        return vtc_data[cond_tpoi, :]

    def scale_vtc(self, vtc_idx, voi_idx):
        mat = self.get_vtc_data(vtc_idx, voi_idx)
        return stats.zscore(mat)

    def find_vtcs_idx(self, vtc_pattern):
        idxlist = []
        for idx in xrange(self.vtcfiles.shape[1]):
            if vtc_pattern in self.get_vtc_filename(idx):
                idxlist.append(idx)
        return idxlist

    def skl_pipeline(self, X_train):
        """
        create preprocessing pipeline and fit to X
        :param X_train:
        :return: pipeline object and fitted X
        """
        scaler  = StandardScaler()
        if self.params.do_pca:
            print("reducing dimensions with PCA")
            pca = decomposition.KernelPCA(kernel='rbf')
            pipe = Pipeline(steps=[('scale', scaler), ('pca', pca)]).fit(X_train)
        else:
            pipe = Pipeline(steps=[('scale', scaler)]).fit(X_train)
        return pipe

    def collapse_lbls(self, y, orig_lbls, target_lbls):
        """
        convert certain labels in y to other labels
        :param y: labels list
        :param orig_lbls: original labels (tuple)
        :param target_lbls: target labels (tuple)
        :return: y converted
        """
        for idx, item in enumerate(y):
            if item in orig_lbls:
                y[idx] = target_lbls[orig_lbls.index(item)]
        return np.asarray(y)

    def split_high_low_sets(self, X, y):
        """
        split high and low samples to different data sets (should be applied to test set)
        :param X:
        :param y:
        :return:2*X and 2*y sets
        """
        X_high = X[np.logical_or(np.logical_or(y==3, y==4), y==5)]
        y_high = y[np.logical_or(np.logical_or(y==3, y==4), y==5)]
        X_low  = X[np.logical_or(np.logical_or(y==0, y==1), y==2)]
        y_low  = y[np.logical_or(np.logical_or(y==0, y==1), y==2)]
        return X_high, y_high, X_low, y_low