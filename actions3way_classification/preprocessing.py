
import configuration
import h5py
from scipy import stats
import numpy as np
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
    # TODO: add print reports and timing
    def do_preprocessing(self):
        voisub_X_tr = [[None]*self.params.sub_num]*self.params.voi_num
        voisub_y_tr = [[None]*self.params.sub_num]*self.params.voi_num
        voisub_X_te = [[None]*self.params.sub_num]*self.params.voi_num
        voisub_y_te = [[None]*self.params.sub_num]*self.params.voi_num
        scale_mdl = [[None]*self.params.sub_num]*self.params.voi_num
        pca_mdl   = [[None]*self.params.sub_num]*self.params.voi_num
        vidx = 0 # voi idx
        sidx = 0 # sub idx
        for voi in xrange(self.voidata.chunks[0]):
            for sub in self.params.sublist:
                # concat each voi-sub to vols train and test sets
                voisub_X_tr[vidx][sidx], voisub_y_tr[vidx][sidx] = \
                    self.get_concat_vols_and_labels(sub+self.params.tr_pattern, voi)
                voisub_X_te[vidx][sidx], voisub_y_te[vidx][sidx] = \
                    self.get_concat_vols_and_labels(sub+self.params.te_pattern, voi)
                # skl preprocessing of tr set:
                scaler, voisub_X_tr[vidx][sidx] = self.skl_scaling(voisub_X_tr[vidx][sidx])
                pca, voisub_X_tr[vidx][sidx] = self.skl_pca(voisub_X_tr[vidx][sidx])
                scale_mdl[vidx][sidx] = scaler    # save the scaling model for each voi-sub
                pca_mdl[vidx][sidx] = pca         # save the pca model for each voi-sub
                vidx += 1
                sidx += 1
        return vidx, sidx
        # save data:
        #

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
            vtc = self.scale_vtc(idx, voi_idx)                     # zscore vtc
            vols, lbls = self.extract_vtc_vols(idx, vtc)           # extract relevant volumes and labels
            cvols.append(vols)
            clbls.append(lbls)
        cvols = np.array(cvols)
        clbls = np.array(clbls)# concat vtc volumes of one subject
        vshape = (cvols.shape[0]*cvols.shape[1], cvols.shape[2])
        lshape = (clbls.shape[0]*clbls.shape[1])
        return np.reshape(cvols, vshape), np.reshape(clbls, lshape)

    def extract_vtc_vols(self, vtc_idx, vtc_data):
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
            vols = vtc_data[prt[cond],:]
            lbls = len(prt[cond]) * [cond]
            cvols.append(vols)                                     # concat conditions volumes of one vtc file
            clbls.append(lbls)                                     # concat conditions labels of one vtc file
        cvols = np.array(cvols)
        clbls = np.array(clbls)
        vshape = (cvols.shape[0]*cvols.shape[1], cvols.shape[2])
        lshape = (clbls.shape[0]*clbls.shape[1])
        return np.reshape(cvols, vshape), np.reshape(clbls, lshape)

    def scale_vtc(self, vtc_idx, voi_idx):
        mat = self.get_vtc_data(vtc_idx, voi_idx)
        return stats.zscore(mat)

    def find_vtcs_idx(self, vtc_pattern):
        idxlist = []
        for idx in xrange(self.vtcfiles.shape[1]):
            if vtc_pattern in self.get_vtc_filename(idx):
                idxlist.append(idx)
        return idxlist

    def skl_scaling(self, X_train):
        scaler  = StandardScaler().fit(X_train)                   # features scaling (should be done only on X train)
        scaled_X_train = scaler.transform(X_train)
        return scaler, scaled_X_train

    def skl_pca(self, X_train):
        if self.params.do_pca:
            print("reducing dimensions with PCA")
            pca = decomposition.KernelPCA(kernel='rbf').fit(X_train)
            decomposed_X_train = pca.transform(X_train)
            return pca, decomposed_X_train
        else:
            return None, X_train