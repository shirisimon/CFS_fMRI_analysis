
import configuration
import h5py
from scipy import stats
import numpy as np


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

    def get_mat(self, vtc_idx, voi_idx):
        return self.file[self.voidata[voi_idx, vtc_idx]][:].transpose()

    def get_vtc_filename(self, vtc_idx):
        s = []
        [s.append(u''.join(unichr(c))) for c in self.file[self.vtcfiles[0,vtc_idx]]]
        return ''.join(s)


class Preprocessing(DataLoader):
    def do_preprocessing(self):
        for voi in xrange(self.voidata.chunks[0]):
            for sub in self.params.sublist:
                # concat each voi-sub to train and test sets
                voisubX_tr = self.concat_vtcs(sub+self.params.tr_pattern, voi)
                voisubX_te = self.concat_vtcs(sub+self.params.te_pattern, voi)

    def concat_vtcs(self, vtc_pattern, voi_idx):
        vtcs_idx = self.find_vtcs_idx(vtc_pattern)
        concat_mat = []   # TODO: insert the size of the voi
        for idx in vtcs_idx:
            mat = self.scale_vtc(idx, voi_idx)
            concat_mat.append(mat)
        concat_mat = np.array(concat_mat)
        shape = (concat_mat.shape[0]*concat_mat.shape[1], concat_mat.shape[2])
        return np.reshape(concat_mat, shape)


    def scale_vtc(self, vtc_idx, voi_idx):
        mat = self.get_mat(vtc_idx, voi_idx)
        return stats.zscore(mat)

    def find_vtcs_idx(self, vtc_pattern):
        idxlist = []
        for idx in range(self.vtcfiles.shape[1]):
            if vtc_pattern in self.get_vtc_filename(idx):
                idxlist.append(idx)
        return idxlist