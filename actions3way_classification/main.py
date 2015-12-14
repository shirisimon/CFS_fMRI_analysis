

import preprocessing as prp
import configuration as cfg

# TODO: if ds is not exist, generate it:
FILE_NAME = 'data2skl.mat'
prep = prp.Preprocessing(FILE_NAME)
X_ds,y_ds = prep.do_preprocessing()
# TODO: else load it

