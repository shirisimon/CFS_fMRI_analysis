import preprocessing as prp
import configuration as cfg
import pickle
import scipy.io as io



def sub_analysis(subs, sub_data_generator):
    for sub in subs:
        print("process sub: %0.0f" % sub)
        X_nm, y_nm, X_ma, y_ma, pipe = sub_data_generator.next()
        yield X_nm, y_nm, X_ma, y_ma

def voi_analysis(vois, subs):
    for voi in vois:
        print("process voi: %0.0f" % voi)
        sub_data_generator = prep.get_sub_data(voi)
        yield zip(*sub_analysis(subs, sub_data_generator))


if __name__ == '__main__':
    FILE_NAME = 'data2skl_high.mat'
    prep = prp.Preprocessing(FILE_NAME)
    vois = xrange(prep.voidata.chunks[0]-1)
    subs = xrange(len(cfg.Params().sublist))

    if prep.params.generate_data:
        X_nm, y_nm, X_ma, y_ma = zip(*voi_analysis(vois, subs))
        pickle.dump([X_nm, y_nm, X_ma, y_ma], open("data2rsa.p", "wb"))
    else:
        dat = pickle.load(open("data2rsa.p", "rb"))
        X_nm = dat[0]
        y_nm = dat[1]
        X_ma = dat[2]
        y_ma = dat[3]

    for voi in vois:
        for sub in subs:
            filename_nm  = prep.params.voilist[voi]+'_s'+prep.params.sublist[sub]+'_'+'nmh'
            io.savemat(filename_nm, {'data': X_nm[voi][sub], 'labels': y_nm[voi][sub]})
            filename_ma  = prep.params.voilist[voi]+'_s'+prep.params.sublist[sub]+'_'+'mah'
            io.savemat(filename_ma, {'data': X_ma[voi][sub], 'labels': y_ma[voi][sub]})



# TRIED:
# with smoothed, with snr, with lasso feature selection
# pca optimization - (best for v1) is with 0.75 variance
# only zscore and zscore + epoch based (% signal change)
# classifers - svm rbf, random forest, lda (none is best)
# separate and collapsed high and low