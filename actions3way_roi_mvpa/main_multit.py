import numpy as np
import preprocessing as prp
import configuration as cfg
import csv
import time
import winsound
import cProfile, pstats, StringIO

from termcolor import colored
pr = cProfile.Profile()
pr.enable()


def sub_analysis(voi, subs, sub_data_generator):
        for sub in subs:
            print("process voi {:02.0f} | sub {:02.0f}".format(voi+1, sub+1))
            X, y, _, _, pipe = sub_data_generator.next()
            X = prep.convert(pipe, X)
            treal = mult.multit_real(X,y)
            if np.isnan(treal):
                raise Warning("treal is nan")
            tshuffle_perms = mult.multit_shuffle(X,y, n_perms=100)
            pval = mult.calc_pval(treal, tshuffle_perms)
            yield treal, tshuffle_perms, pval

def voi_analysis(vois, subs):
    for voi in vois:
        sub_data_generator = prep.get_sub_data(voi)
        treals, tshuffles_perms, sub_pvals = zip(*sub_analysis(voi, subs, sub_data_generator))
        print colored ("done process voi {:02.0f} in {:02.0f}s".format(voi, time.time()-t), 'cyan')
        group_treal = np.mean(treals)
        group_tshuffle = mult.tshuffle_bootstrap(tshuffles_perms, n_perms=100) # TODO: write this function
        group_pval = mult.calc_pval(group_treal, group_tshuffle)
        yield sub_pvals, group_pval

def save_as_csv(filename, var):
    with open(filename, 'w+') as result:
        writer = csv.writer(result, dialect='excel')
        writer.writerows(var)

def termination_alert():
    for i in range(3):
        winsound.Beep(1500,250)
        time.sleep(0.1)
        winsound.Beep(1500,250)




if __name__ == '__main__':
    FILE_NAME = 'data2skl_high.mat'
    prep = prp.Preprocessing(FILE_NAME)
    mult = prp.MultiT()
    vois = xrange(prep.voidata.chunks[0]-1)
    subs = xrange(len(cfg.Params().sublist))
    t = time.time()
    if prep.params.collapse_opacities:
        pass
    else:
        sub_pvals, group_pval = zip(*voi_analysis(vois, subs))
        termination_alert()
        save_as_csv('results_NMH_multit.csv', sub_pvals)
        print group_pval


pr.disable()
s = StringIO.StringIO()
sortby = 'cumulative'
ps = pstats.Stats(pr, stream=s).sort_stats(sortby)
ps.print_stats()


# TRIED:
# with smoothed, with snr, with lasso feature selection
# pca optimization - (best for v1) is with 0.75 variance
# only zscore and zscore + epoch based (% signal change)
# classifers - svm rbf, random forest, lda (none is best)
# separate and collapsed high and low