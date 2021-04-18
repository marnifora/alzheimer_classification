import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
import matplotlib.pyplot as plt
import argparse


def get_nodata_stats(matrix, outfile):
    print('Get nodata stats, input matrix: {}'.format(matrix.shape))
    num_genes = matrix.shape[1]
    print('Number of SNPs: {}'.format(num_genes))
    nodata_ratio = np.apply_along_axis(lambda row: sum(row == -1), 0, matrix) / num_genes
    print('Nodata stats: {}'.format(nodata_ratio.shape))
    with open(outfile, 'w') as f:
        f.write('\n'.join([str(round(float(el), 6)) for el in nodata_ratio]))
    print('Nodata stats saved to {}'.format(outfile))
    return nodata_ratio


def get_importance(X, y, X_test, y_test):
    rf = RandomForestClassifier(n_estimators=500)
    print('Fitting random forest')
    rf.fit(X, y)
    importances = rf.feature_importances_()
    print('Importances of SNPs: {}'.format(importances.shape))
    prob = rf.predict_proba(X_test)
    order = [i for i, c in enumerate(rf.classes_) if c == 1]
    y_score = prob[:, order]
    return importances, rf.score(X, y), rf.score(X_test, y_test), roc_auc_score(y_test, y_score)


def main(datadir, run, perc, chrlist=None, stats=False, forest=False, plot=False):

    if chrlist is not None:
        ch = chrlist[0]
        print('Loading X matrix for chr {}'.format(ch))
        X_train = np.load(os.path.join(datadir, 'X_train_chr{}_{}_{}.npy'.format(ch, perc, run)))
        X_test = np.load(os.path.join(datadir, 'X_test_chr{}_{}_{}.npy'.format(ch, perc, run)))
        for ch in chrlist[1:]:
            print('Loading X matrix for chr {}'.format(ch))
            X_train = np.concatenate([X_train, np.load(os.path.join(datadir, 'X_train_chr{}_{}_{}.npy'.
                                                                   format(ch, perc, run)))], axis=1)
            X_test = np.concatenate([X_test, np.load(os.path.join(datadir, 'X_test_chr{}_{}_{}.npy'.
                                                                   format(ch, perc, run)))], axis=1)
    else:
        print('Loading X matrix for whole genome')
        X_train = np.load(os.path.join(datadir, 'X_train_genome_{}_{}.npy'.format(perc, run)))
        X_test = np.load(os.path.join(datadir, 'X_test_genome_{}_{}.npy'.format(perc, run)))

    print('Loading y matrix for whole genome')
    y_train = np.load(os.path.join(datadir, 'y_train_genome_{}_{}.npy'.format(perc, run)))
    y_test = np.load(os.path.join(datadir, 'y_test_genome_{}_{}.npy'.format(perc, run)))

    if stats:
        outfile = os.path.join(datadir, 'nodata_stats_{}_{}.txt'.format(perc, run))
        nodata = get_nodata_stats(np.concatenate([X_train, X_test], axis=0), outfile)
    if forest:
        importance, train_score, test_score, auc = get_importance(X_train, y_train, X_test, y_test)

    if plot:
        print('Nodata: {}, importances: {}'.format(nodata.shape, importance.shape))
        print('Train score: {}, test score: {}, AUC: {}'.format(train_score, test_score, auc))
        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(nodata, importance)
        plt.savefig(outfile)
        plt.close()


parser = argparse.ArgumentParser('Create scatter plot: nodata ratio of SNPs vs their importance in the classifier')
parser.add_argument('--datadir', type=str, default='', help='Directory with files to use.')
parser.add_argument('--chr', type=str, default='1-23', help='Range of chromosomes to plot')
parser.add_argument('--run', type=int, default=0, help='Run number of analysis to plot')
parser.add_argument('--perc', type=int, default=90, help='Perc value of analysis to plot')
parser.add_argument('--stats', action='store_true', help='Calculate nodata stats')
parser.add_argument('--forest', action='store_true', help='Train random forest and get importances of SNPs')
parser.add_argument('--plot', action='store_true', help='Create scatter plot')
args = parser.parse_args()

fr, to = [int(el) for el in args.chr.split('-')]
chrlist = [el for el in range(fr, to+1)]
main(args.datadir, args.run, args.perc, chrlist, args.stats, args.forest, args.plot)
