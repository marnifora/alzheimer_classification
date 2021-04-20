import numpy as np
import os
from sklearn.ensemble import RandomForestClassifier
from sklearn.metrics import roc_auc_score
import matplotlib.pyplot as plt
import argparse
from scipy.stats import pearsonr


def get_nodata_stats(matrix, outfile):
    print('Get nodata stats, input matrix: {}'.format(matrix.shape))
    num_pat, num_snps = matrix.shape
    print('Number of patients: {}, number of SNPs: {}'.format(num_pat, num_snps))
    nodata = np.apply_along_axis(lambda row: sum(row == -1), 0, matrix)
    print(nodata[:10])
    print(len(nodata))
    nodata_ratio = nodata / num_pat
    print(nodata_ratio[:10])
    print('Nodata stats: {}'.format(nodata_ratio.shape))
    with open(outfile, 'w') as f:
        f.write('\n'.join([str(round(float(el), 6)) for el in nodata_ratio]))
    print('Nodata stats saved to {}'.format(outfile))
    return nodata_ratio


def get_importance(X, y, X_test, y_test, outfile):
    rf = RandomForestClassifier(n_estimators=500)
    print('Fitting random forest')
    rf.fit(X, y)
    importances = rf.feature_importances_
    print('Importances of SNPs: {}'.format(importances.shape))
    prob = rf.predict_proba(X_test)
    order = [i for i, c in enumerate(rf.classes_) if c == 1]
    y_score = prob[:, order]
    train_score, test_score, auc = rf.score(X, y), rf.score(X_test, y_test), roc_auc_score(y_test, y_score)
    with open(outfile, 'w') as f:
        f.write('{:.3f}\t{:.3f}\t{:.3f}\n'.format(train_score, test_score, auc))
        f.write('\n'.join([str(round(float(el), 6)) for el in importances]))
    print('Importances saved to {}'.format(outfile))
    return importances


def main(datadir, run, perc, chrlist=None, stats=False, forest=False, plot=False, title=''):

    if stats or forest:
        if len(chrlist) < 23:
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

    stats_file = os.path.join(datadir, 'nodata_stats_{}_{}.txt'.format(perc, run))
    if stats:
        nodata = get_nodata_stats(np.concatenate([X_train, X_test], axis=0), stats_file)
    elif plot:
        nodata = [float(el) for el in open(stats_file, 'r').read().strip().split('\n')]
        print('Nodata stats loaded from {}, len {}'.format(stats_file, len(nodata)))

    imp_file = os.path.join(datadir, 'importances_{}_{}.txt'.format(perc, run))
    if forest:
        importance = get_importance(X_train, y_train, X_test, y_test, imp_file)
    elif plot:
        importance = open(imp_file, 'r').read().strip().split('\n')
        train_score, test_score, auc = [float(el) for el in importance[0].split('\t')]
        importance = [float(el) for el in importance[1:]]
        print('Importances loaded from {}, len {}'.format(imp_file, len(importance)))

    if plot:
        outfile = os.path.join(datadir, 'nodata_vs_importance_plot_{}_{}.png'.format(perc, run))
        print('Calculating of Pearson correlation')
        corr, pvalue = pearsonr(nodata, importance)
        nodata2 = sum(np.array(nodata) > 0.02) / len(nodata) * 100
        nodata5 = sum(np.array(nodata) > 0.05) / len(nodata) * 100
        title = '{}\nAll SNPs: {}, lack of data is >0.02 for {:.2f}%; >0.05 for {:.2f}%' \
                '\nPearson corr: {:.3f}\n' \
                'Random forest stats: train score: {:.3f}, test score: {:.3f}, AUC: {:.3f}'.\
            format(title, len(nodata), nodata2, nodata5, corr, train_score, test_score, auc)
        fig, ax = plt.subplots(figsize=(10, 8))
        plt.scatter(nodata, importance, marker='.', s=1)
        plt.title(title, fontsize=10)
        plt.xlabel('No-data ratio')
        plt.ylabel('Importance')
        plt.savefig(outfile)
        plt.show()


parser = argparse.ArgumentParser('Create scatter plot: nodata ratio of SNPs vs their importance in the classifier')
parser.add_argument('--datadir', type=str, default='', help='Directory with files to use.')
parser.add_argument('--chr', type=str, default='1-23', help='Range of chromosomes to plot')
parser.add_argument('--run', type=int, default=0, help='Run number of analysis to plot')
parser.add_argument('--perc', type=int, default=90, help='Perc value of analysis to plot')
parser.add_argument('--stats', action='store_true', help='Calculate nodata stats')
parser.add_argument('--forest', action='store_true', help='Train random forest and get importances of SNPs')
parser.add_argument('--plot', action='store_true', help='Create scatter plot')
parser.add_argument('--title', type=str, default='Scatter plot', help='Title of a scatter plot')
args = parser.parse_args()

fr, to = [int(el) for el in args.chr.split('-')]
chrlist = [el for el in range(fr, to+1)]
main(args.datadir, args.run, args.perc, chrlist, args.stats, args.forest, args.plot, args.title)
