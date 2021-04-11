import numpy as np
import pandas as pd
import os
import matplotlib.pyplot as plt
import argparse
import json


def find_weak(ch, indir):
    print('chr {}'.format(ch))
    x = np.genfromtxt(os.path.join(indir, 'matrices/X_chr{}_nodif.csv'.format(ch)), delimiter=',', skip_header=1,
                      dtype=np.int32)
    x = x[:, 1:]  # skip first column
    print('matrix loaded: {}'.format(x.shape))
    weak_ratio = np.apply_along_axis(lambda row: sum(row == -1) / len(row), 0, x)
    print('weak ratio: {}'.format(weak_ratio.shape))
    median = np.median(weak_ratio)
    q1 = np.quantile(weak_ratio, 0.25)
    q3 = np.quantile(weak_ratio, 0.75)
    iqr = q3 - q1
    down_whisker_mask = weak_ratio > (q1 - (1.5 * iqr))
    fliers = []
    num_all_fliers = 0
    if down_whisker_mask.any():
        down_whisker = min(weak_ratio[down_whisker_mask])
        down_fliers = [el for el in weak_ratio if el < down_whisker]
        num_all_fliers += len(down_fliers)
        down_fliers = list(set(down_fliers))
        down_fliers.sort()
    else:
        down_whisker = np.nan
        down_fliers = []
    up_whisker_mask = weak_ratio < (q3 + (1.5 * iqr))
    if up_whisker_mask.any():
        up_whisker = max(weak_ratio[up_whisker_mask])
        up_fliers = [el for el in weak_ratio if el > up_whisker]
        num_all_fliers += len(up_fliers)
        up_fliers = list(set(up_fliers))
        up_fliers.sort()
        fliers += up_fliers
    else:
        up_whisker = np.nan
        up_fliers = []
    num_unique_fliers = len(down_fliers) + len(up_fliers)
    if num_unique_fliers > 1000:
        if len(down_fliers) > 500 and len(up_fliers) > 500:
            down_fliers = down_fliers[:500]
            up_fliers = up_fliers[-500:]
        elif any([len(el) > 500 for el in [down_fliers, up_fliers]]):
            if len(down_fliers) > 500:
                down_fliers = down_fliers[:(1000 - len(up_fliers))]
            else:
                up_fliers = up_fliers[-(1000 - len(up_fliers)):]
    fliers = down_fliers + up_fliers
    if num_unique_fliers > len(fliers):
        print('number of all fliers: {}, unique: {}, limited to: {}'.format(num_all_fliers, num_unique_fliers, len(fliers)))
    else:
        print('number of all fliers: {}, unique: {}'.format(num_all_fliers, num_unique_fliers))
    result = {
        'label': 'chr {}'.format(ch),
        'whislo': down_whisker,
        'q1': q1,
        'med': median,
        'q3': q3,
        'whishi': up_whisker
    }
    print('box for chr {}:\n{}'.format(ch, result))
    result['fliers'] = fliers
    print('IQR: {}'.format(iqr))
    tmp_file = os.path.join(indir, 'boxplot_nodata_chr{}.json'.format(ch))
    with open(tmp_file, 'w') as fp:
        json.dump(result, fp)
    print('box saved to: {}'.format(tmp_file))
    return result


def plot_boxplot(indir, title, outdir=None, fromchr=1, tochr=23, new=False):
    boxplots = []
    for ch in range(fromchr, tochr+1):
        json_file = os.path.join(indir, 'boxplot_nodata_chr{}.json'.format(ch))
        if os.path.exists(json_file) and not new:
            with open(json_file, 'r') as f:
                boxplots.append(json.load(f))
            print('Boxplot of chr {} loaded from file {}'.format(ch, json_file))
        else:
            boxplots.append(find_weak(ch, indir))
    fig, ax = plt.subplots(figsize=(10, 8))
    ax.bxp(boxplots, showfliers=True, flierprops={'marker': '.', 'markersize': 1, 'markerfacecolor': 'green',
                                                  'linestyle': 'none'})
    ax.set_ylabel("Ratio of lack of data")
    ax.set_ylim(-0.1, 1.1)
    ax.set_title(title)
    if outdir is None:
        outdir = '.'
    outfile = os.path.join(outdir, "{}_{}-{}.png".format(title.replace(' ', '-'), fromchr, tochr))
    plt.savefig(outfile)
    plt.close()
    print('Plot for "{}" was saved to {}'.format(title, outfile))


parser = argparse.ArgumentParser('Create boxplot for a given dataset showing ratio of lack of data (number of -1).')
parser.add_argument('indir', type=str, help='Input directory with X matrices for all chromosomes.')
parser.add_argument('title', type=str, help='Title of the plot.')
parser.add_argument('--outdir', type=str, default=None, help='Output directory.')
parser.add_argument('--chr', type=str, default='1-23', help='Range of chromosomes to plot.')
parser.add_argument('--new', action='store_true', help='Do not load boxes from json file(s) even if it is possible.')
args = parser.parse_args()

fromchr, tochr = args.chr.split('-')
plot_boxplot(args.indir, args.title, outdir=args.outdir, fromchr=int(fromchr), tochr=int(tochr), new=args.new)
