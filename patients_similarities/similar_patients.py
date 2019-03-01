import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
import sys
import numpy as np
sys.path.insert(0, '../')
import exceptions


def lower_threshold(thresh, sim):
    linkage = hc.linkage(sp.distance.squareform(sim), method='average')


def upper_threshold(thresh, sim):
    pass


sims = {}

for q in range(len(sys.argv)):

    if sys.argv[q] == '-sims':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            sims[sys.argv[q + 1]] = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')

    if sys.argv[q] == '-outdir':
        if sys.argv[q + 1][0] in ['.', '~', '/']:
            outdir = sys.argv[q + 1]
        else:
            raise exceptions.NoParameterError('outdir',
                                              'After -outdir should appear a directory to output folder.')

    if sys.argv[q] == '-lower':
        lower = float(sys.argv[q+1])

    if sys.argv[q] == '-upper':
        upper = float(sys.argv[q+1])

if 'outdir' not in globals():
    outdir = iter(next(sims.values()))

for name, directory in sims.items():
    sim = np.load('%s%s_similarities.npy' % (directory, name))

if 'lower' in globals():

    toremove = lower_threshold(lower, sim)

