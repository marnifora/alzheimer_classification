import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
import sys
import numpy as np
from collections import deque
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs


def lower_threshold(thresh, linkage, pat):
    toremove = set()
    for row in linkage:
        if row[2] < thresh:
            if row[0] >= pat:
                r = int(row[0]) - pat
                toremove.add(int(linkage[r][0]))
                toremove.add(int(linkage[r][1]))
            else:
                toremove.add(int(row[0]))
        else:
            break
    return toremove


def upper_threshold(thresh, linkage):
    clusters = hc.fcluster(linkage, thresh, criterion='distance')
    best = np.argmax(np.bincount(clusters))
    selected = [i for i, el in enumerate(clusters) if el == best]
    return selected


def diagnoses_dist(dir, patients):
    patients.sort()
    patients = deque(patients)
    pat = patients.popleft()
    diagnoses = {'0': 0, '1': 0}
    with open('%smatrices/Y_chr.csv' % dir, 'r') as file:
        for i, line in enumerate(file):
            if i == pat:
                diagnoses[line.strip().split(',')[-1]] += 1
                try:
                    pat = patients.popleft()
                except IndexError:
                    break
    print('Healthy: %d' % diagnoses['0'])
    print('Ill: %d' % diagnoses['1'])
    return 0


run = None
fixed = False

for q in range(len(sys.argv)):

    if sys.argv[q] == '-dataset':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            name = sys.argv[q + 1]
            dir = sys.argv[q + 2]
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

    if sys.argv[q] == '-run':
        run = int(sys.argv[q+1])

    if sys.argv[q] == '-fixed':
        fixed = True

if 'name' not in globals():
    raise exceptions.NoParameterError('dataset', 'Name and directory of the data set must be given!')

if 'indir' not in globals():
    indir = '%ssimilar/' % dir

if 'outdir' not in globals():
    outdir = '%ssimilar/' % dir

sims = np.load('%s%s_similarities.npy' % (indir, name))
linkage = hc.linkage(sp.distance.squareform(sims), method='average')
pat = sims.shape[0]

if 'lower' in globals():
    toremove = lower_threshold(lower, linkage, pat)
    print('Number of patients to remove (below the lower threshold %.4f): %d' % (lower, len(toremove)))
else:
    print('No lower threshold given')
    toremove = set()

if 'upper' in globals():
    selected = upper_threshold(upper, linkage)
    print('Number of patients in the biggest cluster (the upper threshold %.4f): %d' % (upper, len(selected)))
else:
    print('No upper threshold given')
    selected = [i for i in range(sims.shape[0])]

final = [el for el in selected if el not in toremove]

print('Number of selected patients: %d' % len(final))
diagnoses_dist(dir, final)

run = funcs.establish_run('similar', fixed, outdir, run)
file = open('%ssimilar_patients_%d.txt' % (outdir, run), 'w')
file.write('\n'.join([str(p) for p in final]))
file.close()

run_file = open('%ssimilar_runs.txt' % outdir, 'a')
run_file.write('%d\t%s\t%.4f\t%.4f\t%d\t%d\n' % (run, name, lower, upper, len(selected), pat))
