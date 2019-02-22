import sys
import seaborn as sns
import numpy as np
import matplotlib.pyplot as plt
import scipy.cluster.hierarchy as hc
import scipy.spatial as sp
sys.path.insert(0, '../')
import exceptions


def make_lists(set1, set2, dataset):

    rows1 = []
    rows2 = []
    for s, rows in [[set1, rows1], [set2, rows2]]:  # ss = {set1, set2}
        for el in s:  # s = {['adni', 'healthy'], ['rosmap', 'healthy], ..}
            if len(el) == 1:  # s = {['healthy'], ['case'], ..}
                prevpat = 0
                for d in dataset:
                    l = check_group(d, el[0])
                    rows += [ll+prevpat for ll in l]
                    prevpat += d[2]
            else:
                prevpat = 0
                for d in dataset:
                    if d[0] == el[0]:
                        l = check_group(d, el[1])
                        rows += [ll+prevpat for ll in l]
                        break
                    prevpat += d[2]
    return rows1, rows2


def check_group(dataset, group):

    datadir = dataset[1]
    if group in control_str:
        l = give_rows(datadir, 0)
    elif group in case_str:
        l = give_rows(datadir, 1)
    elif group in all_str:
        l = [i for i in range(dataset[2])]
    else:
        raise exceptions.OtherError('There is no such group of patients as %s' % group)
    return l


def give_rows(datadir, value):

    l = []
    for line in open('%sY_chr.csv' % datadir, 'r'):
        line = line.strip().split(',')
        v = int(line[1])
        if v == value:
            l.append(int(line[0]))
    return l


def get_title(dataset, set1, set2):

    title = []
    for s in [set1, set2]:
        t = [[], []]
        for el in s:
            if len(el) == 1:
                t[0] += [d[0] for d in dataset]
                subset = el[0]
            else:
                t[0].append(el[0])
                subset = el[1]
            if subset in control_str:
                t[1].append('control group')
            elif subset in case_str:
                t[1].append('case group')
            elif subset in all_str:
                t[1].append('all patients')
        title.append(': '.join(list(map(' + '.join, t))))
    return ' - '.join(title)


def add_to_set(args, sets, q):

    try:
        num = int(args[q][-1])
    except ValueError:
        num = 1
    s = []
    for el in args[q + 1:]:
        if not el.startswith('-'):
            s.append(el)
        else:
            break
    if len(s) > 2:
        raise exceptions.WrongValueError(args[q], s,
                                         'First should be name of data set, second name of subset of patients')
    if len(sets) < num:
        sets += [[] for i in range(num-len(sets))]
    sets[num - 1].append(s)
    return sets


seta = []
setb = []
dataset = []
control_str = ['healthy', 'control', 'normal', 'NL', 'CN', '0']
case_str = ['ill', 'case', 'AD', 'alzheimer', '1']
all_str = ['all', 'whole', 'every']
kdeplot = False
clustermap = False
dendrogram = False

for q in range(len(sys.argv)):

    if sys.argv[q].startswith('-dataset'):
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            dataset.append([sys.argv[q+1], sys.argv[q+2]])
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')

    if sys.argv[q] == '-matrix':
        if sys.argv[q+1][0] in ['.', '~', '/']:
            sims = np.load(sys.argv[q+1])
        else:
            raise exceptions.NoParameterError('directory', 'After -matrix should appear a directory to similarity' +
                                              'matrix written in .npy file.')

    if sys.argv[q].startswith('-seta') or sys.argv[q].startswith('-setA'):

        seta = add_to_set(sys.argv, seta, q)

    if sys.argv[q].startswith('-setb') or sys.argv[q].startswith('-setB'):

        setb = add_to_set(sys.argv, setb, q)

    if sys.argv[q] == '-outdir':
        if sys.argv[q + 1][0] in ['.', '~', '/']:
            outdir = sys.argv[q + 1]
        else:
            raise exceptions.NoParameterError('outdir',
                                              'After -outdir should appear a directory to output folder')

    if sys.argv[q] == '-kdeplot':

        kdeplot = True

    if sys.argv[q] == '-clustermap':

        clustermap = True

    if sys.argv[q] == '-dendrogram':

        dendrogram = True

if not dataset:
    raise exceptions.NoParameterError('dataset', 'There must be at least one data set given')

if 'sims' not in globals():
    raise exceptions.NoParameterError('matrix', 'There must be matrix with similarities given')

if 'outdir' not in globals():
    outdir = dataset[0][1]

for sets, name in [[seta, 'setA'], [setb, 'setB']]:
    for i, s in enumerate(sets):
        if not s:
            raise exceptions.NoParameterError(name, '')
        for el in s:
            if len(el) == 2 and el[0] not in [d[0] for d in dataset]:
                    raise exceptions.WrongValueError('%s%d' % (name, i), el, 'There is no such data set as %s' % el[0])

for d in dataset:
    o = open('%sgenome_stats.txt' % d[1], 'r')
    line = o.readline()
    pat = int(line.split()[3])
    d.append(pat)

for set1, set2 in zip(seta, setb):
    rows1, rows2 = make_lists(set1, set2, dataset)
    linkage = hc.linkage(sp.distance.squareform(sims[rows1][:, rows2]), method='average')

    if clustermap:
        cmap = sns.cubehelix_palette(20, start=1, rot=-0.7, dark=0.3, light=1, reverse=True)
        sns.clustermap(sims[rows1][:, rows2], row_linkage=linkage, col_linkage=linkage, cmap=cmap)
    if dendrogram:
        hc.dendrogram(linkage)
    if kdeplot:
        similarities = sims[rows1][:, rows2].ravel()
        sns.kdeplot(similarities, label=get_title(dataset, set1, set2))

font = {'size': 16}
plt.rc('font', **font)
plt.legend()
plt.show()
