import numpy as np
import multiprocessing as multi
import sys
from collections import OrderedDict
import os
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs


def one_process(ec, er, i, matrix, outdir, pat, snp, sc, sr):

    if er == ec:
        er -= 1
        ec = pat
    else:
        ec -= 1
    print('process %d has just started' % i)
    print('start: <%d, %d>' % (sr, sc))
    print('end: <%d, %d>' % (er, ec))
    f = open('%ssimilarities_%d.txt' % (outdir, i), 'w')
    f.write('<%d, %d>\n' % (sr, sc))
    rr = sc
    for r in range(sr, er+1):
        if r == er:
            cc = ec
        else:
            cc = pat
        for c in range(rr, cc+1):
            f.write('%.6f\n' % count_similarity(matrix[r-1], matrix[c-1], snp))
        rr = r + 1

    f.close()
    print('process %d has just finished' % i)


def count_similarity(row1, row2, snp):
    return sum(np.clip(abs(row1 - row2), 0, 1)) / snp


def make_matrix(dataset, outdir, pat, procs):

    sims = np.zeros(shape=(pat, pat), dtype=np.float16)
    oldr = 1
    oldc = 1
    for i in range(procs):
        with open('%ssimilarities_%d.txt' % (outdir, i), 'r') as file:
            r, c = list(map(int, file.readline().strip('<>\n').split(', ')))
            if r != oldr or c != oldc:
                raise exceptions.WrongSubscripts(i - 1, oldr, oldc, i, r, c)
            for line in file:
                value = float(line.strip())
                sims[r - 1, c - 1] = value
                sims[c - 1, r - 1] = value
                c += 1
                if c > pat:
                    r += 1
                    c = r
            oldr, oldc = r, c
    np.save('%s%s_similarities.npy' % (outdir, '-'.join(list(dataset.keys()))), sims)


dataset = OrderedDict()
procs = 1
chrlist = [i for i in range(1, 24)]
snpsubset = None
snpsubset = None
snpruns = None

for q in range(len(sys.argv)):

    if sys.argv[q] == '-dataset':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            dataset[sys.argv[q + 1]] = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')
        continue

    if sys.argv[q] == '-outdir':
        if sys.argv[q + 1][0] in ['.', '~', '/']:
            outdir = sys.argv[q + 1]
        else:
            raise exceptions.NoParameterError('outdir',
                                              'After -outdir should appear a directory to output folder.')
        continue

    if sys.argv[q] == '-procs':
        procs = int(sys.argv[q+1])
        continue

    if sys.argv[q] == '-chr':
        chrlist = funcs.read_chrstr(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-snpsubset':
        snpsubset = sys.argv[q + 1]
        continue

    if sys.argv[q] == '-snprun':
        if sys.argv[q + 1] in dataset.keys():
            if snpruns is None:
                snpruns = OrderedDict()
            try:
                snpruns[sys.argv[q+1]] = int(sys.argv[q+2])
            except ValueError:
                raise exceptions.NoParameterError('snprun',
                                                  'After name of data set should appear number of SNP subset run')
        else:
            try:
                ss = int(sys.argv[q+1])
            except ValueError:
                raise exceptions.NoParameterError('snprun',
                                                  'After -snprun should appear name of data set and its run number ' +
                                                  'or one run number which is the same for every data set.')
        continue

    if sys.argv[q].startswith('-'):
        raise exceptions.WrongParameterName(sys.argv[q])

if not dataset:
    raise exceptions.NoParameterError('dataset', 'data set(s) for which similarities should be counted')

if 'outdir' not in globals():
    if len(dataset) == 1:
        outdir = '%ssimilar/' % next(iter(dataset.values()))
    else:
        raise exceptions.NoParameterError('outdir', 'to what directory output matrix should be written')

if not os.path.isdir(outdir):
    os.mkdir(outdir)

if 'ss' in globals():
    if snpruns is None:
        snpruns = OrderedDict([(n, num) for n, num in zip(dataset.keys(), [ss]*len(dataset))])
    else:
        for name in [nn for nn in dataset.keys() if nn not in snpruns.keys()]:
            snpruns[name] = ss

trainpat = set([i for i in range(sum(funcs.patients(dataset).values()))])
testpat = set()
matrix, *_ = funcs.load_data(chrlist[0], dataset, snpsubset, snpruns, testpat, trainpat)
for ch in chrlist[1:]:
    x, *_ = funcs.load_data(ch, dataset, snpsubset, snpruns, testpat, trainpat)
    matrix = np.concatenate((matrix, x), axis=1)

print('matrix loaded successfully')

pat, snp = matrix.shape
sim = (1+pat)/2 * pat
operations = sim / procs
if operations % 1 != 0.0:
    raise exceptions.WrongValueError('procs', procs, 'It has to divide number of similarities - %d - without rest.' % sim)
else:
    operations = int(operations)
# operations = 28861  # similarities to count = 1 154 440, number of processes = 40, 1 154 440 / 40 = 28 861
pp = []
srow = 1
scol = 1
erow = 1

for i in range(procs):

    suma = pat - scol + 1
    erow += 1
    while suma < operations:
        suma += pat - erow + 1
        erow += 1
    if suma != operations:
        ecol = pat - (suma - operations) + 1
        erow -= 1
    else:
        ecol = erow

    p = multi.Process(target=one_process, args=(ecol, erow, i, matrix, outdir, pat, snp, scol, srow))
    pp.append(p)
    p.start()

    scol = ecol
    srow = erow

for p in pp:
    p.join()

print('Done')
print('Number of similarities to count: %d' % sim)
print('Number of operations in one process: %d' % operations)

make_matrix(dataset, outdir, pat, procs)

print('Similarity structure was written to npy file.')
