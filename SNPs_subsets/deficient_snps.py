import sys
import numpy as np
from collections import OrderedDict, deque
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs


def find_weak(ch, directory, perc, borutarun, thresh):

    weak = []
    x = np.load('%sboruta/X_train_chr%d_%d_%d.npy' % (directory, ch, perc, borutarun))
    for i, col in enumerate(x.T):
        if len([a for a in col if a == -1]) > thresh*len(col):
            weak.append(i)
    return locate_best(ch, directory, weak, perc, borutarun)


def locate_best(ch, directory, snps, perc, borutarun):

    locs = []
    snps = deque(snps)
    snp = snps.popleft()
    with open('%sboruta/bestsnps_chr%d_%d_%d.txt' % (directory, ch, perc, borutarun), 'r') as file:
        all = int(file.readline().strip())
        file.readline()
        for i, line in enumerate(file):
            if i == snp:
                locs.append(int(line.strip()))
                try:
                    snp = snps.popleft()
                except IndexError:
                    return locs, all
    if snps:
        raise exceptions.OtherError('No all SNPs were found - there is no enough SNPs in bestsnps file.')


dataset = OrderedDict()
chrlist = [i for i in range(1, 24)]
fixed = False
run = None
borutaruns = None
perc = 90
thresh = 0.1

for q in range(len(sys.argv)):
    if sys.argv[q] == '-dataset':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            dataset[sys.argv[q + 1]] = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')
    if sys.argv[q] == '-chr':
        chrlist = funcs.read_chrstr(sys.argv[q + 1])
    if sys.argv[q] == '-run':
        run = int(sys.argv[q+1])
    if sys.argv[q] == '-fixed':
        fixed = True
    if sys.argv[q] == '-perc':
        perc = int(sys.argv[q + 1])
    if sys.argv[q] == '-borutarun':
        if sys.argv[q + 1] in dataset.keys():
            if borutaruns is None:
                borutaruns = OrderedDict()
            try:
                borutaruns[sys.argv[q+1]] = int(sys.argv[q+2])
            except ValueError:
                raise exceptions.NoParameterError('borutarun',
                                                  'After name of data set should appear number of boruta run')
        else:
            try:
                pp = int(sys.argv[q+1])
            except ValueError:
                raise exceptions.NoParameterError('borutarun',
                                                  'After -borutarun should appear name of data set and its run number ' +
                                                  'or one run number which is the same for every data set.')
    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q+1]

    if sys.argv[q] == '-thresh':
        thresh = float(sys.argv[q+1])


if 'pp' in globals():
    if borutaruns is None:
        borutaruns = OrderedDict([(n, num) for n, num in zip(dataset.keys(), [pp]*len(dataset))])
    else:
        for name in [nn for nn in dataset.keys() if nn not in borutaruns.keys()]:
            borutaruns[name] = pp

if 'outdir' in globals():
    run = funcs.establish_run('deficient', fixed, outdir, run)
    output = outdir

l = 0
all = 0
for name, directory in dataset.items():
    print('Analysis for %s dataset' % name)
    if 'outdir' not in globals():
        output = '%sdeficient/' % directory
        run = funcs.establish_run('deficient', fixed, output, run)
    for ch in chrlist:
        print('Checking SNPs for chromosome %d has just started!' % ch)
        locs, al = find_weak(ch, directory, perc, borutaruns[name], thresh)
        l += len(locs)
        all += al
        file = open('%sdeficient_snps_chr%d_%d.txt' % (output, ch, run), 'w')
        file.write('\n'.join(list(map(str, locs))))
        file.close()

    if 'outdir' not in globals():
        funcs.runs_file_add('deficient', output, run, '%d\t%s\t%s\t%.2f\t%d\t%s\n' %
                            (run, name, borutaruns[name], thresh, l, '%.1f' % (l/all*100) + '%'))

if 'outdir' in globals():
    funcs.runs_file_add('deficient', output, run, '%d\t%s\t%s\t%.2f\t%d\t%s\n' %
                        (run, '+'.join(dataset.keys()), ','.join(list(map(str, borutaruns.values()))), thresh, l,
                         '%.1f' % (l / all * 100) + '%'))

print('Number of deficient SNPs: %d (%s)' % (l, '%.1f' % (l/all*100) + '%'))
