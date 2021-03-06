import sys
from collections import OrderedDict
import subset_funcs
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs


def find_crossed(dataset, chrlist, fixed, run, borutaruns, perc):
    """
    Searching for crossed SNPs among given data sets, writing them into files.
    :param dataset: (dict) the keys are name of data sets, values are directories to folders with them
    :param chrlist: (list) chromosomes for analysis
    :param fixed: (boolean) if number of run can be overwritten
    :param run: (int or None) number of run given as a parameter - None if not given
    :return: number of crossed SNPs for the given data sets
    """

    crossed_snps = 0
    runs = {}
    for setname in dataset.keys():
        runs[setname] = funcs.establish_run('crossed', fixed, dataset[setname] + 'crossed/', run)

    for ch in chrlist:

        print('Analysis for chromosome %d has started!' % ch)
        crossed, ref = subset_funcs.first_intersection(dataset, ch, borutaruns, perc)

        '''
        for setname in list(dataset.keys())[2:]:
            set = open('%smatrices/snps_chr%d.txt' % (dataset[setname], ch), 'r')
            crossed = subset_funcs.next_intersection(set, crossed, ref, ch)
        '''

        for n, setname in enumerate(dataset.keys()):
            file = open('%scrossed/crossed_snps_chr%d_%d.txt' % (dataset[setname], ch, runs[setname]), 'w')
            for snp in sorted(crossed.keys()):
                file.write('%d\n' % crossed[snp][n])
            file.close()

        crossed_snps += len(crossed)
        print('For chr %d found %d crossed SNPs.' % (ch, len(crossed)))

    print('Run information for every dataset is writing to the file.')

    for setname in dataset.keys():
        funcs.runs_file_add('crossed', dataset[setname] + 'crossed/', runs[setname], '%d\t%s\t%s\t%s\t%d\t%s\t%d\n' %
                            (runs[setname], setname, ', '.join([k for k in dataset.keys() if k != setname]),
                             funcs.make_chrstr(chrlist), crossed_snps, ','.join(list(map(str, borutaruns.values()))),
                             perc))

    return crossed_snps


dataset = OrderedDict()
chrlist = [i for i in range(1, 24)]
fixed = False
run = None
borutaruns = None
perc = 90

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

if 'pp' in globals():
    if borutaruns is None:
        borutaruns = OrderedDict([(n, num) for n, num in zip(dataset.keys(), [pp]*len(dataset))])
    else:
        for name in [nn for nn in dataset.keys() if nn not in borutaruns.keys()]:
            borutaruns[name] = pp

found = find_crossed(dataset, chrlist, fixed, run, borutaruns, perc)
print('%d crossed SNPs found!' % found)
