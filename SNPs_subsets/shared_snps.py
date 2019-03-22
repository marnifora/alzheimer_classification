import sys
import subset_funcs
from collections import OrderedDict
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs

'''
Input:
For each dataset:
- snps_chr{}.txt - list of SNPs from given dataset
Output:
For each dataset:
- shared_runs.txt - list of runs of shared SNPs analysis
- shared_snps_chr{}_{run}.txt - list of shared SNPs from chr {} from run {}
'''


def find_shared(dataset, chrlist, fixed, run):
    """
    Searching for shared SNPs among given data sets, writing them into files.
    :param dataset: (dict) the keys are name of data sets, values are directories to folders with them
    :param chrlist: (list) chromosomes for analysis
    :param fixed: (boolean) if number of run can be overwritten
    :param run: (int or None) number of run given as a parameter - None if not given
    :return: number of shared SNPs for the given data sets
    """

    shared_snps = 0
    runs = {}
    for setname in dataset.keys():
        runs[setname] = funcs.establish_run('shared', fixed, dataset[setname] + 'shared/', run)

    for ch in chrlist:

        print('Analysis for chromosome %d has started!' % ch)
        shared, ref = subset_funcs.first_intersection(dataset, ch)

        '''
        for setname in list(dataset.keys())[2:]:
            set = open('%smatrices/snps_chr%d.txt' % (dataset[setname], ch), 'r')
            shared = subset_funcs.next_intersection(set, shared, ref, ch)
        '''

        print('Writing found shared SNPs from chr %d to the file.' % ch)

        for n, setname in enumerate(dataset.keys()):
            file = open('%sshared/shared_snps_chr%d_%d.txt' % (dataset[setname], ch, runs[setname]), 'w')
            for snp in sorted(shared.keys()):
                file.write('%d\n' % shared[snp][n])
            file.close()

        shared_snps += len(shared)

    print('Run information for every dataset is writing to the file.')

    for setname in dataset.keys():
        funcs.runs_file_add('shared', dataset[setname] + 'shared/', runs[setname], '%d\t%s\t%s\t%s\t%d\n' %
                            (runs[setname], setname, ', '.join([k for k in dataset.keys() if k != setname]),
                             funcs.make_chrstr(chrlist), shared_snps))

    return shared_snps


dataset = OrderedDict()
chrlist = [i for i in range(1, 24)]
fixed = False
run = None

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

found = find_shared(dataset, chrlist, fixed, run)
print('%d shared SNPs found!' % found)
