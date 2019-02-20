import sys
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
    sets = list(dataset.keys())
    runs = {}
    for setname in sets:
        runs[setname] = funcs.establish_run('shared', fixed, dataset[setname], run)

    for ch in chrlist:

        print('Analysis for chromosome %d has started!' % ch)
        shared = {}
        ref = {}
        set0 = open('%ssnps_chr%d.txt' % (dataset[sets[0]], ch), 'r')
        set1 = open('%ssnps_chr%d.txt' % (dataset[sets[1]], ch), 'r')
        line = [set0.readline().split(), set1.readline().split()]
        order = [0, 0]

        while line[0] and line[1]:

            if int(line[0][0]) < int(line[1][0]):

                line[0] = set0.readline().split()
                order[0] += 1

            elif int(line[0][0]) == int(line[1][0]):

                if line[0][1] == line[1][1]:
                    shared[int(line[0][0])] = order
                    ref[int(line[0][0])] = line[0][1]
                    line = [set0.readline().split(), set1.readline().split()]
                    order = [x+1 for x in order]
                else:
                    raise exceptions.SNPReferenceError(ch, int(line[0][0]), line[0][1], line[1][1])
            else:
                line[1] = set1.readline().split()
                order[1] += 1

        print('Found shared SNPs from chr %d for two first sets.' % ch)

        for setname in sets[2:]:
            set = open('%ssnps_chr%d.txt' % (dataset[setname], ch), 'r')
            snps = iter(shared.keys())
            try:
                snp = next(snps)
            except StopIteration:
                break
            for order, line in enumerate(set):
                line = line.split()
                if int(line[0]) == snp:
                    if line[1] == ref[snp]:
                        shared[snp].append(order)
                    else:
                        raise exceptions.SNPReferenceError(ch, snp, ref[snp], line[1])
                elif int(line[0]) > snp:
                    del(shared[snp])
                    try:
                        snp = next(snps)
                    except StopIteration:
                        break

        print('Writing found shared SNPs from chr %d to the file.' % ch)

        for n, setname in enumerate(sets):
            file = open('%sshared_snps_chr%d_%d.txt' % (dataset[setname], ch, runs[setname]), 'w')
            for snp in sorted(shared.keys()):
                file.write('%d\n' % shared[snp][n])
            file.close()

        shared_snps += len(shared)

    print('Run information for every dataset is writing to the file.')

    for setname in sets:
        run_file = open('%sshared_runs.txt' % dataset[setname], 'a')
        run_file.write('%d\t%s\t%s\t%s\t%d\n' % (runs[setname], setname, ', '.join([k for k in dataset.keys() if k != setname]),
                                                 funcs.make_chrstr(chrlist), shared_snps))
        run_file.close()

    return shared_snps


dataset = {}
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
