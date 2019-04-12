import subset_funcs as sfuncs
import sys
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs


def map_rows_to_locs(directory, ch, run, outfile, subsettype, perc, snpsubset, snprun):

    if subsettype == 'best':
        subset = sfuncs.best_snp(directory, ch, run, perc, snpsubset, snprun)
    elif subsettype == 'shared':
        subset = sfuncs.shared_snp(directory, ch, run)
    elif subsettype == 'crossed':
        subset = sfuncs.crossed_snp(directory, ch, run)
    s = next(subset)
    print('Mapping rows to locations for chromosome %d' % ch)
    with open('%smatrices/snps_chr%d.txt' % (directory, ch), 'r') as snpsfile:
        stopped = False
        for i, line in enumerate(snpsfile):
            if i == s:
                snp = line.split()
                outfile.write('chr%d\t%d\t%d\n' % (ch, int(snp[0]), int(snp[0]) + 1))
                try:
                    s = next(subset)
                except StopIteration:
                    stopped = True
                    break
        if not stopped:
            print(next(subset))
            raise exceptions.OtherError('Not all selected SNPs for chr %d were found in the SNP file!' % ch)
    return 0


def make_bedfile(name, directory, chrlist, subsettype, run, perc):

    snpsubset, snprun = None, None
    if subsettype == 'best':
        perc, snpsubset, snprun = sfuncs.check_borutarun(directory, run, perc)
        outfile = open('%sboruta/locs_%s_bestsnps_%d_%d.bed' % (directory, name, perc, run), 'w')
    else:
        outfile = open('%s%s/locs_%s_%s_snps_%d.bed' % (directory, subsettype, name, subsettype, run), 'w')
    for ch in chrlist:
        map_rows_to_locs(directory, ch, run, outfile, subsettype, perc, snpsubset, snprun)
    outfile.close()


chrlist = [i for i in range(1, 24)]
run = None
perc = None

for q in range(len(sys.argv)):
    if sys.argv[q] == '-dataset':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            name = sys.argv[q + 1]
            directory = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')
    if sys.argv[q] == '-chr':
        chrlist = funcs.read_chrstr(sys.argv[q + 1])
    if sys.argv[q] == '-run':
        run = int(sys.argv[q+1])
    if sys.argv[q] == '-perc':
        perc = int(sys.argv[q + 1])
    if sys.argv[q] == '-type':
        type = sys.argv[q + 1]

make_bedfile(name, directory, chrlist, type, run, perc)
