import SNPs_subsets.subset_funcs as sfuncs
import sys
import exceptions
import corporate_funcs as funcs


def map_rows_to_locs(borutadir, directory, ch, run, outfile, subsettype, perc, snpsubset, snprun, rsnumber):

    if subsettype == 'best':
        subset = sfuncs.best_snp(borutadir, directory, ch, run, perc, snpsubset, snprun)
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
                snp = line.strip().split('\t')
                if rsnumber:
                    outfile.write('chr%d\t%d\t%d\t%s\n' % (ch, int(snp[0]) - 1, int(snp[0]), snp[3]))
                else:
                    outfile.write('chr%d\t%d\t%d\n' % (ch, int(snp[0]) - 1, int(snp[0])))
                try:
                    s = next(subset)
                except StopIteration:
                    stopped = True
                    break
        if not stopped:
            print(next(subset))
            raise exceptions.OtherError('Not all selected SNPs for chr %d were found in the SNP file!' % ch)
    return 0


def map_locs_to_rows(directory, infile, run, fixed, name, analysistype='frombed'):

    run = funcs.establish_run(analysistype, fixed, '%s%s/' % (directory, analysistype), run)
    bedfile = open(infile, 'r')
    bedline = bedfile.readline().strip().split()
    ch = int(bedline[0].strip('chr'))
    pos = int(bedline[2])
    chrlist = [ch]
    numsnps = 0
    while bedline:
        out = open('%s%s/%s_snps_chr%d_%d.txt' % (directory, analysistype, analysistype, ch, run), 'w')
        print('Rewriting SNPs for chr %s' % ch)
        with open('%smatrices/snps_chr%d.txt' % (directory, ch), 'r') as snpfile:
            for i, snpline in enumerate(snpfile):
                if snpline.startswith(str(pos)):
                    out.write('%d\n' % i)
                    numsnps += 1
                    bedline = bedfile.readline().strip().split()
                    if not bedline:
                        break
                    pos = int(bedline[2])
                    if int(bedline[0].strip('chr')) != ch:
                        ch = int(bedline[0].strip('chr'))
                        chrlist.append(ch)
                        break
    out.close()
    bedfile.close()
    funcs.runs_file_add(analysistype, '%s%s/' % (directory, analysistype), run, '%d\t%s\t%s\t%d\t%s\t' % (run, infile, name, numsnps, funcs.make_chrstr(chrlist)))
    print('%s SNPs were rewritten to a %s file!' % (numsnps, analysistype))
    return 0


def make_bedfile(name, borutadir, directory, chrlist, subsettype, run, perc, rsnumber):

    snpsubset, snprun = None, None
    if subsettype == 'best':
        perc, snpsubset, snprun = sfuncs.check_borutarun(borutadir, run, perc)
        outfile = open('%sboruta/locs_%s_bestsnps_%d_%d.bed' % (borutadir, name, perc, run), 'w')
    else:
        outfile = open('%s%s/locs_%s_%s_snps_%d.bed' % (directory, subsettype, name, subsettype, run), 'w')
    for ch in chrlist:
        map_rows_to_locs(borutadir, directory, ch, run, outfile, subsettype, perc, snpsubset, snprun, rsnumber)
    outfile.close()


chrlist = [i for i in range(1, 24)]
run = None
perc = None
rsnumber=False
borutadir = None
makebed = False
rewrite = False
fixed = False

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
    if sys.argv[q] == '-rsnumber':
        rsnumber = True
    if sys.argv[q] == '-borutadir':
        borutadir = sys.argv[q + 1]
    if sys.argv[q] == '-infile':
        infile = sys.argv[q + 1]
    if sys.argv[q] == '-make':
        makebed = True
    if sys.argv[q] == '-rewrite':
        rewrite = True
    if sys.argv[q] == '-fixed':
        fixed = True

if borutadir is None and 'directory' in globals():
    borutadir = directory

if makebed:
    make_bedfile(name, borutadir, directory, chrlist, type, run, perc, rsnumber)

if rewrite:
    map_locs_to_rows(directory, infile, run, fixed, name)
