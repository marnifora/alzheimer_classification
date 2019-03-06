import sys

'''
See readme.txt for input, output and possible options.
'''


def makeY(indir, outdir):

    pid = open('%spid_chr.txt' % indir, 'r')
    diag = open('%sdiagnoses.txt' % indir, 'r')

    y = open('%sY_chr.csv' % outdir, 'w')
    o = open('%sdif_chr.txt' % outdir, 'w')

    i = 0
    ad = 0
    nl = 0

    for pline, dline in zip(pid, diag):

        if dline.strip() == 'NL':
            y.write('%d,0\n' % (nl+ad))
            nl += 1
        elif dline.strip() == 'AD':
            y.write('%d,1\n' % (nl+ad))
            ad += 1
        else:
            o.write('%s\t%s\t%d\n' % (pline.strip(), dline.strip(), i))
        i += 1

    pid.close()
    diag.close()
    y.close()
    o.close()

    return ad, nl


def update_genome_stats(indir, pat):
    """
    Writing in addition to genome_stats.txt number of patients with diagnosis NL or AD (who can be used for further
    analysis)
    :param indir: input directory
    :param pat: number of patients with NL or AD diagnosis
    :return: None
    """

    stats = open('%sgenome_stats.txt' % indir, 'r')
    lines = stats.readlines()
    stats.close()
    stats = open('%sgenome_stats.txt' % indir, 'w')
    for line in lines:
        stats.write(line.strip() + '\t%d\n' % pat)
    stats.close()
    return None


dir = './'
for q in range(len(sys.argv)):
    if sys.argv[q] == '-dir':
        dir = sys.argv[q+1]
    if sys.argv[q] == '-indir':
        indir = sys.argv[q+1]
    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q+1]

if 'indir' not in globals():
    indir = '%smatrices/' % dir

if 'outdir' not in globals():
    outdir = '%smatrices/' % dir

ad, nl = makeY(indir, outdir)
update_genome_stats(indir, ad+nl)
print('Number of patients for further analysis: %d (contains AD=%d, NL=%d)' % (ad+nl, ad, nl))
