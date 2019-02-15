import sys
import exceptions

'''
Input:
- map_ascii file
- dbsnp_ascii.txt
- ped file

Output:
- snps_ref.txt
- missing_snps.txt
- genome_stats.txt
- pid_chr.txt

'''


def make_ref(dbsnp, plink, indir, outdir):

    dbsnp = open('%s%s' % (indir, '_ascii.'.join(dbsnp.split('.'))), 'r')

    mapfile = open('%s%s_ascii.map' % (indir, plink), 'r')

    ref = open('%ssnps_ref.txt' % outdir, 'w')
    missing = open('%smissing_snps_ref.txt' % outdir, 'w')

    snps = {}

    snp = mapfile.readline().split()[:2]

    for line in dbsnp:
        if not snp:
            break
        line = line.split()
        while line[0] > snp[1] and snp:
            missing.write('%s\n' % snp[1])
            snp = mapfile.readline().split()[:2]
        if line[0] == snp[1]:
            ref.write('%s\t%s\t%s\n' % (snp[0], snp[1], line[1]))
            snps[snp[0]] = snps.setdefault(snp[0], 0) + 1
            snp = mapfile.readline().split()[:2]

    dbsnp.close()
    mapfile.close()
    ref.close()
    missing.close()

    return snps


def make_pid(plink, indir, outdir):

    pedfile = open('%s%s.ped' % (indir, plink), 'r')
    o = open('%spid_chr.txt' % outdir, 'w')
    pat = 0
    for line in pedfile:
        line = line.split()
        o.write('%s\n' % line[1])
        pat += 1
    o.close()
    pedfile.close()

    return pat


def genome_stats(pat, snps, outdir):

    stats = open('%sgenome_stats.txt' % outdir, 'w')
    chs = list(map(int, snps.keys()))
    chs.sort()
    for ch in chs:
        stats.write('%d\t%d\t%d\n' % (ch, pat, snps[str(ch)]))
    stats.close()


indir = './'
outdir = './'
for q in range(len(sys.argv)):
    if sys.argv[q] == '-dbsnp':
        dbsnp = sys.argv[q+1]
    if sys.argv[q] == '-plink':
        plink = sys.argv[q+1]
    if sys.argv[q] == '-indir':
        indir = sys.argv[q+1]
    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q+1]

if 'plink' not in globals():
    raise exceptions.NoParameterError('plink', 'name of plink files')
if 'dbsnp' not in globals():
    raise exceptions.NoParameterError('dbsnp', 'name of file with list of SNPs assigned to their reference values')

pat = make_pid(plink, indir, outdir)
snps = make_ref(dbsnp, plink, indir, outdir)
genome_stats(pat, snps, outdir)
