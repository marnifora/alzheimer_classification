import sys
import numpy as np
import os
sys.path.insert(0, '../')
import exceptions


'''
Input:
- {plink}.map
- {plink}.ped
- snps_ref.txt
- genome_stats.txt

Output:
- snps_chr{}.txt
- matrix_chr{}.npy

'''


def make_snps_ref(outdir):

    snps_ref = {}
    o = open('%ssnps_ref.txt' % outdir, 'r')
    for line in o:
        line = line.split()
        if line[1] in snps_ref:
            raise exceptions.KeyOverwriting(line[1])
        snps_ref[line[1]] = line[2]
    o.close()

    return snps_ref


def make_snps_count(plink, indir):

    pedfile = open('%s%s.ped' % (indir, plink), 'r')
    snps_count = [{'A': 0, 'C': 0, 'T': 0, 'G': 0} for _ in range(len(pedfile.readline().split()[6:])//2)]
    pedfile.seek(0, 0)
    for i, line in enumerate(pedfile):
        line = line.split()[6:]
        for j, a, aa in zip(list(range(len(line)//2)), line[0::2], line[1::2]):
            if a != '0':
                snps_count[j][a] += 1
            if aa != '0':
                snps_count[j][aa] += 1
    pedfile.close()

    return snps_count


def write_snps_list(plink, indir, outdir, overwrite):

    snps_ref = make_snps_ref(outdir)
    snps_count = make_snps_count(plink, indir)

    mapfile = open('%s%s.map' % (indir, plink), 'r')
    prevch = mapfile.readline()[0]
    mapfile.seek(0, 0)
    filename = '%ssnps_chr%s.txt' % (outdir, prevch)
    if not os.path.isfile(filename) or overwrite:
        file = open(filename, 'w')
    else:
        raise exceptions.FileOverwriteError(filename)
    snps_val = {prevch: []}
    for i, line in enumerate(mapfile):
        line = line.split()
        if line[0] != prevch:
            file.close()
            filename = '%ssnps_chr%s.txt' % (outdir, line[0])
            if not os.path.isfile(filename) or overwrite:
                file = open(filename, 'w')
            else:
                raise exceptions.FileOverwriteError(filename)
            prevch = line[0]
            snps_val[prevch] = []
        snps_val[prevch].append(sorted(snps_count[i], key=snps_count[i].get, reverse=True)[:len([jj for jj in snps_count[i].values() if jj != 0])])
        if not snps_val[prevch][-1]:
            snps_val[prevch][-1].append(snps_ref[line[1]])
        elif snps_val[prevch][-1][0] != snps_ref[line[1]]:
            # print('SNP %s is untypical: reference = %s, its values = %s' % (line[1], snps_ref[line[1]], snps_count[i]))
            try:
                snps_val[prevch][-1].remove(snps_ref[line[1]])
            except ValueError:
                pass
            snps_val[prevch][-1].insert(0, snps_ref[line[1]])
        file.write('%s\t%s\t%s\t%s\n' % (line[3], snps_ref[line[1]], ', '.join(snps_val[prevch][-1][1:]), line[1]))
    mapfile.close()
    file.close()

    return snps_val


def write_matrix(plink, indir, outdir, snps_val):

    stats = open('%sgenome_stats.txt' % outdir, 'r')
    pedfile = open('%s%s.ped' % (indir, plink), 'r')

    done = 0
    for line in stats:
        ch, snp, pat = line.split()
        pat, snp = int(pat), int(snp)
        matrix = np.zeros(shape=(pat, snp, 2), dtype=np.int8)
        pedfile.seek(0, 0)
        for p, row in enumerate(pedfile):
            row = row.split()[6+done:6+done+(2*snp)]
            s = 0
            n = nn = False
            v1 = v2 = -2
            for a, aa in zip(row[0::2], row[1::2]):
                if a == '0':
                    v1 = -1
                    n = True
                if aa == '0':
                    v2 = -1
                    nn = True
                if not n or not nn:
                    for i, base in enumerate(snps_val[ch][s]):
                        if a == base:
                            v1 = i
                        if aa == base:
                            v2 = i
                if v1 == -2 or v2 == -2:
                    raise exceptions.PlinkWrongValue('%d+%d' % (done, s), ch, p, a, aa)
                matrix[p][s] = [v1, v2]
                s += 1
                n = nn = False
                v1 = v2 = -2
        np.save('%smatrix_chr%s.npy' % (outdir, ch), matrix)
        done += 2*snp

    stats.close()
    pedfile.close()


indir = './'
overwrite = False
for q in range(len(sys.argv)):
    if sys.argv[q] == '-plink':
        plink = sys.argv[q+1]
    if sys.argv[q] == '-indir':
        indir = sys.argv[q+1]
    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q+1]
    if sys.argv[q] == '-overwrite':
        overwrite = True

if 'outdir' not in globals():
    outdir = indir

if 'plink' not in globals():
    raise exceptions.NoParameterError('plink', 'name of plink files')

snps_val = write_snps_list(plink, indir, outdir, overwrite)
write_matrix(plink, indir, outdir, snps_val)
