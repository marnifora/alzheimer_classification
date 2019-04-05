import numpy as np
import pandas as pd
import csv
import sys

'''
See readme.txt for input, output and possible options.
'''


def makeX(ch, indir, outdir):

    matrix = np.load('%smatrix_chr%s.npy' % (indir, ch))
    matrix = matrix[:, :, 1]
    matrix = pd.DataFrame(matrix)
    matrix.to_csv('%sX_chr%s.csv' % (outdir, ch))


def makeX_nodif(ch, indir, outdir):

    o = open('%sX_chr%s.csv' % (indir, ch), 'r')
    reader = csv.reader(o, delimiter=',')

    w = open('%sX_chr%s_nodif.csv' % (outdir, ch), 'w')
    writer = csv.writer(w, delimiter=',')

    r = open('%sdif_chr.txt' % indir, 'r')
    dif = int(r.readline().strip().split()[2])

    writer.writerow(next(reader))  # header

    j = 0
    for i, line in enumerate(reader):
        if i == dif:
            try:
                dif = int(r.readline().strip().split('\t')[2])
            except IndexError:
                dif = 0
        else:
            line[0] = j
            writer.writerow(line)
            j += 1

    o.close()
    r.close()
    w.close()
    return 'Chromosome %s, number of lines written into X_chr%s_nodif.csv: %d' % (ch, ch, j)


ch = '1'
indir = './'
for q in range(len(sys.argv)):
    if sys.argv[q] == '-chr':
        ch = sys.argv[q+1]
    if sys.argv[q] == '-indir':
        indir = sys.argv[q+1]
    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q+1]

if 'outdir' not in globals():
    outdir = indir

# makeX(ch, indir, outdir)
print(makeX_nodif(ch, indir, outdir))
