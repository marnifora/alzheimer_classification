import numpy as np
import sys
sys.path.insert(0, '../')
import exceptions

'''
See readme.txt for input, output and possible options.
'''


def vcf_to_matrix(c, inp, outdir):

    o = open(inp, 'r')
    
    header = 0
    for line in o:
        if not line.startswith('##'):
            break
        header += 1
     
    pat = 0
    p = open('%spid_chr%s.txt' % (outdir, c), 'w')
    line = line.split()
    for el in line[9:]:
        p.write(el.strip() + '\n')
        pat += 1
    header += 1
    p.close()

    snp = len(o.readlines())
        
    o.close()

    o = open(inp, 'r')
    
    for i in range(header):
        o.readline()

    s = open('%ssnps_chr%s.txt' % (outdir, c), 'w')
    matrix = np.zeros(shape=(pat, snp, 2), dtype=np.int8)
    for i, line in enumerate(o):
    
        line = line.split()
        s.write('%s\t%s\t%s\n' % (line[1], line[3], line[4]))
        
        for j, e in enumerate(line[9:]):
            e = e.split(':')[0]
            try:
                v1, v2 = map(int, e.split("/"))
            except ValueError:
                try:
                    v1, v2 = map(int, e.split("|"))
                except ValueError:
                    v1, v2 = -1, -1
            matrix[j, i] = [v1, v2]
    s.close()   
    np.save('%smatrix_chr%s.npy' % (outdir, c), matrix)
    
    o.close()
    return "%s\t%d\t%d\n" % (c, snp, pat)


ch = '1'
outdir = './'
for q in range(len(sys.argv)):
    if sys.argv[q] == '-chr':
        ch = sys.argv[q+1]
    if sys.argv[q] == '-input':
        inp = sys.argv[q+1]
    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q+1]

if 'inp' not in globals():
    raise exceptions.NoParameterError('inp', 'name of input file')

print(vcf_to_matrix(ch, inp, outdir))
