import sys

'''
See readme.txt for input, output and possible options.
'''


def run_stats(inp, outp):

    o = open(inp, 'r')

    header = 0
    for line in o:
        if not line.startswith('##'):
            break
        header += 1

    pat = 0

    p = open(outp, 'w')
    string = ''
    line = line.split()
    for el in line[9:]:
        string += el + '\n'
        pat += 1
      
    snp = len(o.readlines())

    p.write('Number of SNPs:\t%d\nNumber of patients:\t%d\nPatients identifiers:\n%s' % (snp, pat, string))
    p.close()
    o.close()


inp = ''
for q in range(len(sys.argv)):
    if sys.argv[q] == '-input':
        inp = sys.argv[q+1]

outp = inp.split('.')[0]+'_stats.txt'
run_stats(inp, outp)
