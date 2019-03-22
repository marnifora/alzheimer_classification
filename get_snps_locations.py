import sys
import os
from collections import OrderedDict, deque
import exceptions
import corporate_funcs as funcs


def snps_locations(ch, directory, analysistype, analysisrun):

    towrite = ''
    numfile = open('%s%s/%s_snps_chr%d_%d.txt' % (directory, analysistype, analysistype, ch, analysisrun), 'r')
    numbers = deque(list(map(int, numfile.readlines())))
    numfile.close()
    try:
        number = numbers.popleft()
        with open('%smatrices/snps_chr%d.txt' % (directory, ch), 'r') as file:
            for i, line in enumerate(file):
                if i == number:
                    towrite += 'chr%d\t%s\t%d\n' % (ch, line.split()[0], int(line.split()[0])+1)
                    number = numbers.popleft()
    except IndexError:
        pass
    return towrite


dataset = OrderedDict()
chrlist = [i for i in range(1, 24)]
fixed = False

for q in range(len(sys.argv)):
    if sys.argv[q] == '-dataset':
        name = sys.argv[q+1]
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            directory = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')
    if sys.argv[q] == '-chr':
        chrlist = funcs.read_chrstr(sys.argv[q + 1])
    if sys.argv[q] == '-analysistype':
        analysistype = sys.argv[q+1]
    if sys.argv[q] == '-analysisrun':
        analysisrun = int(sys.argv[q+1])
    if sys.argv[q] == '-output':
        outfile = sys.argv[q+1]
    if sys.argv[q] == '-fixed':
        fixed = True

if 'outfile' in globals():
    output = open('%s%s' % (directory, outfile), 'w')
else:
    file = '%s%s_%s%d_locations.bed' % (directory, name, analysistype, analysisrun)
    if os.path.isfile(file) and not fixed:
        raise exceptions.FileOverwriteError(file)
    output = open(file, 'w')

for ch in chrlist:
    print('Rewriting for chromosome %d started!' % ch)
    output.write(snps_locations(ch, directory, analysistype, analysisrun))

