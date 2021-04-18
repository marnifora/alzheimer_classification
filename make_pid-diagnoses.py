import random
from collections import deque
import sys
import csv
import exceptions
import os

'''
See readme.txt for input, output and possible options.
For every database there should be different function of mapping diagnoses!!
'''


def rosmap_mapping(files, indir):
    """
    Special for data from Rosmap - there are three sources of diagnoses, with different coding.
    :param files: (list) names of files with diagnoses
    :return: (dict) with diagnoses assigned to each patient
    """

    dd = {}

    for file in files:

        o = open(indir + file, 'r')
        for line in o:
            line = line.strip().split()
            try:
                d = line[1]
            except IndexError:
                break
            if d in ['0', 'Control', '1.0']:
                dd[line[0]] = 'NL'
            elif d in ['AD', '4.0']:
                dd[line[0]] = 'AD'
            else:
                dd[line[0]] = 'DIF'
        o.close()

    return dd


def adni_mapping(files, indir):
    """
    Special for data from ADNI - there is different coding of diagnosis for every edition of ADNI.
    The differences were written in dict 'coding'.
    All information (patients' IDs, editions of ADNI, diagnoses) are written in one file (list 'files').
    :param files: (list) name(s) of file(s) with diagnoses
    :return: (dict) with diagnoses assigned to each patient
    """

    dd = {}
    date = {}

    coding = {'ADNI1': [{'1': 'NL', '2': 'DIF', '3': 'AD'}, 5],
              'ADNI2': [{'179': 'NL', '248': 'DIF', '356': 'AD'}, 4],
              'ADNIGO': [{'179': 'NL', '248': 'DIF', '356': 'AD'}, 4],
              'ADNI3': [{'1': 'NL', '2': 'DIF', '3': 'AD'}, 6]}

    reader = csv.reader(open(indir + files[0], 'r'), delimiter=',')
    next(reader)  # header
    for row in reader:
        pid = '0'*(4-len(row[2])) + row[2]
        try:
            vdate = int(''.join(row[3].split('-')))
        except ValueError:
            vdate = 0
        if pid not in date:
            date[pid] = vdate
        if date[pid] <= vdate:
            for k in coding.keys():
                if row[1] == k:
                    n = coding[k][-1]
                    for kk in coding[k][0].keys():
                        if row[n] in kk:
                            dd[pid] = coding[k][0][kk]
                            break
                    break
    return dd


def test_mapping(files, indir):
    assert len(files) == 1, files
    dd = {}
    with open(os.path.join(indir, files[0]), 'r') as f:
        for line in f:
            l = line.strip().split('\t')
            dd[l[0]] = l[1]
    return dd


def check_pidfiles(indir):
    """
    Checking if number of patients and their order in vcf files for every chromosomes is the same.
    """
    pid_files = [os.path.join(indir, el) for el in os.listdir(indir) if el.startswith('pid_chr') and el.endswith('.txt') and el != 'pid_chr.txt']
    if pid_files and os.path.exists(pid_files[0]):
        stat = open(pid_files[0], 'r')
    else:
        try:
            open("%spid_chr.txt" % outdir, 'r')
        except FileNotFoundError:
            raise exceptions.NoFileError('pid_chr1.txt')
        return 'pid_chr.txt file already exist!'

    p = 0
    for i, line in enumerate(stat):
        if line == '':
            p -= 1

    p += i+1

    stat.seek(0, 0)

    sample = random.sample(range(0, p), 100)
    sample.sort()
    ss = deque(sample)

    order = {k: 0 for k in sample}
    s = ss.popleft()

    for i, line in enumerate(stat):

        if i == s:
            order[s] = line.strip()
            try:
                s = ss.popleft()
            except IndexError:
                break

    stat.close()

    for file in pid_files[1:]:
        c = file.split('pid_chr')[-1].replace('.txt', '')
        stat = open(file, 'r')

        ss = deque(sample)
        s = ss.popleft()

        for i, line in enumerate(stat):

            if i == s:
                if line.strip() != order[s]:
                    raise exceptions.FilesError('id', i, c)
                try:
                    s = ss.popleft()
                except IndexError:
                    pass

        if i+1 != p:
            raise exceptions.FilesError('number', i, c)

        stat.close()

    return "All checks pass!"


def write_files(dataset, dd, indir, outdir):
    """
    Writing established diagnoses and patients order into a file.
    :param dd: (dict) keys - IDs of patients, values - diagnoses
    :return: (str) Info about written diagnoses
    """
    pid_files = [os.path.join(indir, el) for el in os.listdir(indir) if
                 el.startswith('pid_chr') and el.endswith('.txt')]
    try:
        o = open(pid_files[0], 'r')
        fp = open("%spid_chr.txt" % outdir, 'w')
    except FileNotFoundError:
        try:
            o = open("%spid_chr.txt" % outdir, 'r')
        except FileNotFoundError:
            raise exceptions.NoFileError('pid_chr1.txt')

    fd = open("%sdiagnoses.txt" % outdir, 'w')

    n = 0
    i = 0
    for line in o:
        if 'fp' in locals():
            fp.write(line)
        try:
            if dataset == 'adni':
                pid = line.strip().split('_')[-1]
            else:
                pid = line.strip()
            fd.write('%s\n' % dd[pid])
            i += 1
        except KeyError:
            fd.write('NN\n')
            n += 1

    o.close()
    if 'fp' in globals():
        fp.close()
    fd.close()

    return 'Number of written diagnoses: %d (NN=%d)' % (i+n, n)


dir = './'
for q in range(len(sys.argv)):
    if sys.argv[q] == '-dir':
        dir = sys.argv[q+1]
    if sys.argv[q] == '-indir':
        indir = sys.argv[q+1]
    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q+1]
    if sys.argv[q] == '-dataset':
        dataset = sys.argv[q+1]
    if sys.argv[q] == '-diagdir':
        diagdir = sys.argv[q+1]

if 'dataset' not in globals():
    raise exceptions.NoParameterError('dataset', 'e.g. adni or rosmap')

if 'diagdir' not in globals():
    diagdir = '%sfiles/' % dir

if 'indir' not in globals():
    indir = '%smatrices/' % dir

if 'outdir' not in globals():
    outdir = '%smatrices/' % dir

if dataset == 'test':
    dfiles = ['test_diagnoses.tsv']
    dd = test_mapping(dfiles, diagdir)
if dataset == 'adni':
    dfiles = ['dxsum.csv']
    dd = adni_mapping(dfiles, diagdir)
elif dataset == 'rosmap':
    dfiles = ['diagnoses_Mayo.txt', 'diagnoses_MSBB.txt', 'diagnoses_Rosmap.txt']
    dd = rosmap_mapping(dfiles, diagdir)

print(check_pidfiles(indir))
print(write_files(dataset, dd, indir, outdir))
