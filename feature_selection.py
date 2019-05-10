import ast
import os
import boruta
from collections import Counter, OrderedDict, deque
import csv
import exceptions
import multiprocessing
import numpy as np
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score
import sys
import corporate_funcs as funcs
import argparse

'''
See readme.txt for input, output and possible options.
'''


def pooling(chrlist, classperc, dataset, outdir, pat, perc, r, run, snpsubset, snpruns, testpat, trainpat):
    """
    Running function one_process for every chr on chrlist.
    """
    procs = []

    ytrain, _ = build_y_matrices(dataset, run, outdir, pat, testpat, trainpat)

    for ch in chrlist:

        p = multiprocessing.Process(target=one_process,
                                    args=(ch, classperc, dataset, outdir, perc, r, run, snpsubset, snpruns,
                                          testpat, trainpat, ytrain))
        procs.append(p)
        p.start()

    for p in procs:
        p.join()

    return 0


def one_process(ch, classperc, dataset, outdir, perc, r, run, snpsubset, snpruns, testpat, trainpat, ytrain):
    """
    Loading data to matrices - function load_data.
    Selecting best SNPs subsets for every perc - function best_snps.
    Writing best SNPs for every chromosome into a file.
    """

    print('Analysis for chromosome %d started\n' % ch)

    Xtrain, Xtest, snp = funcs.load_data(ch, dataset, snpsubset, snpruns, testpat, trainpat)

    print('matrices X and y for chromosome %d have been loaded\n' % ch)

    snps = best_snps(perc, r, snp, Xtrain, ytrain)

    print('best SNPs for chromosome %d have been selected by Boruta\n' % ch)

    for p in classperc:
        np.save('%sX_train_chr%d_%d_%d.npy' % (outdir, ch, p, run), Xtrain[:, snps[p]])
        print('X train matrix for chr %d for perc %d was saved to file' % (ch, p))
        if testpat:
            np.save('%sX_test_chr%d_%d_%d.npy' % (outdir, ch, p, run), Xtest[:, snps[p]])

    for p in perc:
        lista = open('%sbestsnps_chr%d_%d_%d.txt' % (outdir, ch, p, run), 'w')
        lista.write('%d\n\n' % len(snps[p]))
        for el in snps[p]:
            lista.write('%d\n' % el)
        lista.close()

    print('process for chr %d finished\n' % ch)


def best_snps(perc, r, snp, X, y):
    s = snp // r
    snps = {a: [] for a in perc}

    for n in range(s + 1):

        if n != s:
            xx = X[:, n * r:n * r + r]
        elif n == s and snp % r != 0:
            xx = X[:, n * r:]

        for p in perc:
            result = run_boruta(xx, y, p)
            if not result:
                break
            else:
                snps[p] += [el + n * r for el in result]

    return snps


def run_boruta(X, y, p):
    rf = RandomForestClassifier(n_jobs=-1, class_weight='balanced', max_depth=5)
    feat_selector = boruta.BorutaPy(rf, n_estimators='auto', random_state=1, perc=p)
    feat_selector.fit(X, y)
    chosen = []
    for i, value in enumerate(feat_selector.support_):
        if value:
            chosen.append(i)
    return chosen


def build_y_matrices(dataset, run, outdir, pat, testpat, trainpat, testing=False):

    if not testing:
        if not testpat:
            try:
                with open('%stestpat_%d.txt' % (outdir, run), 'r') as ts:
                    testpat = set([int(el.strip()) for el in ts.readlines()])
            except FileNotFoundError:
                pass

    y_train = np.zeros(shape=(len(trainpat),), dtype=np.int8)
    if testpat:
        y_test = np.zeros(shape=(len(testpat),), dtype=np.int8)

    train_row = 0
    test_row = 0
    done = 0

    for name, directory in dataset.items():

        y = pd.read_csv('%smatrices/Y_chr.csv' % directory, header=None, index_col=0).values
        y = y.ravel()

        for i in range(pat[name]):
            if (done + i) in trainpat:
                y_train[train_row] = y[i]
                train_row += 1
            elif (done + i) in testpat:
                y_test[test_row] = y[i]
                test_row += 1

        done += i + 1

    if not testing:
        np.save('%sy_train_%d.npy' % (outdir, run), y_train)
        if testpat:
            np.save('%sy_test_%d.npy' % (outdir, run), y_test)

    if testpat:
        return y_train, y_test
    else:
        return y_train, None


def read_patlist(dataset, outdir, pat, patsubset, patruns, runs, testsize):

    for name, directory in dataset.items():
        if patsubset is not None:
            with open('%s%s/%s_patients_%d.txt' % (directory, patsubset, patsubset, patruns[name])) as file:
                patients = [int(line.strip()) for line in file.readlines()]
        else:
            patients = [i for i in range(pat[name])]
        if testsize != 0:
            testpat[name] = set(random.sample(patients, int(len(patients) * testsize)))
            trainpat[name] = set([p for p in patients if p not in testpat[name]])
            with open('%s/%s_testpat_%d.txt' % (name, outdir, runs[name]), 'w') as ts:
                ts.write('\n'.join([str(s) for s in sorted(testpat)]))
        else:
            testpat[name] = set()
            trainpat[name] = patients

    return testpat, trainpat


def read_boruta_params(chrlist, continuation, dataset, fixed, outdir, pat, run):

    for name, directory in dataset.items():
        file = '%sboruta_runs.txt' % (directory + 'boruta/')
        run_file = open(file, 'r')
        lines = run_file.readlines()
        towrite = ''
        occur = False
        for line in lines:
            if line.startswith(str(run) + '\t'):
                line = line.strip().split('\t')
                sets_order = line[1].strip().split('+')
                if list(dataset.keys()) != sets_order:
                    if len(dataset.keys()) != len(sets_order):
                        raise exceptions.WrongValueError('dataset', dataset,
                                                         'Other data sets were used in the given boruta run')
                    else:
                        for s in sets_order:
                            try:
                                dataset.move_to_end(s)
                            except KeyError:
                                raise exceptions.WrongValueError('dataset', dataset,
                                                                 'Data set named %s was not used in the given boruta run'
                                                                 % s)
                if line[3] == 'None':
                    patsubset = None
                else:
                    patsubset = line[3].strip()
                if line[4] == '-':
                    patruns = None
                else:
                    patruns = list(map(int, line[4].split('+')))
                    patruns = OrderedDict([(name, number) for name, number in zip(sets_order, patruns)])

                if line[5] == 'None':
                    snpsubset = None
                else:
                    snpsubset = line[5].strip()
                if line[6] == '-':
                    snpruns = None
                else:
                    snpruns = list(map(int, line[6].split('+')))
                    snpruns = OrderedDict([(name, number) for name, number in zip(sets_order, snpruns)])

                testsize = float(line[7])
                perc = list(map(int, line[8].split(',')))
                r = int(line[9])
                try:
                    outdir = line[11].strip()
                except IndexError:
                    pass
                if continuation:
                    line[10] = update_chrlist(fixed, line[10], chrlist)
                    strin = ''
                    for el in line:
                        strin += str(el) + '\t'
                    strin += '\n'
                    line = strin
                else:
                    chrlist = funcs.read_chrstr(line[10])
                occur = True
            if continuation:
                towrite += line
        run_file.close()

    if not occur:
        raise exceptions.WrongValueError('-run', str(run),
                                         'Boruta run number %d has not been conducted yet.' % run)

    patients = set()
    if patsubset is not None:
        done = 0
        for name in dataset.keys():
            with open('%s%s/%s_patients_%d.txt' % (dataset[name], patsubset, patsubset, patruns[name])) as file:
                selected = [int(line.strip()) for line in file.readlines()]
            patients = patients.union([pp + done for pp in selected])
            done += pat[name]
    else:
        patients = set([i for i in range(sum(pat.values()))])

    if testsize != 0:
        with open('%stestpat_%d.txt' % (outdir, run), 'r') as ts:
            testpat = set([int(el.strip()) for el in ts.readlines()])
        trainpat = set([p for p in patients if p not in testpat])
    else:
        testpat = set()
        trainpat = patients

    return chrlist, dataset, patruns, perc, r, snpsubset, snpruns, testpat, testsize, towrite, trainpat


def update_chrlist(fixed, linechrs, chrlist):

    chrs = funcs.read_chrstr(linechrs) + chrlist
    for key, value in Counter(chrs).items():
        if value > 1:
            if not fixed:
                print("WARNING: chromosome %d has already been analysed in this run, so it was omitted. " % key +
                      "If you want to analyse it anyway, please add '-fixed' attribute")
                chrlist.remove(key)
                if not chrlist:
                    raise exceptions.WrongValueError('chrlist', chrlist,
                                                     'There are no chromosomes to analyze!!!')
    chrs = list(set(chrs))
    chrs.sort()
    return funcs.make_chrstr(chrs)


perc = [90]
perc_given = False
r = 5000
chrlist = [i for i in range(1, 24)]

dataset = OrderedDict()
testsize = 0.1
testsize_given = False
runs = None
fixed = False
patsubset = None
patruns = None
snpsubset = None
snpruns = None
continuation = False


for q in range(len(sys.argv)):

    if sys.argv[q] == '-dataset':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            dataset[sys.argv[q+1]] = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')
        continue

    if sys.argv[q] == '-test':
        testsize = float(sys.argv[q + 1])
        testsize_given = True
        continue

    if sys.argv[q] == '-perc':
        perc = list(map(int, sys.argv[q+1].split(',')))
        perc_given = True
        continue

    if sys.argv[q] == '-snpsubset':
        snpsubset = sys.argv[q + 1]
        continue

    if sys.argv[q] == '-snprun':
        if sys.argv[q + 1] in dataset.keys():
            if snpruns is None:
                snpruns = OrderedDict()
            try:
                snpruns[sys.argv[q+1]] = int(sys.argv[q+2])
            except ValueError:
                raise exceptions.NoParameterError('snprun',
                                                  'After name of data set should appear number of SNP subset run')
        else:
            try:
                ss = int(sys.argv[q+1])
            except ValueError:
                raise exceptions.NoParameterError('snprun',
                                                  'After -snprun should appear name of data set and its run number ' +
                                                  'or one run number which is the same for every data set.')
        continue

    if sys.argv[q] == '-patsubset':
        patsubset = sys.argv[q + 1]
        continue

    if sys.argv[q] == '-patrun':
        if sys.argv[q + 1] in dataset.keys():
            if patruns is None:
                patruns = OrderedDict()
            try:
                patruns[sys.argv[q+1]] = int(sys.argv[q+2])
            except ValueError:
                raise exceptions.NoParameterError('patrun',
                                                  'After name of data set should appear number of patients subset run')
        else:
            try:
                pp = int(sys.argv[q+1])
            except ValueError:
                raise exceptions.NoParameterError('patrun',
                                                  'After -patrun should appear name of data set and its run number ' +
                                                  'or one run number which is the same for every data set.')
        continue

    if sys.argv[q] == '-chr':
        chrlist = funcs.read_chrstr(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-run':
        if sys.argv[q + 1] in dataset.keys():
            if runs is None:
                runs = OrderedDict()
            try:
                runs[sys.argv[q + 1]] = int(sys.argv[q + 2])
            except ValueError:
                raise exceptions.NoParameterError('run',
                                                  'After name of data set should appear number of boruta run')
        else:
            try:
                rr = int(sys.argv[q + 1])
            except ValueError:
                raise exceptions.NoParameterError('run',
                                                  'After -run should appear name of data set and its run number ' +
                                                  'or one run number which is the same for every data set.')
        continue

    if sys.argv[q] == '-fixed':
        fixed = True
        continue

    if sys.argv[q] == '-outdir':
        outdir = sys.argv[q + 1]
        continue

    if sys.argv[q] == '-cont':
        continuation = True
        continue

    if sys.argv[q].startswith('-'):
        raise exceptions.WrongParameterName(sys.argv[q])


if 'outdir' not in globals():
    if len(dataset) == 1:
        outdir = next(iter(dataset.values())) + 'boruta/'
    elif len(dataset) > 1:
        raise exceptions.NoParameterError('outdir', 'There is more than one dataset, but outdir was not given!')
    else:
        raise exceptions.NoParameterError('outdir', 'No dataset or outdir was given!')

if not os.path.isdir(outdir):
    os.mkdir(outdir)

if perc_given:
    perc.sort()

if 'ss' in globals():
    if snpruns is None:
        snpruns = OrderedDict([(n, num) for n, num in zip(dataset.keys(), [ss]*len(dataset))])
    else:
        for name in [nn for nn in dataset.keys() if nn not in snpruns.keys()]:
            snpruns[name] = ss

if 'pp' in globals():
    if patruns is None:
        patruns = OrderedDict([(n, num) for n, num in zip(dataset.keys(), [pp]*len(dataset))])
    else:
        for name in [nn for nn in dataset.keys() if nn not in patruns.keys()]:
            patruns[name] = pp

if 'rr' in globals():
    if runs is None:
        runs = OrderedDict([(n, num) for n, num in zip(dataset.keys(), [rr]*len(dataset))])
    else:
        for name in [nn for nn in dataset.keys() if nn not in runs.keys()]:
            runs[name] = rr

# determination number of patient in given data sets
pat = funcs.patients(dataset)

if not continuation:
    for name, directory in dataset:
        runs[name] = funcs.establish_run('boruta', fixed, outdir, runs[name])
    testpat, trainpat = read_patlist(dataset, outdir, pat, patsubset, patruns, testsize)

else:
    dataset, outdir, patruns, perc, r, snpsubset, snpruns, testpat, testsize, towrite, trainpat = \
        read_boruta_params(chrlist, continuation, dataset, fixed, pat, runs)

    # running Boruta analysis
    pooling(chrlist, dataset, outdir, pat, perc, r, borutarun, snpsubset, snpruns, testpat, trainpat)

    # saving information about done run to boruta_runs file
    if not continuation:
        if patruns is None:
            patruns_string = '-'
        else:
            patruns_string = '+'.join(list(map(str, patruns.values())))
        if snpruns is None:
            snpruns_string = '-'
        else:
            snpruns_string = '+'.join(list(map(str, snpruns.values())))

        funcs.runs_file_add('boruta', outdir, borutarun, '%d\t%s\t%d\t%s\t%s\t%s\t%s\t%.2f\t%s\t%d\t%s\n' %
                            (borutarun, '+'.join(dataset.keys()), len(trainpat) + len(testpat), patsubset,
                             patruns_string, snpsubset, snpruns_string, testsize, ','.join(list(map(str, perc))), r,
                             funcs.make_chrstr(chrlist)))
    else:
        funcs.runs_file_rewrite('boruta', outdir, towrite)