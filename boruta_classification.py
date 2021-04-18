import ast
import os
import boruta
from collections import Counter, OrderedDict, deque
import math
import exceptions
import multiprocessing
import numpy as np
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
from sklearn.tree import DecisionTreeClassifier
from sklearn.linear_model import LogisticRegression
from sklearn.model_selection import KFold
from sklearn.metrics import roc_auc_score
import sys
import corporate_funcs as funcs

'''
See readme.txt for input, output and possible options.
'''


def pooling(num_cores, chrlist, classperc, dataset, outdir, pat, perc, r, run, snpsubset, snpruns, testpat, trainpat):
    """
    Running function one_process for every chr on chrlist.
    """
    if num_cores is None:
        num_cores = len(chrlist)

    ytrain, _ = build_y_matrices(dataset, run, outdir, pat, testpat, trainpat)

    for chrlist_subset in [chrlist[i*num_cores:(i+1)*num_cores] for i in range(math.ceil(len(chrlist)/num_cores))]:

        procs = []
        for ch in chrlist_subset:

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


def read_typedata(chrlist, outdir, p, run, type):

    ch = chrlist[0]
    X = np.load('%sX_%s_chr%d_%d_%d.npy' % (outdir, type, ch, p, run))

    for ch in chrlist[1:]:

        X = np.concatenate((X, np.load('%sX_%s_chr%d_%d_%d.npy' % (outdir, type, ch, p, run))), axis=1)

    y = np.load('%sy_%s_%d.npy' % (outdir, type, run))

    return X, y


def build_data(borutarun, chrlist, classrun, dataset, frombed, newforest_notcv, outdir, p, snpsubset, snpruns, testset, testsize):

    patients = sum(funcs.patients(testset).values())
    if (newforest_notcv or frombed) and testsize != 0:
        case, control = funcs.patients_diagnoses(dataset, set([i for i in range(patients)]))
        half = round(patients * testsize) / 2
        testpat = set(random.sample(case, max(math.floor(half), 1)) + random.sample(control, max(math.ceil(half), 1)))
        with open('%stestpat_class_%d.txt' % (outdir, classrun), 'w') as ts:
            ts.write('\n'.join([str(s) for s in sorted(testpat)]))
        trainpat = [i for i in range(patients) if i not in testpat]
    else:
        testpat = None
        trainpat = list(range(patients))

    X, X_test = None, None
    for ch in chrlist:
        if frombed:
            selected_snps = {
                next(iter(testset.keys())):
                    [int(el) for el in
                     open(os.path.join(next(iter(dataset.values())), 'frombed', 'frombed_snps_chr{}_{}.txt'.
                                       format(ch, borutarun)), 'r').read().strip().split('\n')]
            }
        else:
            selected_snps = read_selected_snps(ch, dataset, frombed, outdir, p, borutarun, snpsubset, snpruns, testset)

        xx, xx_test, snp = funcs.read_Xs(ch, testset, len(next(iter(selected_snps.values()))), selected_snps, testpat, trainpat)

        if X is None:
            X = xx.copy()
            if xx_test is not None:
                X_test = xx_test.copy()
        else:
            X = np.concatenate((X, xx), axis=1)
            if xx_test is not None:
                X_test = np.concatenate((X_test, xx_test), axis=1)

    y, y_test = build_y_matrices(testset, None, outdir, pat, testpat, trainpat, testing=True)
    return X, y, X_test, y_test, patients


def classify(X, y, X_test, y_test, cv, classifier=None):

    if classifier == 'tree':
        rf = DecisionTreeClassifier()
    elif classifier == 'logreg':
        rf = LogisticRegression()
    else:
        rf = RandomForestClassifier(n_estimators=500)

    if not cv:
        rf.fit(X, y)
        prob = rf.predict_proba(X_test)
        order = [i for i, c in enumerate(rf.classes_) if c == 1]
        y_score = prob[:, order]
        return rf.score(X, y), rf.score(X_test, y_test), roc_auc_score(y_test, y_score)
    else:
        kf = KFold(n_splits=cv)
        scores = [[], [], []]
        for train_index, test_index in kf.split(X):
            rf.fit(X[train_index], y[train_index])
            prob = rf.predict_proba(X[test_index])
            order = [i for i, c in enumerate(rf.classes_) if c == 1]
            y_score = prob[:, order]
            try:
                [s.append(el) for s, el in
                 zip(scores, [rf.score(X[train_index], y[train_index]), rf.score(X[test_index], y[test_index]),
                              roc_auc_score(y[test_index], y_score)])]
            except ValueError:
                scores[0].append(rf.score(X[train_index], y[train_index]))
                scores[1].append(rf.score(X[test_index], y[test_index]))
        return list(map(np.mean, scores)), list(map(np.std, scores))


def classify_cv_both(X_dataset, y_dataset, X_testset, y_testset, cv=10, classifier=None):

    if classifier == 'tree':
        rf = DecisionTreeClassifier()
    elif classifier == 'logreg':
        rf = LogisticRegression()
    else:
        rf = RandomForestClassifier(n_estimators=500)

    kf = KFold(n_splits=cv)
    scores = [[], [], []]
    scores_testset = [[], []]
    for train_index, test_index in kf.split(X_dataset):
        rf.fit(X_dataset[train_index], y_dataset[train_index])
        order = [i for i, c in enumerate(rf.classes_) if c == 1]

        prob = rf.predict_proba(X_dataset[test_index])
        y_score = prob[:, order]
        scores[0].append(rf.score(X_dataset[train_index], y_dataset[train_index]))
        scores[1].append(rf.score(X_dataset[test_index], y_dataset[test_index]))
        try:
            scores[2].append(roc_auc_score(y_dataset[test_index], y_score))
        except ValueError:
            pass

        prob = rf.predict_proba(X_testset)
        y_score = prob[:, order]
        scores_testset[0].append(rf.score(X_testset, y_testset))
        try:
            scores_testset[1].append(roc_auc_score(y_testset, y_score))
        except ValueError:
            pass
    return list(map(np.mean, scores)), list(map(np.std, scores)), \
           list(map(np.mean, scores_testset)), list(map(np.std, scores_testset))


def first_run(dataset, fixed, outdir, pat, patsubset, patruns, run, testsize):

    run = funcs.establish_run('boruta', fixed, outdir, run)

    patients = set()
    if patsubset is not None:
        done = 0
        for name in dataset.keys():
            with open('%s%s/%s_patients_%d.txt' % (dataset[name], patsubset, patsubset, patruns[name])) as file:
                selected = [int(line.strip()) for line in file.readlines()]
            patients = patients.union([p+done for p in selected])
            done += pat[name]
    else:
        patients = set([i for i in range(sum(pat.values()))])

    case, control = funcs.patients_diagnoses(dataset, patients)
    if testsize != 0:
        half = round(len(patients) * testsize) / 2
        testpat = set(random.sample(case, max(math.floor(half), 1)) + random.sample(control, max(math.ceil(half), 1)))
        trainpat = set([p for p in patients if p not in testpat])
        with open('%stestpat_%d.txt' % (outdir, run), 'w') as ts:
            ts.write('\n'.join([str(s) for s in sorted(testpat)]))
    else:
        testpat = set()
        trainpat = patients

    return run, testpat, trainpat


def read_boruta_params(chrlist, continuation, dataset, fixed, outdir, pat, run):

    file = '%sboruta_runs.txt' % outdir
    funcs.correct_boruta_runs_file(file)
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
            r = int(line[-2])
            if continuation:
                line = update_chrlist(fixed, line, chrlist)
                strin = ''
                for el in line:
                    strin += str(el) + '\t'
                strin += '\n'
                line = strin
            else:
                chrlist = funcs.read_chrstr(line[-1])
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


def update_chrlist(fixed, line, chrlist):

    chrs = funcs.read_chrstr(line[-1]) + chrlist
    for key, value in Counter(chrs).items():
        if value > 1:
            if not fixed:
                print("WARNING: chromosome %d has already been analysed in this run, so it was omited. " % key +
                      "If you want to analyse it anyway, please add '-fixed' attribute")
                chrlist.remove(key)
                if not chrlist:
                    raise exceptions.WrongValueError('chrlist', chrlist,
                                                     'There are no chromosomes to analyze!!!')
    chrs = list(set(chrs))
    chrs.sort()
    line[-1] = funcs.make_chrstr(chrs)

    return line


def read_selected_snps(ch, dataset, frombed, outdir, p, run, snpsubset, snpruns, testset):

    selected_snps = deque()

    if frombed:
        bestfile = open('%sfrombed/frombed_snps_chr%d_%d.txt' % (next(iter(dataset.values())), ch, run))
    else:
        bestfile = open('%sbestsnps_chr%d_%d_%d.txt' % (outdir, ch, p, run), 'r')
        for _ in range(2):  # header
            bestfile.readline()
    if snpsubset is not None:
        snprun = next(iter(snpruns.values()))
        subsetfile = open('%s%s/%s_snps_chr%d_%d.txt' % (next(iter(dataset.values())), snpsubset, snpsubset, ch, snprun), 'r')
        best, last = best_from_subset(bestfile, subsetfile, -1)
    else:
        best = int(bestfile.readline())
    done = 0
    for directory in dataset.values():
        snpfile = open('%smatrices/snps_chr%d.txt' % (directory, ch))
        for i, line in enumerate(snpfile):
            if i == (best - done):
                selected_snps.append(int(line.split()[0]))
                try:
                    if snpsubset is not None:
                        best, last = best_from_subset(bestfile, subsetfile, last)
                    else:
                        best = int(bestfile.readline())
                except ValueError:
                    best = 0
        snpfile.close()
        done += i
    ssnp = selected_snps.popleft()
    snps_test = {s: [] for s in testset.keys()}
    for name, directory in testset.items():
        snpfile = open('%smatrices/snps_chr%d.txt' % (directory, ch))
        for i, line in enumerate(snpfile):
            if int(line.split()[0]) == ssnp:
                snps_test[name].append(i)
                try:
                    ssnp = selected_snps.popleft()
                except IndexError:
                    break
        snpfile.close()

    return snps_test


def best_from_subset(bestfile, subsetfile, last):

    b = int(bestfile.readline())
    for _ in range(b-last-1):
        subsetfile.readline()
    return int(subsetfile.readline()), b


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


perc = [90]
perc_given = False
r = 5000
chrlist = [i for i in range(1, 24)]

dataset = OrderedDict()
testset = OrderedDict()
testsize = 0.1
testsize_given = False
class_only = False
boruta_only = False
borutarun = None
classrun = None
frombedrun = None
run = None
fixed = False
patsubset = None
patruns = None
snpsubset = None
snpruns = None
continuation = False
makey = False
cv = None
newforest = False
frombed = False
num_cores = None
method = 'rforest'

for q in range(len(sys.argv)):

    if sys.argv[q] == '-dataset':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            dataset[sys.argv[q+1]] = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of data set should appear a directory to folder with it.')
        continue

    if sys.argv[q] == '-testset':
        if sys.argv[q + 2][0] in ['.', '~', '/']:
            testset[sys.argv[q+1]] = sys.argv[q + 2]
        else:
            raise exceptions.NoParameterError('directory',
                                              'After name of test set should appear a directory to folder with it.')
        continue

    if sys.argv[q] == '-test':
        testsize = float(sys.argv[q + 1])
        testsize_given = True
        continue

    if sys.argv[q] == '-perc':
        perc = list(map(int, sys.argv[q+1].split(',')))
        perc_given = True
        continue

    if sys.argv[q] == '-classperc':
        classperc = list(map(int, sys.argv[q+1].split(',')))
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

    if sys.argv[q] == '-r':
        r = int(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-class':
        class_only = True
        continue

    if sys.argv[q] == '-boruta':
        boruta_only = True
        continue

    if sys.argv[q] == '-chr':
        chrlist = funcs.read_chrstr(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-run':
        run = int(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-borutarun':
        borutarun = int(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-classrun':
        classrun = int(sys.argv[q + 1])
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

    if sys.argv[q] == '-makeY':
        makey = True
        continue

    if sys.argv[q] == '-cv':
        cv = int(sys.argv[q+1])
        testsize = 0
        continue

    if sys.argv[q] == '-newforest':
        newforest = True
        continue

    if sys.argv[q] == '-frombed':
        frombed = True
        continue

    if sys.argv[q] == '-frombedrun':
        frombedrun = int(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-num_cores':
        num_cores = int(sys.argv[q + 1])
        continue

    if sys.argv[q] == '-method':
        method = sys.argv[q + 1]
        continue

    if sys.argv[q].startswith('-'):
        raise exceptions.WrongParameterName(sys.argv[q])


if 'outdir' not in globals():
    if dataset:
        outdir = next(iter(dataset.values())) + 'boruta/'
    elif testset:
        outdir = next(iter(testset.values())) + 'boruta/'
    else:
        raise exceptions.NoParameterError('outdir', 'No dataset, testset or outdir was given!')
    if not os.path.isdir(outdir):
        os.mkdir(outdir)


if 'classperc' in globals():
    classperc.sort()
if perc_given:
    perc.sort()
    if 'classperc' not in globals():
        classperc = perc

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

# determination number of patient in given data set
pat = funcs.patients({el: la for el, la in zip(list(dataset.keys())+list(testset.keys()),
                                               list(dataset.values())+list(testset.values()))})

if not class_only:

    if borutarun is None and run is not None:
        borutarun = run

    # determination of some parameters
    if not continuation:
        borutarun, testpat, trainpat = first_run(dataset, fixed, outdir, pat, patsubset, patruns, borutarun, testsize)
    else:
        dataset, patruns, perc, r, snpsubset, snpruns, testpat, testsize, towrite, trainpat = \
            read_boruta_params(chrlist, continuation, dataset, fixed, outdir, pat, borutarun)

    if 'classperc' not in globals():
        classperc = perc
    # running Boruta analysis
    pooling(num_cores, chrlist, classperc, dataset, outdir, pat, perc, r, borutarun, snpsubset, snpruns, testpat, trainpat)

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

if not boruta_only:

    if classrun is None and run is not None:
        classrun = run

    # determination of number of class run
    classrun = funcs.establish_run('class', fixed, outdir, classrun)
    scores_file = open('%sclass_scores_%d.txt' % (outdir, classrun), 'w', 1)
    if frombed and cv:
        scores_file.write('perc\tSNPs\tdataset_train_score\tdataset_test_score\tdataset_AUC\ttestset_score\ttestset_AUC\n')
    else:
        scores_file.write('perc\tSNPs\ttrain_score\ttest_score\tAUC\n')  # writing heading to class_scores file

    if makey:
        build_y_matrices(dataset, borutarun, outdir, funcs.patients(dataset), testpat, trainpat)

    if cv and dataset and testset:
        assert frombed
        X_test, y_test, _, _, testpat_val = build_data(borutarun, chrlist, classrun, dataset, frombed, newforest,
                                                       outdir, None, None, None, testset, testsize)
        X_train, y_train = read_typedata(chrlist, outdir, perc[0], borutarun, 'train')
        print('X train shape: {}, X test shape: {}'.format(X_train.shape, X_test.shape))
        print('Data loaded!')
        if X_train.shape[1] > 0:
            ((score_train, score_test, score_auc), (std_train, std_test, std_auc),
             (score_test_testset, score_auc_testset), (std_score_testset, std_auc_testset)) = \
                classify_cv_both(X_train, y_train, X_test, y_test, cv, classifier=method)
            scores_file.write(
                '%s\t%d\t%.3f +- %.3f\t%.3f +- %.3f\t%.3f +- %.3f\t%.3f +- %.3f\t%.3f +- %.3f\n' %
                ('cv-frombed', X_train.shape[1], score_train, std_train, score_test, std_test, score_auc,
                 std_auc, score_test_testset, std_score_testset, score_auc_testset, std_auc_testset))

    if frombed:

        if not testset and dataset:
            testset = dataset

        if testsize != 0:
            X_train, y_train, X_test, y_test, testpat_val = build_data(frombedrun, chrlist, classrun, dataset, True,
                                                                       newforest, outdir, None, None, None, testset, testsize)
        else:
            X_test, y_test, _, _, testpat_val = build_data(borutarun, chrlist, classrun, dataset, True, newforest,
                                                           outdir, None, None, None, testset, testsize)
            X_train, y_train = read_typedata(chrlist, outdir, perc[0], borutarun, 'train')

        print('X train shape: {}, X test shape: {}'.format(X_train.shape, X_test.shape))
        print('Data loaded!')
        if X_train.shape[1] > 0:
            if cv:
                ((score_train, score_test, score_auc), (std_train, std_test, std_auc)) = \
                    classify(X_train, y_train, X_test, y_test, cv, classifier=method)
                scores_file.write(
                    '%s\t%d\t%.3f +- %.3f\t%.3f +- %.3f\t%.3f +- %.3f\n' % ('frombed', X_train.shape[1], score_train,
                                                                            std_train, score_test, std_test, score_auc,
                                                                            std_auc))
            else:
                score_train, score_test, score_auc = classify(X_train, y_train, X_test, y_test, cv, classifier=method)
                scores_file.write('%s\t%d\t%.3f\t%.3f\t%.3f\n' % ('frombed', X_train.shape[1], score_train, score_test,
                                                                  score_auc))
            # saving matrices (on which was based the classification) to file
            for name in ['X_train', 'y_train', 'X_test', 'y_test']:
                np.save('%s%s_genome_%s_%d.npy' % (outdir, name, 'frombed', classrun), eval(name))
        else:
            print('No SNPs were chosen!')
            scores_file.write('frombed\t0\t-\t-\n')

        print('Scores saved to file')

        scores_file.close()
        print('Classification done!')

        if cv:
            teststr = 'CV-%d: %s' % (cv, '+'.join(testset))
            trainpat_val = 0
        else:
            teststr = '+'.join(testset)
            trainpat_val = sum(funcs.patients(dataset).values())
        trainstr = '+'.join(dataset)

        funcs.runs_file_add('class', outdir, classrun, '%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n' %
                            (classrun, teststr, testpat_val, frombedrun, trainstr, trainpat_val,
                             'frombed', funcs.make_chrstr(chrlist), method))

    else:

        if dataset:
            chrlist, dataset, patruns, perc, r, snpsubset, snpruns, testpat, _, towrite, trainpat = \
                read_boruta_params(chrlist, False, dataset, False, outdir, pat, borutarun)
        elif testset:
            chrlist, testset, patruns, perc, r, snpsubset, snpruns, testpat, _, towrite, trainpat = \
                read_boruta_params(chrlist, False, testset, False, outdir, pat, borutarun)

        if 'classperc' not in globals():
            classperc = perc

        for p in classperc:

            # establishing testing and training data based on given test set(s) or test subset of patients
            if not cv and not newforest:
                if not testset:
                    if testsize == 0:
                        raise exceptions.NoParameterError(
                            'testset', 'Test size was not given - define what set should be used as a testset.')
                    else:
                        print('Standard classification after Boruta')
                        X_test, y_test = read_typedata(chrlist, outdir, p, borutarun, 'test')
                else:
                    print('Classification based on the testset and given forest')
                    X_test, y_test, _, _, testpat_val = build_data(borutarun, chrlist, classrun, dataset, False, newforest, outdir,
                                                                   p, snpsubset, snpruns, testset, testsize)
            elif cv:
                if newforest:
                    print('CV classification based on testset and list of selected SNPs')
                    X_train, y_train, _, _, _ = build_data(borutarun, chrlist, classrun, dataset, False, False, outdir, p,
                                                           snpsubset, snpruns, testset, 0)
                else:
                    print('CV classification based on earlier-built matrices for given set')
                X_test = y_test = None
            elif newforest:
                print('Classification based on testset and list of selected SNPs')
                X_train, y_train, X_test, y_test, testpat_val = build_data(borutarun, chrlist, classrun, dataset, False, True,
                                                                           outdir, p, snpsubset, snpruns, testset, testsize)

            if not newforest:
                X_train, y_train = read_typedata(chrlist, outdir, p, borutarun, 'train')

            # running classification and saving scores to class_scores file
            print('Data loaded!')
            if X_train.shape[1] > 0:
                if cv:
                    ((score_train, score_test, score_auc), (std_train, std_test, std_auc)) = \
                        classify(X_train, y_train, X_test, y_test, cv, classifier=method)
                    scores_file.write('%s\t%d\t%.3f +- %.3f\t%.3f +- %.3f\t%.3f +- %.3f\n' %
                                      (p, X_train.shape[1], score_train, std_train, score_test, std_test, score_auc,
                                       std_auc))
                else:
                    score_train, score_test, score_auc = \
                        classify(X_train, y_train, X_test, y_test, cv, classifier=method)
                    scores_file.write(
                        '%s\t%d\t%.3f\t%.3f\t%.3f\n' % (p, X_train.shape[1], score_train, score_test,
                                                        score_auc))
                # saving matrices (on which was based the classification) to file
                for name in ['X_train', 'y_train', 'X_test', 'y_test']:
                    np.save('%s%s_genome_%d_%d.npy' % (outdir, name, p, classrun), eval(name))
            else:
                print('No SNPs were chosen for perc %d' % p)
                scores_file.write('%d\t0\t-\t-\n' % p)

            print('Scores for perc %d saved to file' % p)

        scores_file.close()
        print('Classification done!')

        # writing information about class run to class_run file
        if newforest:
            trainstr = 'only-SNPs: %s' % '+'.join(dataset.keys())
        else:
            trainstr = '+'.join(dataset.keys())
        if 'testpat_val' not in globals():
            testpat_val = len(testpat)
        if 'trainpat' in globals():
            trainpat_val = len(trainpat)
        if cv:
            teststr = 'CV-%d: %s' % (cv, trainstr)
            testpat_val = 0
            trainpat_val = sum(pat.values())
        elif not testset:
            teststr = '%.2f*(%s)' % (testsize, trainstr)
        else:
            teststr = '+'.join(testset.keys())

        funcs.runs_file_add('class', outdir, classrun, '%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n' %
                            (classrun, teststr, testpat_val, borutarun, trainstr, trainpat_val,
                             ','.join(list(map(str, classperc))), funcs.make_chrstr(chrlist), method))
