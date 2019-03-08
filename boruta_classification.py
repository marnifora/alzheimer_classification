import ast
import boruta
from collections import Counter, OrderedDict
import csv
import exceptions
import multiprocessing
import numpy as np
import pandas as pd
import random
from sklearn.ensemble import RandomForestClassifier
import sys
import corporate_funcs as funcs

'''
See readme.txt for input, output and possible options.
'''


def pooling(chrlist, classperc, dataset, outdir, pat, perc, r, run, snpsubset, snpruns, testpat, trainpat):
    """
    Running function one_process for every chr on chrlist.
    """
    procs = []

    ytrain = build_y_matrices(dataset, run, outdir, pat, testpat, trainpat)

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

    Xtrain, Xtest, snp = load_data(ch, dataset, snpsubset, snpruns, testpat, trainpat)

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


def load_data(ch, dataset, snpsubset, snpruns, testpat, trainpat):
    """
    Loading data from files into X and y matrices.
    Selection of patients not being in test set and SNPs present in subset-file.
    """

    snplist = {name: [] for name in dataset.keys()}
    if snpsubset is not None:
        for name in dataset.keys():
            cc = open('%s%s/%s_snps_chr%d_%d.txt' % (dataset[name], snpsubset, snpsubset, ch, snpruns[name]), 'r')
            for line in cc:
                snplist[name].append(int(line.split()[0]))
            cc.close()
        snp = len(snplist[name])
    else:
        if len(dataset) > 1:
            raise exceptions.NoParameterError('subset', 'There is more than one given data set, but subset of SNPs ' +
                                                        'is not given.')
        else:
            cc = open('%smatrices/genome_stats.txt' % list(dataset.values())[0], 'r')
            snp = None
            for line in cc:
                if line.startswith('%d\t' % ch):
                    snp = int(line.split()[1])
                    break
            cc.close()
            if 'snp' is None:
                raise exceptions.OtherError('There is no information about chromosome %d in %sgenome_stats.txt file'
                                            % (ch, list(dataset.values())[0]))
            snplist[next(iter(dataset.keys()))] = list(range(snp))

    train_row = 0
    test_row = 0
    done = 0

    X_train = np.zeros(shape=(len(trainpat), snp), dtype=np.int8)

    X_test = np.zeros(shape=(len(testpat), snp), dtype=np.int8)

    for name in dataset.keys():

        o = open('%smatrices/X_chr%d_nodif.csv' % (dataset[name], ch), 'r')
        reader = csv.reader(o, delimiter=',')
        next(reader)  # header

        # writing values from file to matrix X

        for i, line in enumerate(reader):
            if (done + i) in trainpat:
                for j, s in enumerate(snplist[name]):
                    X_train[train_row][j] = line[s + 1]
                train_row += 1
            elif (done + i) in testpat:
                for j, s in enumerate(snplist[name]):
                    X_test[test_row][j] = line[s + 1]
                test_row += 1

        o.close()
        done += i + 1

    if testpat:
        return X_train, X_test, snp
    else:
        return X_train, None, snp


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


def build_testdata(chrlist, selected_snps, testset):

    pat = patients(testset)
    for name in testset.keys():

        xx = np.zeros(shape=(pat[name], sum([len(x) for x in [selected_snps[ch] for ch in chrlist]])), dtype=np.int8)
        col = 0
        for ch in chrlist:
            o = open('%smatrices/X_chr%d_nodif.csv' % (testset[name], ch), 'r')
            reader = csv.reader(o, delimiter=',')
            next(reader)

            snps = selected_snps[ch]
            for i, line in enumerate(reader):
                xx[i][col:col+len(snps)] = [line[1:][ii] for ii in snps]

            col += len(snps)
            o.close()

        yy = pd.read_csv('%smatrices/Y_chr.csv' % testset[name], header=None, index_col=0).values

        if 'X_test' not in locals() or 'y_test' not in locals():
            X_test = xx
            y_test = yy
        else:
            X_test = np.concatenate((X_test, xx), axis=0)
            y_test = np.concatenate((y_test, yy), axis=0)

    return X_test, y_test


def classify(X_train, y_train, X_test, y_test):

    rf = RandomForestClassifier(n_estimators=500)
    rf.fit(X_train, y_train)
    return rf.score(X_train, y_train), rf.score(X_test, y_test)


def first_run(dataset, fixed, outdir, pat, patsubset, patruns, run, testsize):

    run = funcs.establish_run('boruta', fixed, outdir, run)

    patients = []
    if patsubset is not None:
        done = 0
        for name in dataset.keys():
            with open('%s%s/%s_patients_%d.txt' % (dataset[name], patsubset, patsubset, patruns[name])) as file:
                selected = [int(line.strip()) for line in file.readlines()]
            patients += [p+done for p in selected]
            done += pat[name]
    else:
        patients = [i for i in range(sum(pat.values()))]

    if testsize != 0:
        testpat = random.sample(patients, int(len(patients)*testsize))
        trainpat = [p for p in patients if p not in testpat]
        ts = open('%stestpat_%d.txt' % (outdir, run), 'w')
        for el in testpat:
            ts.write('%d\n' % el)
        ts.close()
    else:
        testpat = None
        trainpat = patients

    return run, testpat, trainpat


def read_boruta_params(chrlist, continuation, dataset, fixed, outdir, pat, run):

    file = '%sboruta_runs.txt' % outdir
    funcs.correct_boruta_runs_file(file)
    run_file = open(file, 'r+')
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
            if line[4] == '-':
                patruns = None
            else:
                patruns = list(map(int, line[4].split('+')))
                patruns = OrderedDict([(name, number) for name, number in zip(sets_order, patruns)])

            if line[5] == 'None':
                snpsubset = None
            if line[6] == '-':
                snpruns = None
            else:
                snpruns = list(map(int, line[6].split('+')))
                snpruns = OrderedDict([(name, number) for name, number in zip(sets_order, snpruns)])

            testsize = float(line[7])
            perc = ast.literal_eval(line[8])
            if isinstance(perc, int):
                perc = [perc]
            r = int(line[-2])
            if continuation:
                line = update_chrlist(fixed, line, chrlist)
                strin = ''
                for el in line:
                    strin += str(el) + '\t'
                strin += '\n'
                line = strin
            occur = True
        if continuation:
            towrite += line
    run_file.close()

    if not occur:
        raise exceptions.WrongValueError('-run', str(run),
                                         'Boruta run number %d has not been conducted yet.' % run)

    patients = []
    if patsubset is not None:
        done = 0
        for name in dataset.keys():
            with open('%s%s/%s_patients_%d.txt' % (dataset[name], patsubset, patsubset, patruns[name])) as file:
                selected = [int(line.strip()) for line in file.readlines()]
            patients += [pp + done for pp in selected]
            done += pat[name]
    else:
        patients = [i for i in range(sum(pat.values()))]

    if testsize != 0:
        ts = open('%stestpat_%d.txt' % (outdir, run), 'r')
        testpat = []
        for line in ts:
            testpat.append(int(line.strip()))
        ts.close()
        trainpat = [p for p in patients if p not in testpat]
    else:
        testpat = None
        trainpat = patients

    return dataset, patruns, perc, r, snpsubset, snpruns, testpat, testsize, towrite, trainpat


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


def patients(dataset):

    pat = {name: 0 for name in dataset.keys()}
    for name in dataset.keys():
        g = open('%smatrices/genome_stats.txt' % dataset[name], 'r')
        line = g.readline()
        p = int(line.split()[3])
        for line in g:
            if int(line.split()[3]) != p:
                raise exceptions.OtherError('Error: there is different number of patients for different chromosomes!')
        pat[name] = p
        g.close()
    return pat


def read_selected_snps(chrlist, outdir, p, run):

    selected_snps = {ch: [] for ch in chrlist}
    for ch in chrlist:
        o = open('%sbestsnps_chr%d_%d_%d.txt' % (outdir, ch, p, run), 'r')
        for i in range(2):  # header
            o.readline()
        for line in o:
            selected_snps[ch].append(int(line.strip()))
        o.close()

    return selected_snps


def build_y_matrices(dataset, run, outdir, pat, testpat, trainpat):

    if testpat is None:
        testpat = []
        try:
            ts = open('%stestpat_%d.txt' % (outdir, run), 'r')
            for line in ts:
                testpat.append(int(line.strip()))
            ts.close()
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

    np.save('%sy_train_%d.npy' % (outdir, run), y_train)
    if testpat:
        np.save('%sy_test_%d.npy' % (outdir, run), y_test)

    return y_train


perc = [90]
r = 5000
chrlist = [i for i in range(1, 24)]

dataset = OrderedDict()
testset = OrderedDict()
testsize = 0
class_only = False
boruta_only = False
borutarun = None
classrun = None
fixed = False
patsubset = None
patruns = None
snpsubset = None
snpruns = None
continuation = False
makey = False

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
        continue

    if sys.argv[q] == '-perc':
        perc = list(map(int, sys.argv[q+1].split(',')))
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
        borutarun = int(sys.argv[q + 1])
        classrun = int(sys.argv[q + 1])
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

    if sys.argv[q].startswith('-'):
        raise exceptions.WrongParameterName(sys.argv[q])

perc.sort()

if 'outdir' not in globals():
    outdir = next(iter(dataset.values())) + 'boruta/'

if 'classperc' not in globals():
    classperc = perc
else:
    classperc.sort()

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
pat = patients(dataset)

if not class_only:

    # determination of some parameters
    if not continuation:
        borutarun, trainpat, testpat = first_run(dataset, fixed, outdir, pat, patsubset, patruns, borutarun, testsize)
    else:
        dataset, patruns, perc, r, snpsubset, snpruns, testpat, testsize, towrite, trainpat = \
            read_boruta_params(chrlist, continuation, dataset, fixed, outdir, pat, borutarun)

    # running Boruta analysis
    pooling(chrlist, classperc, dataset, outdir, pat, perc, r, borutarun, snpsubset, snpruns, testpat, trainpat)

    # saving information about done run to boruta_runs file
    if not continuation:
        run_file = open('%sboruta_runs.txt' % outdir, 'a')
        run_file.write('%d\t%s\t%d\t%s\t%s\t%s\t%s\t%.1f\t%s\t%d\t%s\n' % (borutarun, '+'.join(dataset.keys()),
                                                                           len(trainpat) + len(testpat), patsubset,
                                                                           '+'.join(list(map(str, patruns.values()))),
                                                                           snpsubset,
                                                                           '+'.join(list(map(str, snpruns.values()))),
                                                                           testsize, ','.join(list(map(str, perc))), r,
                                                                           funcs.make_chrstr(chrlist)))
    else:
        run_file = open('%sboruta_runs.txt' % outdir, 'w')
        run_file.write(towrite)

    run_file.close()


if not boruta_only:

    # determination of number of class run
    classrun = funcs.establish_run('class', fixed, outdir, classrun)
    dataset, patruns, perc, r, snpsubset, snpruns, testpat, testsize, towrite, trainpat = \
        read_boruta_params(chrlist, False, dataset, False, outdir, pat, borutarun)
    scores_file = open('%sclass_scores_%d.txt' % (outdir, classrun), 'w', 1)
    scores_file.write('perc\tSNPs\ttrain_score\ttest_score\n')  # writing heading to class_scores file
    if makey:
        build_y_matrices(dataset, borutarun, outdir, patients(dataset), testpat, trainpat)

    for p in classperc:
        # reading training data from given run of boruta analysis
        X_train, y_train = read_typedata(chrlist, outdir, p, borutarun, 'train')

        # establishing testing data based on given test set(s) or test subset of patients
        if not testset:
            if testsize == 0:
                raise exceptions.NoParameterError(
                    'testset', 'Test size was not given - define what set should be used as a testset.')
            X_test, y_test = read_typedata(chrlist, outdir, p, borutarun, 'test')

        else:
            selected_snps = read_selected_snps(chrlist, outdir, p, borutarun)
            X_test, y_test = build_testdata(chrlist, selected_snps, testset)

        # running classification and saving scores to class_scores file
        if X_train.shape[1] > 0:
            score_train, score_test = classify(X_train, y_train, X_test, y_test)
            scores_file.write('%d\t%d\t%.3f\t%.3f\n' % (p, X_train.shape[1], score_train, score_test))
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
    trainstr = ' + '.join(dataset.keys())
    if not testset:
        teststr = '%.2f*(%s)' % (testsize, trainstr)
    else:
        teststr = ' + '.join(testset.keys())
    run_file = open('%sclass_runs.txt' % outdir, 'a')
    'run\ttest_set\ttest_pat\ttrain_run\ttrain_set\ttrain_pat\tperc\tchromosomes\n'
    run_file.write('%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\n' % (classrun, teststr, len(testpat), borutarun, trainstr,
                                                         len(trainpat), ','.join(classperc), funcs.make_chrstr(chrlist)))
    run_file.close()
