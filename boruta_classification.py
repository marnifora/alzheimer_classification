import ast
import boruta
from collections import Counter
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


def pooling(chrlist, classperc, dataset, outdir, pat, perc, r, run, subset, subsetrun, testpat):
    """
    Running function one_process for every chr on chrlist.
    Getting list of selected SNPs for every chromosome (returned by every process).
    Returning selected SNPs for every perc in class_perc as a dictionary named selected_snps (keys - chromosomes,
    values - dict (keys - perc, values -numbers of SNPs)).
    Writing number of selected SNPs for every chromosome for every perc to file named all_snps<run>.txt.
    """
    procs = []
    q = multiprocessing.Queue()
    qytrain = multiprocessing.Queue()
    if testpat:
        qytest = multiprocessing.Queue()
    else:
        qytest = None

    for ch in chrlist:

        p = multiprocessing.Process(target=one_process,
                                    args=(ch, classperc, dataset, outdir, pat, perc, q, qytest, qytrain, r, run,
                                          subset, subsetrun, testpat))
        procs.append(p)
        p.start()

    for p in procs:
        p.join()

    selected_snps = {ch: None for ch in chrlist}
    all_snps = {p: 0 for p in perc}
    while q.qsize():
        qq = q.get()
        selected_snps[qq[0]] = qq[2]
        for j, p in enumerate(perc):
            all_snps[p] += qq[1][p]

    if qytest is None:
        params = [[qytrain, 'train']]
    else:
        params = [[qytrain, 'train'], [qytest, 'test']]
    for qy, type in params:
        yt = qy.get()
        while qy.qsize():
            vec = qy.get()
            if np.array_equal(vec, yt):
                pass
            else:
                raise exceptions.OtherError('Y %s matrix is different for different chromosomes!' % type)
        np.save('%sy_%s_%d.npy' % (outdir, type, run), yt)

    a = open('%sall_snps%d.txt' % (outdir, run), 'w')
    for p in perc:
        a.write('%d\t%d\n' % (p, all_snps[p]))
    a.close()

    return selected_snps


def one_process(ch, classperc, dataset, outdir, pat, perc, q, qytest, qytrain, r, run, subset, subsetrun, testpat):
    """
    Loading data to matrices - function load_data.
    Selecting best SNPs subsets for every perc - function best_snps.
    Writing best SNPs for every chromosome into a file.
    Adding to multiprocessing-queue an element: [number of chromosome,
        dictionary - key: perc, value: selected SNPs by Boruta with perc.
    """

    print('Analysis for chromosome %d started\n' % ch)

    X, y, snp, Xtest, ytest = load_data(ch, dataset, pat, subset, subsetrun, testpat)

    print('matrices X and y for chromosome %d have been loaded\n' % ch)

    snps = best_snps(perc, r, snp, X, y)

    print('best SNPs for chromosome %d have been selected by Boruta\n' % ch)

    qytrain.put(y)
    if testpat:
        qytest.put(ytest)
    for p in classperc:
        np.save('%sX_train_chr%d_%d_%d.npy' % (outdir, ch, p, run), X[:, snps[p]])
        print('X train matrix for chr %d for perc %d was saved to file' % (ch, p))
        if testpat:
            np.save('%sX_test_chr%d_%d_%d.npy' % (outdir, ch, p, run), Xtest[:, snps[p]])

    ll = {p: 0 for p in perc}
    for p in perc:

        ll[p] = len(snps[p])

        lista = open('%sbestsnps_chr%d_%d_%d.txt' % (outdir, ch, p, run), 'w')
        lista.write('%d\n\n' % len(snps[p]))
        for el in snps[p]:
            lista.write('%d\n' % el)
        lista.close()

    print('process for chr %d finished\n' % ch)

    q.put([ch, ll, snps])


def load_data(ch, dataset, pat, subset, subsetrun, testpat):
    """
    Loading data from files into X and y matrices.
    Selection of patients not being in test set and SNPs present in subset-file.
    """

    snplist = {name: [] for name in dataset.keys()}
    if subset is not None:
        for name in dataset.keys():
            if subsetrun is None:
                cc = open('%s%s_snps_chr%d.txt' % (dataset[name], subset, ch), 'r')
            else:
                cc = open('%s%s_snps_chr%d_%d.txt' % (dataset[name], subset, ch, subsetrun), 'r')
            for line in cc:
                snplist[name].append(int(line.split()[0]))
            cc.close()
        snp = len(snplist[name])
    else:
        if len(dataset) > 1:
            raise exceptions.NoParameterError('subset', 'There is more than one given data set, but subset of SNPs ' +
                                                        'is not given.')
        else:
            cc = open('%sgenome_stats.txt' % list(dataset.values())[0], 'r')
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

    if testpat is None:
        testpat = []
    test = len(testpat)

    train_row = 0
    test_row = 0
    done = 0

    X_train = np.zeros(shape=(sum(pat.values()) - test, snp), dtype=np.int8)
    y_train = np.zeros(shape=(sum(pat.values()) - test,), dtype=np.int8)

    X_test = np.zeros(shape=(test, snp), dtype=np.int8)
    y_test = np.zeros(shape=(test,), dtype=np.int8)

    for name in dataset.keys():

        o = open('%sX_chr%d_nodif.csv' % (dataset[name], ch), 'r')
        reader = csv.reader(o, delimiter=',')
        next(reader)  # header

        y = pd.read_csv('%sY_chr.csv' % dataset[name], header=None, index_col=0).values
        y = y.ravel()

        # writing values from file to matrix X and Y

        for i, line in enumerate(reader):
            if (done + i) not in testpat:
                y_train[train_row] = y[i]
                for j, s in enumerate(snplist[name]):
                    X_train[train_row][j] = line[s + 1]
                train_row += 1
            else:
                y_test[test_row] = y[i]
                for j, s in enumerate(snplist[name]):
                    X_test[test_row][j] = line[s + 1]
                test_row += 1

        o.close()
        done += i + 1

    if testpat:
        return X_train, y_train, snp, X_test, y_test
    else:
        return X_train, y_train, snp, None, None


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


def build_testdata(chrlist, p, selected_snps, testset):

    pat = patients(testset)
    for name in testset.keys():

        xx = np.zeros(shape=(pat[name], sum([len(x) for x in [selected_snps[ch][p] for ch in chrlist]])), dtype=np.int8)
        col = 0
        for ch in chrlist:
            o = open('%sX_chr%d_nodif.csv' % (testset[name], ch), 'r')
            reader = csv.reader(o, delimiter=',')
            next(reader)

            snps = selected_snps[ch][p]
            for i, line in enumerate(reader):
                xx[i][col:col+len(snps)] = [line[1:][ii] for ii in snps]

            col += len(snps)
            o.close()

        yy = pd.read_csv('%sY_chr.csv' % testset[name], header=None, index_col=0).values
        try:
            X_test = np.concatenate((X_test, xx), axis=0)
            y_test = np.concatenate((y_test, yy), axis=0)
        except NameError:
            X_test = xx
            y_test = yy

    return X_test, y_test


def classify(X_train, y_train, X_test, y_test):

    rf = RandomForestClassifier(n_estimators=500)
    rf.fit(X_train, y_train)
    return rf.score(X_train, y_train), rf.score(X_test, y_test)


def first_run(fixed, outdir, pat, run, testsize):

    run = funcs.establish_run('boruta', fixed, outdir, run)
    p = sum(pat.values())

    if testsize != 0:
        testpat = random.sample(range(p), int(p*testsize))
        ts = open('%stestpat_%d.txt' % (outdir, run), 'w')
        for el in testpat:
            ts.write('%d\n' % el)
        ts.close()
    else:
        testpat = None

    return run, testpat


def cont_run(chrlist, fixed, outdir, run):

    run_file = open('%sboruta_runs.txt' % outdir, 'r+')
    lines = run_file.readlines()
    towrite = ''
    occur = False
    for line in lines:
        if line.startswith(str(run) + '\t'):
            line = line.strip().split('\t')
            if len(line[3].split('-run')) == 2:
                subset, subsetrun = line[3].strip('-run')
                subsetrun = int(subsetrun)
            else:
                subset = line[3]
                subsetrun = None
            testsize = float(line[4])
            perc = ast.literal_eval(line[5])
            if isinstance(perc, int):
                perc = [perc]
            r = int(line[6])
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
            strin = ''
            for el in line:
                strin += str(el) + '\t'
            strin += '\n'
            line = strin
            occur = True
        towrite += line
    run_file.close()

    if not occur:
        raise exceptions.WrongValueError('-run', str(run),
                                         'You set that it is a continuation, but this run has not been conducted yet.')

    if testsize != 0:
        ts = open('%stestpat_%d.txt' % (outdir, run), 'r')
        testpat = []
        for line in ts:
            testpat.append(int(line.strip()))
        ts.close()
    else:
        testpat = None

    return perc, r, subset, subsetrun, testpat, testsize, towrite


def patients(dataset):

    pat = {name: 0 for name in dataset.keys()}
    for name in dataset.keys():
        g = open('%sgenome_stats.txt' % dataset[name], 'r')
        line = g.readline()
        p = int(line.split()[3])
        for line in g:
            if int(line.split()[3]) != p:
                raise exceptions.OtherError('Error: there is different number of patients for different chromosomes!')
        pat[name] = p
        g.close()
    return pat


perc = [90]
r = 5000
chrlist = [i for i in range(1, 24)]

dataset = {}
testset = {}
testsize = 0
class_only = False
boruta_only = False
borutarun = None
classrun = None
subsetrun = None
fixed = False
snp_subset = None
continuation = False

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
        perc = ast.literal_eval(sys.argv[q + 1])
        if isinstance(perc, int):
            perc = [perc]
        continue

    if sys.argv[q] == '-classperc':
        classperc = ast.literal_eval(sys.argv[q + 1])
        if isinstance(classperc, int):
            classperc = [classperc]
        continue

    if sys.argv[q] == '-subset':
        snp_subset = sys.argv[q + 1]
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

    if sys.argv[q] == '-subsetrun':
        subsetrun = int(sys.argv[q + 1])
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
    outdir = next(iter(dataset.values()))

if 'classperc' not in globals():
    classperc = perc


if not class_only:

    # determination number of patient in given data set
    pat = patients(dataset)

    # determination of some parameters
    if not continuation:
        borutarun, testpat = first_run(fixed, outdir, pat, borutarun, testsize)
    else:
        perc, r, subset, subsetrun, testpat, testsize, towrite = cont_run(chrlist, fixed, outdir, borutarun)

    # running Boruta analysis
    selected_snps = pooling(chrlist, classperc, dataset, outdir, pat, perc, r, borutarun, snp_subset,
                            subsetrun, testpat)

    # saving information about done run to boruta_runs file
    if continuation:
        run_file = open('%sboruta_runs.txt' % outdir, 'w')
        run_file.write(towrite)
    else:
        run_file = open('%sboruta_runs.txt' % outdir, 'a')
        if subsetrun is not None:
            subsetstr = '%s-run%d' % (snp_subset, subsetrun)
        else:
            subsetstr = snp_subset
        run_file.write('%d\t%s\t%d\t%s\t%.1f\t%s\t%d\t%s\n' % (borutarun, ' + '.join(dataset.keys()), sum(pat.values()),
                                                               subsetstr, testsize, ','.join(list(map(str, perc))), r,
                                                               funcs.make_chrstr(chrlist)))
    run_file.close()


if not boruta_only:

    # determination of number of class run
    classrun = funcs.establish_run('class', fixed, outdir, classrun)
    scores_file = open('%sclass_scores_%d.txt' % (outdir, classrun), 'w', 1)
    num_snps_perc = []
    selected_snps = {ch: {} for ch in chrlist}
    for i, p in enumerate(classperc):
        # reading training data from given run of boruta analysis
        X_train, y_train = read_typedata(chrlist, outdir, p, borutarun, 'train')

        # establishing testing data based on given test set(s) or test subset of patients
        if not testset:
            if testsize == 0:
                raise exceptions.NoParameterError(
                    'testset', 'Test size was not given - define what set should be used as a testset.')
            X_test, y_test = read_typedata(chrlist, outdir, p, borutarun, 'test')

        else:
            if class_only:
                for ch in chrlist:
                    o = open('%sbestsnps_chr%d_%d_%d.txt' % (outdir, ch, p, borutarun), 'r')
                    for i in range(2):  # header
                        o.readline()
                    for line in o:
                        selected_snps[ch][p] = selected_snps[ch].setdefault(p, []).append(int(line.strip()))
                    o.close()
            X_test, y_test = build_testdata(chrlist, p, selected_snps, testset)

        # writing heading to class_scores file
        if i == 0:
            trainpat = X_train.shape[0]
            testpat = y_test.shape[0]
            trainstr = ' + '.join(dataset.keys())
            if not testset:
                teststr = '%.1f*(%s)' % (testsize, trainstr)
            else:
                teststr = ' + '.join(testset.keys())
            scores_file.write('Random forest classification\n' +
                              'TRAINING DATA:\nData set =  %s\nPatients = %d\ntrain run = %d\n'
                              % (trainstr, trainpat, borutarun) +
                              'TESTING DATA:\nData set = %s\nPatients = %d\n'
                              % (teststr, testpat))
            scores_file.write('RESULT of ANALYSIS:\nperc\tSNPs\ttrain_score\ttest_score\n')

        num_snps = X_train.shape[1]
        num_snps_perc.append(num_snps)

        # running classification and saving scores to class_scores file
        if num_snps > 0:
            score_train, score_test = classify(X_train, y_train, X_test, y_test)
            print('Classification done!')
            scores_file.write('%d\t%d\t%.3f\t%.3f\n' % (p, num_snps, score_train, score_test))
            # saving matrices (on which was based the classification) to file
            for name in ['X_train', 'y_train', 'X_test', 'y_test']:
                np.save('%s%s_genome_%d_%d.npy' % (outdir, name, p, classrun), eval(name))
        else:
            print('No SNPs were chosen for perc %d' % p)
            scores_file.write('%d\t%d\t-\t-\n' % (p, num_snps))

        print('Scores for perc %d saved to file' % p)

    scores_file.close()

    # writing information about class run to class_run file
    run_file = open('%sclass_runs.txt' % outdir, 'a')
    'run\ttest_set\ttest_pat\ttrain_run\ttrain_set\ttrain_pat\tperc\tSNPs\tchromosomes\n'
    run_file.write('%d\t%s\t%d\t%d\t%s\t%d\t%s\t%s\t%s\n' % (classrun, teststr, testpat, borutarun, trainstr, trainpat,
                                                             classperc, num_snps_perc, funcs.make_chrstr(chrlist)))
    run_file.close()

