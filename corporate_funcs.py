import exceptions
import numpy as np
import csv
import os


def establish_run(analysistype, fixed, outdir, run):

    directory = outdir
    name = '-'
    if not os.path.isdir(directory):
        os.mkdir(directory)
    try:
        file = '%s%s_runs.txt' % (directory, analysistype)
        run_file = open(file, 'r+')
        if run is None:
            run = 0
            runchanged = False
        else:
            runchanged = True
        lines = run_file.readline()  # header
        rr = []
        rewrite = False
        for line in run_file:
            try:
                val = int(line.split()[0])
            except ValueError:
                continue
            rr.append(val)
            if val != run:
                lines += line
            else:
                rewrite = True
        if rr:
            d = [i for i in range(1, max(rr)+2)]
            for el in rr:
                d.remove(el)
            if not runchanged:
                run = min(d)
                print('%s run number for %s has been established! Run = %d' % (analysistype, name, run))
            elif rewrite:
                if not fixed:
                    raise exceptions.WrongValueError('-run', run,
                                                     "Run number %d has already been conducted (%s analysis for %s)! "
                                                     % (run, analysistype, name) +
                                                     "If you want to overwrite it, please add '-fixed' atribute.")
                else:
                    run_file.seek(0)
                    run_file.write(lines)
                    run_file.truncate()
        else:
            if not runchanged:
                run = 1
                print('%s run number for %s has been established! Run = %d' % (analysistype, name, run))

    except FileNotFoundError:
        run = 1
        run_file = open('%s%s_runs.txt' % (directory, analysistype), 'w')
        if analysistype == 'boruta':
            run_file.write('run\tdata_set\tpatients\tpat_subset\tpat_runs\tSNPs_subset\tSNPs_runs\ttest_size\tperc\t' +
                           'window_size\tchromosomes\toutdir\n')
        elif analysistype == 'class':
            run_file.write('run\ttest_set\ttest_pat\ttrain_run\ttrain_set\ttrain_pat\tperc\tchromosomes\toutdir\tclassifier\n')
        elif analysistype == 'shared':
            run_file.write('run\thome_set\tcompared_set(s)\tchromosomes\tnumber_of_shared_SNPs\n')
        elif analysistype == 'crossed':
            run_file.write('run\thome_set\tcompared_set(s)\tchromosomes\tnumber_of_shared_SNPs\tboruta_run(s)\tperc\n')
        elif analysistype == 'similar':
            run_file.write('run\tdata_set(s)\tlower_thresh\tupper_thresh\tsimilar_pat\tall_pat\n')
        elif analysistype == 'deficient':
            run_file.write(
                'run\tdata_set(s)\tboruta_run(s)\tperc\tthreshold\tnumber_of_deficient_SNPs\tpercent_of_all\n')
        elif analysistype == 'frombed':
            run_file.write('run\tinfile\ttemplate_dataset\tnumber_of_SNPs\tchromosomes\n')
        else:
            raise exceptions.OtherError('First line for %s run file is not defined!' % analysistype)
        print('%s run file for %s has been made! Run number has been established! Run = %d' % (analysistype, name, run))

    run_file.write('%d\tin progress\n' % run)
    run_file.close()

    return run


def correct_boruta_runs_file(file):
    o = open(file, 'r+')
    line = o.readline().split()
    if len(line) <= 8:
        towrite = '\t'.join(line[:3] + ['pat_subset', 'pat_runs', 'SNPs_subset', 'SNPs_runs'] + line[4:]) + '\n'
        for line in o:
            line = line.split()
            if len(line[3].split('-run')) == 2:
                subset, subsetrun = line[3].strip().split('-run')
            else:
                subset, subsetrun = 'None', '-'
            towrite += '\t'.join(line[:3] + ['None', '-', subset, subsetrun] + line[4:]) + '\n'
        o.seek(0)
        o.write(towrite)
        o.truncate()
    o.close()


def runs_file_add(analysistype, directory, run, toadd):
    run_file = open('%s%s_runs.txt' % (directory, analysistype), 'r+')
    lines = run_file.readlines()
    towrite = ''
    for line in lines:
        if not line.startswith('%d\tin progress' % run):
            towrite += line
        else:
            towrite += toadd
    run_file.seek(0)
    run_file.write(towrite)
    run_file.truncate()
    run_file.close()


def runs_file_rewrite(analysistype, directory, towrite):
    run_file = open('%s%s_runs.txt' % (directory, analysistype), 'w')
    run_file.write(towrite)
    run_file.close()


def make_chrstr(chrlist):

    cl = chrlist.copy()
    cl.append(0)
    chrstr = ''
    first = cl[0]
    och = first
    for ch in cl:
        if ch == och:
            och += 1
        elif first != och-1:
            if len(chrstr) != 0:
                chrstr += ', '
            chrstr += '%d-%d' % (first, och-1)
            first = ch
            och = ch+1
        else:
            if len(chrstr) != 0:
                chrstr += ', '
            chrstr += '%d' % first
            first = ch
            och = ch+1

    return chrstr


def read_chrstr(chrstr):

    chrstr = chrstr.strip('[]')
    c = chrstr.split(',')
    chrlist = []
    for el in c:
        el = el.split('-')
        if len(el) == 1:
            chrlist.append(int(el[0]))
        else:
            chrlist += [i for i in range(int(el[0]), int(el[1])+1)]
    chrlist.sort()

    return chrlist


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

    return read_Xs(ch, dataset, snp, snplist, testpat, trainpat)


def read_Xs(ch, dataset, snp, snplist, testpat, trainpat):

    train_row = 0
    test_row = 0
    done = 0

    X_train = np.zeros(shape=(len(trainpat), snp), dtype=np.int8)

    if testpat:
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


def patients_diagnoses(dataset, patients):
    case, control = [], []
    done = 0
    for name in dataset.keys():
        with open('%smatrices/Y_chr.csv' % dataset[name], 'r') as file:
            for i, line in enumerate(file):
                if i in patients:
                    line = list(map(int, line.split(',')))
                    if line[1] == 0:
                        control.append(done + line[0])
                    else:
                        case.append(done + line[0])
        done += i
    return case, control
