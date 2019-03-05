import exceptions


def establish_run(analysistype, fixed, outdir, run):

    try:
        file = '%s%s_runs.txt' % (outdir, analysistype)
        if analysistype == 'boruta':
            correct_boruta_runs_file(file)
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
                print('%s run number has been established! Run = %d' % (analysistype, run))
            elif rewrite:
                if not fixed:
                    raise exceptions.WrongValueError('-run', run,
                                                     "Run number %d has already been conducted (%s analysis)! "
                                                     % (run, analysistype) +
                                                     "If you want to overwrite it, please add '-fixed' atribute.")
                else:
                    run_file.seek(0)
                    run_file.write(lines)
                    run_file.truncate()
        else:
            if not runchanged:
                run = 1
                print('%s run number has been established! Run = %d' % (analysistype, run))

    except FileNotFoundError:
        run = 1
        run_file = open('%s%s_runs.txt' % (outdir, analysistype), 'w')
        if analysistype == 'boruta':
            run_file.write('run\tdata_set\tpatients\tpat_subset\tpat_runs\tSNPs_subset\tSNPs_runs\ttest_size\tperc\t' +
                           'window_size\tchromosomes\n')
        elif analysistype == 'class':
            run_file.write('run\ttest_set\ttest_pat\ttrain_run\ttrain_set\ttrain_pat\tperc\tSNPs\tchromosomes\n')
        elif analysistype == 'shared':
            run_file.write('run\thome_set\tcompared_set(s)\tchromosomes\tnumber_of_shared_SNPs\n')
        elif analysistype == 'similar':
            run_file.write('run\tdata_set(s)\tlower_thresh\tupper_thresh\tsimilar_pat\tall_pat\n')
        else:
            raise exceptions.OtherError('First line for %s run file is not defined!' % analysistype)
        print('%s run file has been made! Run number has been established! Run = %d' % (analysistype, run))

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
                subset, subsetrun = line[3].strip('-run')
            else:
                subset, subsetrun = 'None', '-'
            towrite += '\t'.join(line[:3] + ['None', '-', subset, subsetrun] + line[4:]) + '\n'
        o.seek(0)
        o.write(towrite)
        o.truncate()
    o.close()


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
