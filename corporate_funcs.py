import exceptions


def establish_run(analysistype, fixed, outdir, run):

    try:
        run_file = open('%s%s_runs.txt' % (outdir, analysistype), 'r+')
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
            run_file.write('run\tdata_set\tpatients\tsnps_subset\ttest_size\tperc\twindow_size\tchromosomes\n')
        elif analysistype == 'class':
            run_file.write('run\ttest_set\ttest_pat\ttrain_run\ttrain_set\ttrain_pat\tperc\tSNPs\tchromosomes\n')
        elif analysistype == 'shared':
            run_file.write('run\thome_set\tcompared_set(s)\tchromosomes\tnumber_of_shared_SNPs\n')
        else:
            raise exceptions.OtherError('First line for %s run file is not defined!' % analysistype)
        print('%s run file has been made! Run number has been established! Run = %d' % (analysistype, run))

    run_file.close()

    return run


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
