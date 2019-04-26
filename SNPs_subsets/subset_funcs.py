import sys
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs


def best_snp(borutadir, directory, ch, borutarun, perc, snpsubset, snprun):

    if snpsubset is not None:
        subfile = open('%s%s/%s_snps_chr%d_%d.txt' % (directory, snpsubset, snpsubset, ch, snprun), 'r')

    with open('%sboruta/bestsnps_chr%d_%d_%d.txt' % (borutadir, ch, perc, borutarun), 'r') as file:
        for _ in range(2):
            file.readline()
        if snpsubset is not None:
            done = 0
            for line in file:
                for _ in range(int(line.strip())-1-done):
                    subfile.readline()
                done = int(line.strip())
                yield int(subfile.readline().strip())
        else:
            for line in file:
                yield int(line.strip())


def shared_snp(directory, ch, sharedrun):
    with open('%sshared/shared_snps_chr%d_%d.txt' % (directory, ch, sharedrun), 'r') as file:
        for line in file:
            yield int(line.strip())


def crossed_snp(directory, ch, crossedrun):
    with open('%scrossed/crossed_snps_chr%d_%d.txt' % (directory, ch, crossedrun), 'r') as file:
        for line in file:
            yield int(line.strip())


def snp_list(directory, ch):
    with open('%smatrices/snps_chr%d.txt' % (directory, ch), 'r') as file:
        for i, line in enumerate(file):
            yield line.strip().split() + [i]


def wind_to_best(iterator, snp, best):

    while snp[-1] < best:
        snp = next(iterator)
    if snp[-1] > best:
        raise exceptions.OtherError('There is no SNP with number %d' % best)
    return iterator, snp


def first_intersection(dataset, ch, borutarun=None, perc=None):

    shared = {}
    ref = {}
    dirs = iter(dataset.values())
    iter_snps = [snp_list(next(dirs), ch), snp_list(next(dirs), ch)]
    if borutarun:
        iter_best = []
        for name, directory in dataset.items():
            perc, snpsubset, snprun = check_borutarun(directory, borutarun[name], perc)
            iter_best.append(best_snp(directory, ch, borutarun[name], perc, snpsubset, snprun))
    snps = list(map(next, iter_snps))

    if borutarun:
        bests = list(map(next, iter_best))
        iter_snps, snps = list(map(list, zip(*[wind_to_best(iterator, snp, best) for iterator, snp, best in
                                               zip(iter_snps, snps, bests)])))

    while len(snps) == 2:

        if int(snps[0][0]) < int(snps[1][0]):

            try:
                snps[0] = next(iter_snps[0])
                if borutarun:
                    bests[0] = next(iter_best[0])
                    iter_snps[0], snps[0] = wind_to_best(iter_snps[0], snps[0], bests[0])
            except StopIteration:
                break

        elif int(snps[0][0]) == int(snps[1][0]):

            if snps[0][1] == snps[1][1]:
                shared[int(snps[0][0])] = [s[-1] for s in snps]
                ref[int(snps[0][0])] = snps[0][1]
                snps = list(map(next, iter_snps))
                if borutarun:
                    bests = list(map(next, iter_best))
                    if not bests:
                        break
                    iter_snps, snps = list(map(list, zip(*[wind_to_best(iterator, snp, best)
                                                           for iterator, snp, best in zip(iter_snps, snps, bests)])))
            else:
                raise exceptions.SNPReferenceError(ch, int(snps[0][0]), snps[0][1], snps[1][1])
        else:
            try:
                snps[1] = next(iter_snps[1])
                if borutarun:
                    bests[1] = next(iter_best[1])
                    iter_snps[1], snps[1] = wind_to_best(iter_snps[1], snps[1], bests[1])
            except StopIteration:
                break
    return shared, ref


def next_intersection(set, shared, ref, ch):

    snps = iter(shared.keys())
    snp = next(snps)
    for order, line in enumerate(set):
        line = line.split()
        if int(line[0]) == snp:
            if line[1] == ref[snp]:
                shared[snp].append(order)
            else:
                raise exceptions.SNPReferenceError(ch, snp, ref[snp], line[1])
        elif int(line[0]) > snp:
            del (shared[snp])
            try:
                snp = next(snps)
            except StopIteration:
                break

    return shared


def check_borutarun(directory, run, perc):
    with open('%sboruta/boruta_runs.txt' % directory, 'r') as file:
        for line in file:
            if line.startswith(str(run) + '\t'):
                line = line.split()
                snpsubset = line[5]
                if snpsubset == 'None':
                    snpsubset = None
                    snprun = None
                else:
                    snprun = int(line[6].split('+')[0])
                perc_list = line[8].split(',')
                if len(perc_list) > 1 and perc is None:
                    raise exceptions.NoParameterError('perc',
                                                      'There is more than one perc value for given boruta run.')
                elif perc is None:
                    perc = int(perc[0])
                break
    if 'snpsubset' not in locals():
        raise exceptions.WrongValueError('run', run, 'Run number %d was not conducted' % run)
    return perc, snpsubset, snprun
