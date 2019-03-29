import sys
sys.path.insert(0, '../')
import exceptions
import corporate_funcs as funcs


def best_snp(directory, ch, borutarun, perc):
    with open('%sboruta/bestsnps_chr%d_%d_%d.txt' % (directory, ch, perc, borutarun), 'r') as file:
        for _ in range(2):
            file.readline()
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
    return snp, iterator


def first_intersection(dataset, ch, borutarun=None, perc=None):

    shared = {}
    ref = {}
    names = list(dataset.keys())[:2]
    iter_snps = [snp_list(dataset[names[0]], ch), snp_list(dataset[names[1]], ch)]
    if borutarun:
        iter_best = [best_snp(dataset[names[0]], ch, borutarun[names[0]], perc), best_snp(dataset[names[1]], ch,
                                                                                          borutarun[names[1]], perc)]

    order = [0, 0]

    try:
        snps = list(map(next, iter_snps))

        if borutarun:
             bests = list(map(next, iter_best))
             snps, iter_snps = list(map(list, zip(*[wind_to_best(iterator, snp, best) for iterator, snp, best in
                                                       zip(iter_snps, snps, bests)])))

        while len(snps) == 2:

            if int(snps[0][0]) < int(snps[1][0]):

                snps[0] = next(iter_snps[0])
                if borutarun:
                    bests[0] = next(iter_best[0])
                    snps[0], iter_snps[0] = wind_to_best(iter_snps[0], snps[0], bests[0])
                order[0] += 1

            elif int(snps[0][0]) == int(snps[1][0]):

                if snps[0][1] == snps[1][1]:
                    shared[int(snps[0][0])] = order
                    ref[int(snps[0][0])] = snps[0][1]
                    snps = list(map(next, iter_snps))
                    if borutarun:
                        bests = list(map(next, iter_best))
                        snps, iter_snps = list(map(list, zip(*[wind_to_best(iterator, snp, best) for iterator, snp, best
                                                               in zip(iter_snps, snps, bests)])))
                    order = [x + 1 for x in order]
                else:
                    raise exceptions.SNPReferenceError(ch, int(snps[0][0]), snps[0][1], snps[1][1])
            else:
                snps[1] = next(iter_snps[1])
                if borutarun:
                    bests[1] = next(iter_best[1])
                    snps[1], iter_snps[1] = wind_to_best(iter_snps[1], snps[1], bests[1])
                order[1] += 1
    except StopIteration:
        pass

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


def map_rows_to_locs(dataset, ch, perc, borutarun, outfile):

    directory = next(iter(dataset.values()))
    snps_loc = snp_list(directory, ch)
    snp = next(snps_loc)
    print('Mapping rows to locations for chromosome %d' % ch)
    with open('%sboruta/bestsnps_chr%d_%d_%d.txt' % (directory, ch, perc, borutarun), 'r') as file:
        for _ in range(2):
            file.readline()
        for i, line in enumerate(file):
            if snp[-1] == int(line.strip()):
                outfile.write('chr%d\t%d\t%d\n' % (ch, int(snp[0]), int(snp[0])+1))
            try:
                snp = next(snps_loc)
            except IndexError:
                break
    return 0
