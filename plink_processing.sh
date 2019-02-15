#!/bin/bash
set -e

: '
See readme.txt for input, output and possible options.

dbsnp.txt in format:
SNP_ID  REF_BASE
'

binary=0
step1=0
step2=0
nosort=0

indir="./"
outdir="./"
plinkdir="./"

while [[ "$1" != "" ]]; do
        case $1 in

        -plink )                shift
                                plink=$1
                                ;;
        -dbsnp )                shift
                                dbsnp=$1
                                ;;
        -binary )               binary=1
                                ;;
        -step1 )                step1=1
                                ;;
        -step2 )                step2=1
                                ;;
        -indir )                shift
                                indir=$1
                                ;;
        -outdir )               shift
                                outdir=$1
                                ;;
        -nosort )               nosort=1
                                ;;
        -plinkdir )             shift
                                plinkdir=$1
                                ;;
        *)                      usage
                                exit 1
    esac
    shift
done

# convert binary to text format
if [[ ${binary} -eq 1 ]]; then
    ${plinkdir}plink-1.9/plink --bfile ${indir}${plink} --recode --not-chr 0 --out ${indir}${plink}
    echo Convertion of binary plink files into text format done!
fi

if [[ ${nosort} -eq 0 ]]; then
    echo Sorting has just began!
    export LC_ALL=C
    sort -k1 ${indir}${dbsnp} >${indir}${dbsnp/./_ascii.}
    sort -k2 ${indir}${plink}".map" >${indir}${plink}"_ascii.map"
    echo Files were sorted succesfully!
fi

# compare SNPs with reference, make pid_chr.txt and genome_stats.txt, filter plink files
if [[ ${step1} -eq 1 ]]; then
    echo Step one has just began!
    python3 plink_step_one.py -indir ${indir} -outdir ${outdir} -plink ${plink} -dbsnp ${dbsnp}
    echo Step one done!
    ${plinkdir}plink-1.9/plink --file ${indir}${plink} --exclude ${outdir}"missing_snps_ref.txt" --recode --out ${indir}${plink}"_filtered"
    echo SNPs were excluded, filtered plink files have been made.
fi

# make files with matrices and snps lists
if [[ ${step2} -eq 1 ]]; then
    echo Step two has just began!
    python3 plink_step_two.py -indir ${indir} -outdir ${outdir} -plink ${plink}"_filtered"
    echo Step two done!
fi