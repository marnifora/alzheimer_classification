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
sort=0
overwrite=0

indir="./"
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
        -sort )                 sort=1
                                ;;
        -plinkdir )             shift
                                plinkdir=$1
                                ;;
        -overwrite )            overwrite=1
                                ;;
        *)                      usage
                                exit 1
    esac
    shift
done

if [[ ! -v outdir ]]; then
    outdir=${indir}
fi

# convert binary to text format
if [[ ${binary} -eq 1 ]]; then
    ${plinkdir}plink-1.9/plink --bfile ${indir}${plink} --recode --not-chr 0 --out ${indir}${plink}
    echo "Convertion of binary plink files into text format done!"
fi

# sort SNPs in dbsnp and in map file
if [[ ${sort} -eq 1 ]]; then
    echo "Sorting has just began!"
    export LC_ALL=C
    if [[ ! -f ${dbsnp/./_ascii.} ]]; then
        echo "Sorting dbsnp!"
        sort -k1 ${dbsnp} >${dbsnp/./_ascii.}
    fi
    echo "Sorting map file!"
    sort -k2 ${indir}${plink}".map" >${indir}${plink}"_ascii.map"
    echo "File(s) were sorted successfully!"
fi

# compare SNPs with reference, make pid_chr.txt and genome_stats.txt, filter plink files
if [[ ${step1} -eq 1 ]]; then
    echo "Step one has just began!"
    python3 plink_step_one.py -indir ${indir} -outdir ${outdir} -plink ${plink} -dbsnp ${dbsnp}
    echo "Step one done!"
    ${plinkdir}plink-1.9/plink --file ${indir}${plink} --exclude ${outdir}"missing_snps_ref.txt" --recode --out ${indir}${plink}"_filtered"
    echo "SNPs were excluded, filtered plink files have been made."
fi

# make files with matrices and snps lists
if [[ ${step2} -eq 1 ]]; then
    echo "Step two has just began!"
    if [[ ${overwrite} -eq 1 ]]; then
        python3 plink_step_two.py -indir ${indir} -outdir ${outdir} -plink ${plink}"_filtered" -overwrite
    else
        python3 plink_step_two.py -indir ${indir} -outdir ${outdir} -plink ${plink}"_filtered"
    fi
    echo "Step two done!"
fi