#!/bin/bash

. job_pool.sh

: '
See readme.txt for input, output and possible options.
'


function vcf {

        if [[ $5 -eq 1 ]]; then

                tend=".tar"
                tarfile=$1${tend}
                tar -xvf ${tarfile} -C ${13}
                echo "tar file for chr=$4 decompressed!"
        fi

        if [[ $6 -eq 1 ]]; then

                gend=".gz"
                gzfile=$1${gend}
                gunzip -f ${gzfile} > $1
                echo "gz file for chr=$4 decompressed!"
        fi
        

        if [[ $7 -eq 1 ]]; then
                echo $8           
                echo "running SelectVariants for chr=$4"
                ${12}gatk SelectVariants -R $3 -V $1 -O $2 -select-type-to-include SNP
                echo "$4 SNPs.vcf done!"
        fi

        if [[ $8 -eq 1 ]]; then
                python3 vcf_stats.py -input $2
                echo "stats file for chr=$4 done!"
        fi

        if [[ $9 -eq 1 ]]; then
                echo "running vcf_to_matrix for chr $4"
                python3 vcf_to_matrix.py -chr $4 -input $2 -outdir ${11} >> ${10}
                sort -k1 -n -o ${10} ${10}
                echo "matrices for chr $4 done!"
        fi
        # after this, run only once each function:
        # > make_pid-diagnoses.py (change rules for establishing diagnoses)
        # > makeY.py
        # > makeX_pooling.sh
        
}


chr=0
all=0
from=1
to=24
x=0

tar=0
gz=0
snp=0
stats=0
matrix=0

#name of input and output files
istart='default_chr'
iend='.vcf'

# name of database
base="dataset1"

# directory to the input files
gatkdir="./gatk-4.0.11.0/"
dir="./"

while [[ "$1" != "" ]]; do
        case $1 in

        -chr | -c )             shift
                                chr=$1
                                ;;
        -all | -a )             all=1
                                ;;
        -x )                    x=1
                                ;;
        -from )                 shift
                                all=1
                                from=$1
                                ;;
        -to )                   shift
                                all=1
                                to=$1
                                ;;
        -tar )                  tar=1
                                ;;
        -gz )                   gz=1
                                ;;
        -snp )                  snp=1
                                ;;
        -stats )                stats=1
                                ;;
        -matrix )               matrix=1
                                ;;
        -base )                 shift
                                base=$1
                                ;;
        -dir )                  shift
                                dir=$1
                                ;;
        -indir )                shift
                                indir=$1
                                ;;
        -outdir )               shift
                                outdir=$1
                                ;;
        -gatkdir )              shift
                                gatkdir=$1
                                ;;
        -reference )            shift
                                reference=$1
                                ;;
        -name )                 shift
                                s=(${1//\{chr\}/ })
                                istart=${s[0]}
                                iend=${s[1]}
                                ;;
        *)                      usage
                                exit 1
    esac
    shift
done

if [[ ! -v ${indir} ]]; then
    indir=${dir}"files/"
fi

if [[ ! -d ${indir} ]]; then
    mkdir ${indir}
fi

if [[ ! -v ${outdir} ]]; then
    outdir=${dir}"matrices/"
fi

if [[ ! -d ${outdir} ]]; then
    mkdir ${outdir}
fi

# prefix and sufix of name of vcf files containing WGS data, between them there should be only number of chromosome
# (it will be added automatically)

# name of downloaded data for rosmap: NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{chr}.recalibrated_variants.vcf
# name of dowloaded data for adni: ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr{chr}.vcf

# name of fasta file with reference genome based on which WGS files were made
# reference for rosmap: human_g1k_v37.fasta
# reference for adni: Homo_sapiens_assembly19.fasta

# name of output file with stats for whole genome
genomestats=${outdir}"genome_stats.txt"

# prefix and sufix of output vcf files
ostart=${base}"_chr"
oend="_SNPs.vcf"


if [[ $all -eq 1 ]]; then
        job_pool_init $((to-from+1)) 0

        for (( i = $from; i <= $to; i++ )); do                

                if [[ $i -eq 23 ]] && [[ $x -eq 1 ]]; then

                        input=${indir}${istart}X${iend}

                elif [[ $i -eq 24 ]]; then

                        input=${indir}${istart}Y${iend}

                else
                        input=${indir}${istart}${i}${iend}
                fi
                        
                output=${indir}${ostart}${i}${oend}

                job_pool_run vcf ${input} ${output} ${reference} ${i} ${tar} ${gz} ${snp} ${stats} ${matrix} ${genomestats} ${outdir} ${gatkdir} ${indir}

        done

        job_pool_shutdown
        echo "job_pool_nerrors: ${job_pool_nerrors}"
else
        input=${indir}${istart}${chr}${iend}
        output=${indir}${ostart}${chr}${oend}
        vcf ${input} ${output} ${reference} ${chr} ${tar} ${gz} ${snp} ${stats} ${matrix} ${genomestats} ${outdir} ${gatkdir} ${indir}
fi
