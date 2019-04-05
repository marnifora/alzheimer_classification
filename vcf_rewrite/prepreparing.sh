#!/bin/bash

. job_pool.sh

: '
See readme.txt for input, output and possible options.
'

function pooling {

# tar 1
# gz 2
# filtered_snp 3
# do stats 4
# do matrix 5
# chromosome number 6
# vcf file name 7
# filtered vcf file name 8
# directory of files (indir) 9
# directory of matrices (outdir) 10
# dir to gatk 11
# genome reference 12

        vcf=${9}${7}
        vcfsnps=${9}${8}"_chr"${6}"_SNPs.vcf"

        if [[ $1 -eq 1 ]]; then

                tend=".tar"
                tarfile=${vcf}${tend}
                tar -xvf ${tarfile} -C $9
                echo "tar file for chr=$6 decompressed!"
        fi

        if [[ $2 -eq 1 ]]; then

                gend=".gz"
                gzfile=${vcf}${gend}
                gunzip -f ${gzfile} > ${vcf}
                echo "gz file for chr=$6 decompressed!"
        fi
        

        if [[ $3 -eq 1 ]]; then
                echo "running SelectVariants for chr=$4"
                ${11}gatk SelectVariants -R ${12} -V ${vcf} -O ${vcfsnps} -select-type-to-include SNP
                echo "$6 SNPs.vcf done!"
        fi

        if [[ $4 -eq 1 ]]; then
                echo "stats for chr=$6 started!"
                python3 vcf_stats.py -input ${vcfsnps}
                echo "stats file for chr=$6 done!"
        fi

        if [[ $5 -eq 1 ]]; then
                echo "running vcf_to_matrix for chr $6"
                python3 vcf_to_matrix.py -chr $6 -input ${vcfsnps} -outdir ${10} >> ${10}'genome_stats.txt'
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

#name of input vcf files
vcfstart='default_chr'
vcfend='.vcf'

# name of database
base="dataset1"

# directory to the input files
gatk_dir="./gatk-4.0.11.0/"
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
        -indir | -files_dir)    shift
                                files_dir=$1
                                ;;
        -outdir | -matrices_dir) shift
                                matrices_dir=$1
                                ;;
        -gatk_dir )             shift
                                gatk_dir=$1
                                ;;
        -reference )            shift
                                reference=$1
                                ;;
        -vcf | -vcffile )       shift
                                s=(${1//\{chr\}/ })
                                vcfstart=${s[0]}
                                vcfend=${s[1]}
                                ;;
        *)                      usage
                                exit 1
        esac
        shift
done

if [[ ! -v ${files_dir} ]]; then
    files_dir=${dir}"files/"
fi

if [[ ! -d ${files_dir} ]]; then
    mkdir ${files_dir}
fi

if [[ ! -v ${matrices_dir} ]]; then
    matrices_dir=${dir}"matrices/"
fi

if [[ ! -d ${matrices_dir} ]]; then
    mkdir ${matrices_dir}
fi

echo "Parameters are established!"

# prefix and sufix of name of vcf files containing WGS data, between them there should be only number of chromosome
# (it will be added automatically)

# name of downloaded data for rosmap: NIA_JG_1898_samples_GRM_WGS_b37_JointAnalysis01_2017-12-08_{chr}.recalibrated_variants.vcf
# name of dowloaded data for adni: ADNI.808_indiv.minGQ_21.pass.ADNI_ID.chr{chr}.vcf

# name of fasta file with reference genome based on which WGS files were made
# reference for rosmap: human_g1k_v37.fasta
# reference for adni: Homo_sapiens_assembly19.fasta

if [[ $all -eq 1 ]]; then
        job_pool_init $((to-from+1)) 0

        for (( ch = $from; ch <= $to; ch++ )); do

                if [[ $ch -eq 23 ]] && [[ $x -eq 1 ]]; then

                        vcfname=${vcfstart}X${vcfend}

                elif [[ $ch -eq 24 ]]; then

                        vcfname=${vcfstart}Y${vcfend}

                else
                        vcfname=${vcfstart}${ch}${vcfend}
                fi

                job_pool_run pooling ${tar} ${gz} ${snp} ${stats} ${matrix} ${ch} ${vcfname} ${base} ${files_dir}
                ${matrices_dir} ${gatk_dir} ${ref}

        done

        job_pool_shutdown
        echo "job_pool_nerrors: ${job_pool_nerrors}"
else
        vcfname=${vcfstart}${ch}${vcfend}
        pooling ${tar} ${gz} ${snp} ${stats} ${matrix} ${ch} ${vcfname} ${base} ${files_dir} ${matrices_dir} ${gatk_dir} ${ref}
fi

if [[ $matrix -eq 1 ]]; then

        sort -k1 -n -o ${outdir}'genome_stats.txt' ${outdir}'genome_stats.txt'

fi
