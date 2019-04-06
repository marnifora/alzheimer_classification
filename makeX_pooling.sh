#!/bin/bash

. job_pool.sh

: '
See readme.txt for input, output and possible options.
'

function make(){

	echo "running makeX.py for chr $1"
	python3 makeX.py -chr $1 -indir $2 -outdir $3
	echo "X matrix for chr $1 done!"

}

ch=0
all=0
from=1
to=24
dir='./'

while [[ "$1" != "" ]]; do
        case $1 in

        -chr | -c )             shift
                                ch=$1
                                ;;
        -all | -a )             all=1
                                ;;
        -from )                 shift
                                all=1
                                from=$1
                                ;;
        -to )                   shift
                                all=1
                                to=$1
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
        *)                      usage
                                exit 1
        esac
    shift
done

if [[ ! -v indir ]]; then
    indir=${dir}"matrices/"
fi

if [[ ! -v outdir ]]; then
    outdir=${dir}"matrices/"
fi

if [[ ${all} -eq 1 ]]; then
        job_pool_init $((to-from+1)) 0

        for (( i=$from; i<=$to; i++ )); do                

                job_pool_run make ${i} ${indir} ${outdir}

        done

        job_pool_shutdown
        echo "job_pool_nerrors: ${job_pool_nerrors}"
else
        make ${ch} ${indir} ${outdir}
fi
