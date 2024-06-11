#!/bin/bash


fasta=""
outdir="."
id=80.0
cov=80.0
threads=1
Memory="4G"

Help(){
    echo "Build a diamond database and runs a clustering"
    echo
    echo "usage: clustering.sh [-h] -f [-o] [-t] [-i] [-c] [-m]"
    echo
    echo "options:"
    echo "  -h    show this help message"
    echo "  -f    multi fasta file"
    echo "  -o    output directory ['./]"
    echo "  -t    number of threads [defaults: 1]"
    echo "  -i    identity cutoff [default: 80.0]"
    echo "  -c    minimum coverage of cluster member [default: 80.0]"
    echo "  -m    memory limit for the diamond process [default: '4G']"
}

export -f Help

while getopts "hf:o:i:c:t:m:" options
do
    case $options in
        h)
            Help
            exit 0
            ;;
        f)
            fasta=$OPTARG
            ;;
        o)
            outdir=$OPTARG
            ;;
        i)
            id=$OPTARG
            ;;
        c)
            cov=$OPTARG
            ;;
        t)
            threads=$OPTARG
            ;;
        m)
            Memory=$OPTARG
            ;;
        :)
            exit 1
            ;;
        \?)
            exit 1
            ;;
        *)
            exit 1
            ;;
    esac
done

if [ ! -f $outdir ] && [ $outdir != "." ]
then
    mkdir $outdir
fi

if [ ! -f $fasta ] || [ -z $fasta ]
then
    echo "Error: argument -f is required and must be an existing file"
    echo
    Help
    exit 1
fi

echo $(date)
SECONDS=0

diamond makedb --threads $threads --in $fasta -d $outdir/candidate_db.dmnd
diamond cluster --threads $threads -d $outdir/candidate_db -o $outdir/clusters.tsv --approx-id $id --member-cover $cov -M $Memory

ELAPSED=$SECONDS
printf 'Elapsed time: %02dh:%02dm:%02ds\n'  $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
