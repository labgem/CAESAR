#!/bin/bash

echo $(date)
SECONDS=0

outdir=""

Help(){
    echo "Script used only with a slurm sbatch command"
    echo
    echo "usage: phylo.sh [-h] [-o] -f [-c] <path>/phylo.py"
    echo
    echo "options:"
    echo "  -h    show this help message"
    echo "  -o    output directory [default: './']"
    echo "  -f    fasta file containing all candidates sequences"
    echo "  -c    clusters tsv file"
}

export -f Help

while getopts ":o:f:c:" option
do
    case $option in
        o)
            outdir=$OPTARG
            ;;
        f)
            fasta=$OPTARG
            ;;
        c)
            cluster=$OPTARG
            ;;
        :)
            Help
            exit 1
            ;;
        \?)
            Help
            exit 1
            ;;
        *)
            Help
            exit 1
            ;;
    esac
done

shift $((OPTIND-1))
script_path=$1

command=$(echo $script_path -f $fasta)

if [ ! -f $outdir ] && [ $outdir != "./" ]
then
    command=$(echo $command -o $outdir)
    mkdir $outdir
fi

if [ ! -z $cluster ]
then
    command=$(echo $command --clusters $cluster)
fi

python $command