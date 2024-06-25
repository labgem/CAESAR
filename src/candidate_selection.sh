#!/bin/bash

echo $(date)
SECONDS=0

script_path=$1
shift 1

outdir="."


Help(){
    echo "Script used only with a slurm sbatch command"
    echo
    echo "usage: candidate_selection.sh <path>/candidate_selection.py [-h] [-o] -C -f [-d] -c [-g] [-n] [-v] -s [-u]"
    echo
    echo "options:"
    echo "  -h    show this help message"
    echo "  -o    output directory [default: './']"
    echo "  -C    the configuration file in yaml format"
    echo "  -f    fasta file containing all candidates sequences"
    echo "  -d    fltered blast output"
    echo "  -c    clusters tsv file"
    echo "  -g    target GC percentage to decide between candidates [default: 50.0]"
    echo "  -n    maximum number of candidate per cluster [default: 1]"
    echo "  -v    selects a fraction of the cluster rather than a specific number"
    echo "        e.g: 10 (10%)"
    echo "  -s    sources file indicating the sources database of each sequences"
    echo "  -u    file containing a list of identifiers, the groups in which they"
    echo "        are found will be excluded"
}

export -f Help

while getopts ":o:C:f:d:c:g:n:v:s:u:" option
do
    case $option in
        o)
            outdir=$OPTARG
            ;;
        C)
            config=$OPTARG
            ;;
        f)
            fasta=$OPTARG
            ;;
        d)
            data=$OPTARG
            ;;
        c)
            clusters=$OPTARG
            ;;
        g)
            gc=$OPTARG
            ;;
        n)
            nb=$OPTARG
            ;;
        v)
            cov=$OPTARG
            ;;
        s)
            sources=$OPTARG
            ;;
        u)
            update=$OPTARG
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

command=$(echo $script_path -o $outdir -c $config -f $fasta --clusters $clusters --sources $sources)

if [ ! -z $data ]
then
    command=$(echo $command -d $data)
fi

if [ ! -z $gc ]
then
    command=$(echo $command --gc $gc)
fi

if [ ! -z $nb ]
then
    command=$(echo $command -n $nb)
fi

if [ ! -z $cov ]
then
    command=$(echo $command --cov-per-cluster $cov)
fi

if [ ! -z $update ]
then
    command=$(echo $command -u $update)
fi

python $command
