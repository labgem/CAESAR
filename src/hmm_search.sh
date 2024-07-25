#!/bin/bash

Help(){
    echo "Search one HMM against one or several database with hmmsearch"
    echo
    echo "usage: hmm_search.sh [-h] -q [-o] [-t] [-p] <db.fasta>"
    echo
    echo "options:"
    echo "  -h    show this help message"
    echo "  -q    the query hmm"
    echo "  -o    output directory [default: './']"
    echo "  -t    number of threads [default: 1]"
    echo "  -p    uses gnu parallel to run 2 jobs, each job uses n threads"
    echo "        corresponding to the value of the -t option divide by 2"
}

export -f Help

Run_hmmsearch(){

    echo $1 $2
    if [ $1 == 1 ]
    then
        n=$(($2 / 2))
    else
        n=$2
    fi
    db=$3
    query=$4
    outdir=$5
    
    source=$(echo $db | grep -oE '([A-Za-z0-9_]+)\.fasta' | cut -d . -f 1)
    source=$(echo "${source,,}")

    hmmsearch --cpu $n --noali --domtblout $outdir/hits_${source}.domtbl $query $db
}

export -f Run_hmmsearch

Parallel_run(){
    threads=$1
    query=$2
    outdir=$3
    shift 3

    parallel --link --jobs 2 Run_hmmsearch ::: 1 ::: $threads ::: $@ ::: $query ::: $outdir
}

par=0
outdir="./"

while getopts "ht:q:o:p" option
do
    case $option in
        h)
            Help
            exit 0
            ;;
        t)
            threads=$OPTARG
            ;;
        q)
            query=$OPTARG
            ;;
        o)
            outdir=$OPTARG
            ;;
        p)
            par=1
            ;;
        \?)
            Help
            exit 1
            ;;
        :)
            Help
            exit 1
            ;;
        *)
            Help
            exit 1
            ;;
    esac
done

if [ $outdir == "./" ] || [ $outdir == "." ]
then
    outdir=$(pwd)
elif [ ! -d  $outdir ]
then
    mkdir $outdir
fi

shift $((OPTIND-1))

if [ $# == 0 ]
then
    echo "Error: requires at least one positional argument"
    echo
    Help
    exit 1
fi

if [ $par == 1 ]
then
    if [ `expr $threads % 2` != 0 ]
    then
        echo "Error: if -p is set, the -t value must be even"
        echo
        Help
        exit 1
    else
        echo $(date)
        SECONDS=0
        Parallel_run $threads $query $outdir $@
    fi
else
    echo $(date)
    SECONDS=0
    for db in $@
    do
        Run_hmmsearch $par $threads $db $query $outdir
    done
fi

ELAPSED=$SECONDS
printf 'Elapsed time: %02dh:%02dm:%02ds\n'  $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
