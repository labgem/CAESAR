#!/bin/bash

Help(){
    echo "Runs a blastp with diamond software on one or several diamond database"
    echo
    echo "usage: blastp.sh [-h] -q [-o] [-t] [-i] [-c] [-p] <db.dmnd>"
    echo
    echo "options:"
    echo "  -h    show this help message"
    echo "  -q    set of reference sequences"
    echo "  -o    output directory [default: './']"
    echo "  -t    number of threads [default: 1]"
    echo "  -i    minimum %identity [default: 30.0]"
    echo "  -c    minimum %coverage [default: 80.0]"
    echo "  -p    uses gnu parallel to run 2 jobs, each job uses n threads"
    echo "        corresponding to the value of the -t option divide by 2"
}

export -f Help

Diamond_blastp(){

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
    id=$6
    cov=$7

    echo $n

    source=$(echo $db | grep -oE '([A-Za-z0-9_]+)\.dmnd' | cut -d . -f 1)
	source=$(echo "${source,,}")

    diamond blastp --threads $n --very-sensitive -d $db -q $query -o $outdir/matches_${source}.tsv --id $id --query-cover $cov -k0 --outfmt 6 qseqid qlen sseqid slen length pident qcovhsp positive mismatch gaps evalue

}

export -f Diamond_blastp

Parallel_run(){
    threads=$2
    query=$3
    outdir=$4
    id=$5
    cov=$6
    shift 6

    parallel --link --jobs 2 Diamond_blastp ::: 1 ::: $threads ::: $@ ::: $query ::: $outdir ::: $id ::: $cov
}

export -f Parallel_run

threads=1
query=""
outdir="."
id=30.0
cov=80.0
par=0

while getopts "ht:q:o:i:c:p" option
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
        i)
            id=$OPTARG
            ;;
        c)
            cov=$OPTARG
            ;;
        p)
            par=1
            ;;
        \?) # Invalid option
            Help
            exit 1
            ;;
        :) # value required
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

if [ ! -f $outdir ] && [ $outdir != "." ]
then
    mkdir $outdir
fi

if [ ! -f $query ] || [ -z $query ]
then
    echo "Error: argument -q is required and muste be an existing file"
    echo
    Help
    exit 1
fi

if [ $# == 0 ]
then
    echo "requires at least one positional argument"
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
        Parallel_run 1 $threads $query $outdir $id $cov $@
    fi
else
    echo $(date)
    SECONDS=0
    for db in $@
    do
        Diamond_blastp 0 $threads $db $query $outdir $id $cov
    done
fi

ELAPSED=$SECONDS
printf 'Elapsed time: %02dh:%02dm:%02ds\n'  $((ELAPSED/3600)) $((ELAPSED%3600/60)) $((ELAPSED%60))
