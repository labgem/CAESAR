#!/bin/bash

echo $(date)
SECONDS=0

script_path=$1
shift 1

outdir="."
query=""
config=""
blast_path=""
id=30.0
cov=80.0
min_len=200
max_len=1000
tax="ABE"

Help(){
    echo "Script used only with a slurm sbatch command"
    echo
    echo "usage: filter.sh <path>/filter.py [-h] -q -C -d [-i] [-c] [-l] [-L] [-t]"
    echo
    echo "options:"
    echo "  -h    show this help message"
    echo "  -q    set of reference sequences"
    echo "  -C    the configuration file in yaml format"
    echo "  -d    directory containing a blast output or a blast output"
    echo "  -i    percentage identity threshold [default: 30.0]"
    echo "  -c    percentage coverage threshold [default: 80.0]"
    echo "  -l    minimum candidate sequence length [default:200]"
    echo "  -L    maximum candidate sequence length [default: 1000]"
    echo "  -t    selected superkingdoms to filter candidates sequences "
    echo "        A=Archae, B=Bacteria, E=Eukaryota [default: ABE]"
}

export -f Help

while getopts ":o:q:C:d:i:c:l:L:t:" option
do
    case $option in
        o)
            outdir=$OPTARG
            ;;
        q)
            query=$OPTARG
            ;;
        C)
            config=$OPTARG
            ;;
        d)
            blast_path=$OPTARG
            ;;
        i)
            id=$OPTARG
            ;;
        c)
            cov=$OPTARG
            ;;
        l)
            min_len=$OPTARG
            ;;
        L)
            max_len=$OPTARG
            ;;
        t)
            tax=$OPTARG
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

if [ ! -f $outdir ] && [ $outdir != "." ]
then
    mkdir $outdir
fi

python $script_path -o $outdir -q $query -c $config -d $blast_path --id $id --cov $cov --min $min_len --max $max_len --tax $tax

if [ $? != 0 ];
then
    echo an error has occured
	exit 1
fi
