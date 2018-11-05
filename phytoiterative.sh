#!/bin/bash

# Version 0.2

usage()
{
cat << EOF

Usage:

phytoiterative.sh [-skipref] [-ref REFERENCE_GENOME.FASTA | -ref REFERENCE_GENOME.FASTQ ] (-ref2 REFERENCE_GENOME.FASTQ_2 ) [-skipreads] (-readsfa READS_CONTIGS.FASTA ) [-reads READS.FASTQ ] (-reads2 READS.FASTQ_2 ) [-threads THREADS ]

This script runs Phytoassembly with increasing cutoffs until it finds an optimum.

OPTIONS:
   -help        show this message

   -skipref     don't assemble the reference genome fastq (you already have an
                assembled fasta file).

   -ref         input the reference genome file as a fasta file if you already
                have assembled it, or as an intelaced fastq file if you didn't.
                If you have paired-end fastq files, put the forward reads file
                here.

   -ref2        input the reverse reference reads file if you have a paired-end
                fastq reference genome.

   -skipreads   don't assemble the diseased plant reads fastq (you already have
                an assembled fasta file).

   -readsfa     input the assembled diseased plant reads as a fasta file, if
                you already have assembled it with a5_pipeline.pl (you'll also
                need to input the *.ec.fastq file with the -reads flag)

   -reads       input the diseased plant reads intelaced fastq, or the forward
                reads file you have paired-end reads files.

   -reads2      input the reverse diseased plant reads file if you have
                paired-end fastq files.

   -threads     (optional) use N threads for the assembly. Default: 1


EXAMPLES:

phytoiterative.sh -ref REFERENCE_GENOME.FASTQ_1 -ref2 REFERENCE_GENOME.FASTQ_2 -reads READS.FASTQ -threads 4

phytoiterative.sh -skipref -ref REFERENCE_GENOME.FASTA -skipreads -readsfa READS_CONTIGS.FASTA -reads READS.FASTQ_1 -reads2 READS.FASTQ_2



EOF
}

if [[ -z $1 ]] ; then usage ; exit 1 ; fi

for arg in "$@" ; do
    shift
    case "$arg" in
        "-help")      set -- "$@" "-h" ;;
        "-skipref")   set -- "$@" "-a" ;;
        "-ref")       set -- "$@" "-b" ;;
        "-ref2")      set -- "$@" "-c" ;;
        "-skipreads") set -- "$@" "-d" ;;
        "-readsfa")   set -- "$@" "-e" ;;
        "-reads")     set -- "$@" "-f" ;;
        "-reads2")    set -- "$@" "-g" ;;
        "-threads")   set -- "$@" "-t" ;;
        *)            set -- "$@" "$arg"
    esac
done

SKIPREF=0
SKIPREADS=0
READS=""
READS2=""
THREADS=0
CUTESTIM=0
OPTIND=1

while getopts "hab:c:de:f:g:t:j:k:l:" opt ; do
    case "$opt" in
        "h") usage ; exit 1;;
        "a") SKIPREF=1;;
        "b") REF=$OPTARG;;
        "c") REF2=$OPTARG;;
        "d") SKIPREADS=1;;
        "e") READSFA=$OPTARG;;
        "f") READS=$OPTARG;;
        "g") READS2=$OPTARG;;
        "t") THREADS=$OPTARG;;
        "?") usage ; exit;;
    esac
done
shift $(($OPTIND - 1))

if [[ "$SKIPREF" -eq 1 ]] ; then
    if [[ "$SKIPREADS" -eq 1 ]] ; then
        if [[ -z $READS2 ]] ; then
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -skipref -ref $REF -skipreads -readsfa $READSFA -reads $READS -min 0 -max 0
            else
                phytoassembly.sh -skipref -ref $REF -skipreads -readsfa $READSFA -reads $READS -threads $THREADS -min 0 -max 0
            fi
        else
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -skipref -ref $REF -skipreads -readsfa $READSFA -reads $READS -reads2 $READS2 -min 0 -max 0
            else
                phytoassembly.sh -skipref -ref $REF -skipreads -readsfa $READSFA -reads $READS -reads2 $READS2 -threads $THREADS -min 0 -max 0
            fi
        fi
    else
        if [[ -z $READS2 ]] ; then
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -skipref -ref $REF -reads $READS -min 0 -max 0
            else
                phytoassembly.sh -skipref -ref $REF -reads $READS -threads $THREADS -min 0 -max 0
            fi
        else
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -skipref -ref $REF -reads $READS -reads2 $READS2 -min 0 -max 0
            else
                phytoassembly.sh -skipref -ref $REF -reads $READS -reads2 $READS2 -threads $THREADS -min 0 -max 0
            fi
        fi
    fi
else
    if [[ "$SKIPREADS" -eq 1 ]] ; then
        if [[ -z $REF2 ]] ; then
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -ref $REF -skipreads -readsfa $READSFA -reads $READS -reads2 $READS2 -min 0 -max 0
            else
                phytoassembly.sh -ref $REF -skipreads -readsfa $READSFA -reads $READS -reads2 $READS2 -threads $THREADS -min 0 -max 0
            fi
        else
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -ref $REF -ref2 $REF2 -skipreads -readsfa $READSFA -reads $READS -reads2 $READS2 -min 0 -max 0
            else
                phytoassembly.sh -ref $REF -ref2 $REF2 -skipreads -readsfa $READSFA -reads $READS -reads2 $READS2 -threads $THREADS -min 0 -max 0
            fi
        fi
    else
        if [[ -z $REF2 ]] ; then
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -ref $REF -reads $READS -reads2 $READS2 -min 0 -max 0
            else
                phytoassembly.sh -ref $REF -reads $READS -reads2 $READS2 -threads $THREADS -min 0 -max 0
            fi
        else
            if [[ -z $THREADS ]] ; then
                phytoassembly.sh -ref $REF -ref2 $REF2 -reads $READS -reads2 $READS2 -min 0 -max 0
            else
                phytoassembly.sh -ref $REF -ref2 $REF2 -reads $READS -reads2 $READS2 -threads $THREADS -min 0 -max 0
            fi
        fi
    fi
fi

gunzip -k Results_0/Stage3.0.contigs.phyto.fasta.gz
nt_size=$(countSequences.pl -r Results_0/Stage3.0.contigs.phyto.fasta)
let nt_size_pre=$nt_size

let nt_limit=$nt_size_pre/200
cutopt=0
sizediff=0

while [ $sizediff -lt $nt_limit ] ; do
    cp Healthy_bkp/Healthy.contigs.fasta .
    cp Diseased_bkp/Diseased.contigs.fasta .
    cp Diseased_bkp/Diseased.ec.fastq.gz .
    gunzip Diseased.ec.fastq.gz

    let cutopt++
    #let nt_size_pre=$nt_size

    if [[ -z $THREADS ]] ; then
        phytoassembly.sh -skipref -ref Healthy.contigs.fasta -skipreads -readsfa Diseased.contigs.fasta -reads Diseased.ec.fastq -min $cutopt -max $cutopt
    else
        phytoassembly.sh -skipref -ref Healthy.contigs.fasta -skipreads -readsfa Diseased.contigs.fasta -reads Diseased.ec.fastq -threads $THREADS -min $cutopt -max $cutopt
    fi

    gunzip Results_$cutopt/Stage3.$cutopt.contigs.phyto.fasta.gz
    nt_size=$(countSequences.pl -r Results_$cutopt/Stage3.$cutopt.contigs.phyto.fasta)
    echo "$sizediff    $nt_size_pre    $nt_size" >> stats.txt
    let sizediff=$nt_size_pre-$nt_size
    #let sizediff=$nt_size_pre-$nt_size
    echo "There is a difference of $sizediff reads between cutoff $cutopt and $[cutopt-1]"
done

cp Results_$[cutopt-1]/Stage3.$[cutopt-1].contigs.phyto.fasta .
echo "The optimal cutoff determined through iteration is $[cutopt-1]"

