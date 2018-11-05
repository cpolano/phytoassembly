#!/bin/bash

# Version 0.9.2

# This pipeline analyzes a mixed plant/phytoplasma genome and separates
# the phytoplasma genes from those of the plant.

# Runs on 64-bit linux kernel 2.6.15 or later; requirements are working
# installations of:
#     BioPerl
#     BWA
#     Samtools
#     NCBI Blast+
#     A5 pipeline (https://sourceforge.net/projects/ngopt/)

# Subprograms are:
#     phytocount.pl
#     phytofilter.pl
#     phytoblast.pl

usage()
{
cat << EOF

Usage:

phytoassembly.sh [-skipref] [-ref REFERENCE_GENOME.FASTA | -ref REFERENCE_GENOME.FASTQ ] (-ref2 REFERENCE_GENOME.FASTQ_2 ) [-skipreads] (-readsfa READS_CONTIGS.FASTA ) [-reads READS.FASTQ ] (-reads2 READS.FASTQ_2 ) [-threads THREADS ] [-min CUTOFFMIN ] [-max CUTOFFMAX ] [-step CUTOFFSTEP ]

This procedure analyzes a mixed plant/phytoplasma genome and separates the phytoplasma genes from those of the plant.

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

   -threads     (optional) use N threads for the assembly.
  
   -min         (optional) cutoff min value.
  
   -max         (optional) cutoff max value.
  
   -step        (optional) cutoff step.


EXAMPLES:

phytoassembly.sh -ref REFERENCE_GENOME.FASTQ_1 -ref2 REFERENCE_GENOME.FASTQ_2 -reads READS.FASTQ -threads 4 -min 0 -max 20 -step 5

phytoassembly.sh -skipref -ref REFERENCE_GENOME.FASTA -skipreads -readsfa READS_CONTIGS.FASTA -reads READS.FASTQ_1 -reads2 READS.FASTQ_2 -min 5 -max 25 -step 10



EOF
}

#--------------------------------------------------------------------------------
extractPhyto () {
    let folderids[$count]=$1
    echo "3. Calculating the coverage with cutoff at $1..."
    echo "
3. Calculating the coverage with cutoff at $1" >> $infotxt
    datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt
    phytocount.pl Diseased.sorted.csv Diseased.contigs.fasta $1 2>> $infotxt
    cutoffasta=Diseased.cutoff.$1.fasta
    err=$? ; datenow=$(date +"%H:%M:%S") ;
    echo "END: $datenow with error status $err" >> $infotxt
    echo "4. Extracting mapping reads..."
    echo "
4. Extracting mapping reads" >> $infotxt
    datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt
    bwa index $cutoffasta 2>> $infotxt
    
    if [[ $THREADS -eq 0 ]] ; then
        bwa mem $cutoffasta Diseased.ec.fastq > Stage1.$1.match.sam 2>> $infotxt
    else
        bwa mem -t $THREADS $cutoffasta Diseased.ec.fastq > Stage1.$1.match.sam 2>> $infotxt
    fi
    samtools view -bS Stage1.$1.match.sam > Stage1.$1.match.bam 2>> $infotxt

    # The following step is equivalent to "samtools fastq", but A5_pipeline.pl is bundled
    # with an older version of samtools, so just to be safe this ad-hoc method is used
    phytofilter.pl -aq Stage1.$1.match.sam 2>> $infotxt
    err=$? ; datenow=$(date +"%H:%M:%S");
    echo "END: $datenow with error status $err" >> $infotxt

    echo "5. Re-aligning mapping reads..."
    echo "
5. Re-aligning mapping reads" >> $infotxt
    datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt
    bwa index Healthy.contigs.fasta 2>> $infotxt
    if [[ $THREADS -eq 0 ]] ; then
        bwa mem Healthy.contigs.fasta Stage1.$1.match.fastq > Stage2.$1.match.sam 2>> $infotxt
    else
        bwa mem -t $THREADS Healthy.contigs.fasta Stage1.$1.match.fastq > Stage2.$1.match.sam 2>> $infotxt
    fi
    samtools view -bS Stage2.$1.match.sam > Stage2.$1.match.bam 2>> $infotxt
    phytofilter.pl -aq Stage2.$1.match.sam 2>> $infotxt
    
    err=$? ; datenow=$(date +"%H:%M:%S");
    echo "END: $datenow with error status $err" >> $infotxt
    
    if [[ $CUTESTIM -eq 1 ]] ; then return 0 ; fi

    echo "6. Assembling the selected reads..."
    echo "
6. Assembling the selected reads" >> $infotxt
    datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt
    if [[ $THREADS -eq 0 ]]; then
        a5_pipeline.pl --end=2 Stage2.$1.nonmatch.fastq Stage3.$1 2>> $infotxt
    else
        a5_pipeline.pl --threads=$THREADS --end=2 Stage2.$1.nonmatch.fastq Stage3.$1 2>> $infotxt
    fi
    err=$? ; datenow=$(date +"%H:%M:%S");
    echo "END: $datenow with error status $err" >> $infotxt

    echo "7. Blast-ing non-mapping reads..."
    echo "
7. Blast-ing non-mapping reads" >> $infotxt
    datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt
    phytoblast.pl Healthy.contigs.fasta Stage3.$1.contigs.fasta $THREADS 2>> $infotxt
    err=$? ; datenow=$(date +"%H:%M:%S");
    echo "END: $datenow with error status $err" >> $infotxt

    let count++
    return 0
}
#--------------------------------------------------------------------------------

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
        "-min")       set -- "$@" "-j" ;;
        "-max")       set -- "$@" "-k" ;;
        "-step")      set -- "$@" "-l" ;;
        *)            set -- "$@" "$arg"
    esac
done

SKIPREF=0
SKIPREADS=0
READS=""
READS2=""
THREADS=0
CUTMIN="estimated"
CUTMAX="estimated"
CUTSTEP="estimated"
CUTESTIM=0
OPTIND=1

while getopts "hab:c:de:f:g:t:j:k:l:" opt ; do
    case "$opt" in
        "h") usage; exit 1;;
        "a") SKIPREF=1;;
        "b") REF=$OPTARG;;
        "c") REF2=$OPTARG;;
        "d") SKIPREADS=1;;
        "e") READSFA=$OPTARG;;
        "f") READS=$OPTARG;;
        "g") READS2=$OPTARG;;
        "t") THREADS=$OPTARG;;
        "j") CUTMIN=$OPTARG;;
        "k") CUTMAX=$OPTARG;;
        "l") CUTSTEP=$OPTARG;;
        "?") usage; exit;;
    esac
done
shift $(($OPTIND - 1))


refext=`echo "${REF##*.}"`
ref2ext=`echo "${REF2##*.}"`
readext=`echo "${READSFA##*.}"`
read1ext=`echo "${READS##*.}"`
read2ext=`echo "${READS2##*.}"`

if [[ "$SKIPREF" -eq 1 ]] && [[ "$refext" != "fasta" ]] && [[ "$refext" != "fas" ]] && [[ "$refext" != "fna" ]] && [[ "$refext" != "fa" ]]; then
    echo "The reference genome file must be input as a fasta file!"; exit 1
elif [[ "$SKIPREF" -eq 0 ]] && [[ "$refext" != "fastq" ]] && [[ "$refext" != "fq" ]]; then
    echo "The reference genome file must be input as a fastq file!"; exit 1
elif [[ -n "$ref2ext" ]] && [[ "$ref2ext" != "fastq" ]] && [[ "$ref2ext" != "fq" ]]; then
    echo "The reference genome right-hand file must be input as a fastq file!"; exit 1
fi

if [[ "$SKIPREADS" -eq 1 ]] && [[ "$readext" != "fasta" ]] && [[ "$readext" != "fas" ]] && [[ "$readext" != "fna" ]] && [[ "$readext" != "fa" ]]; then
    echo "Contigs file must be input as a .fasta file!"; exit 1
fi

if [[ -n "$read1ext" ]] && [[ "$read1ext" != "fastq" ]] && [[ "$read1ext" != "fq" ]]; then
    echo "Contigs file must be input as a .fastq file!"; exit 1
fi

if [[ "$CUTMIN" != "estimated" ]] ; then
    if [[ "$CUTMAX" != "estimated" ]] ; then
        if [[ "$CUTMAX" -lt "$CUTMIN" ]]; then
            let CUTSWAP=$CUTMAX; let CUTMAX=$CUTMIN; let CUTMIN=$CUTSWAP
        fi
        let CUTDIFF=$CUTMAX-$CUTMIN
        if [[ "$CUTSTEP" == "estimated" ]]; then
            if [[ "$CUTDIFF" -lt 1 ]]; then
                let CUTSTEP=1
            else
                let CUTSTEP=$CUTDIFF
            fi
        elif [[ "$CUTSTEP" -lt 1 ]]; then
            echo "Cutoff steps must be at least 1!"
            let CUTSTEP=1
        elif [[ "$CUTMIN" -lt 0 ]] || [[ "$CUTMAX" -lt 0 ]]; then
            echo "Cutoff must be greater than zero!"; exit 1
        fi
    else
        echo "Maximum cutoff is not defined!"; exit 1
    fi
else
    if [[ "$CUTMAX" != "estimated" ]] ; then
        echo "Minimum cutoff is not defined!"; exit 1
    else
        let CUTESTIM=1
    fi
fi

infotxt=info_`date +"%y%m%d"`.txt
echo "
*** PhytoAssembly: start ***
"

echo "*****" >> $infotxt
echo "Pipeline: MiSeq A5 pipeline" >> $infotxt
echo "Reference files: $REF $REF2" >> $infotxt
echo "Reads files: $READSFA $READS $READS2" >> $infotxt
echo "*****" >> $infotxt
if [[ $SKIPREF -eq 1 ]]; then
    echo "Assemble reference genome: no" >> $infotxt
else
    echo "Assemble reference genome: yes" >> $infotxt
fi
if [[ $SKIPREADS -eq 1 ]]; then
    echo "Assemble reads: no" >> $infotxt
else
    echo "Assemble reads yes" >> $infotxt
fi
echo -n "Cutoffs:" >> $infotxt
if [[ $CUTESTIM -eq 1 ]]; then
    echo -n " estimated" >> $infotxt
else
    for (( cutoffval=$CUTMIN; cutoffval<=$CUTMAX; cutoffval+=$CUTSTEP )); do
        echo -n " $cutoffval" >> $infotxt
    done
fi
echo >> $infotxt
echo "*****" >> $infotxt

# 0. Assemble the healthy plant genome
if [[ $SKIPREF -eq 0 ]]; then
    echo "0. Assembling the healthy plant genome..."
    datenow=$(date +"%H:%M:%S"); echo "
Assembling the healthy plant genome: $datenow" >> $infotxt
    if [[ $THREADS -eq 0 ]]; then
        a5_pipeline.pl --end=2 $REF $REF2 Healthy 2>> $infotxt
    else
        a5_pipeline.pl --threads=$THREADS --end=2 $REF $REF2 Healthy 2>> $infotxt
    fi
    err=$? ; datenow=$(date +"%H:%M:%S") ; echo "END: $datenow with error status $err" >> $infotxt
    countSequences.pl Healthy.contigs.fasta >> $infotxt
    mkdir Healthy_bkp
    mv Healthy.* Healthy_bkp
    cp Healthy_bkp/Healthy.contigs.fasta .
else
    if [[ $REF != "Healthy.contigs.fasta" ]]; then
        cp $REF Healthy.contigs.fasta
    fi
fi

if [[ $SKIPREADS -eq 0 ]]; then
    echo "0. Assembling the diseased plant genome..."
    echo "
0. Assembling the diseased plant genome" >> $infotxt
    datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt
    if [[ $THREADS -eq 0 ]]; then
        a5_pipeline.pl --end=2 $READS $READS2 Diseased 2>> $infotxt
    else
        a5_pipeline.pl --threads=$THREADS --end=2 $READS $READS2 Diseased 2>> $infotxt
    fi
    err=$? ; datenow=$(date +"%H:%M:%S") ; echo "END: $datenow with error status $err" >> $infotxt
    countSequences.pl Diseased.contigs.fasta >> $infotxt
    mkdir Diseased_bkp
    mv Diseased.* Diseased_bkp
    cp Diseased_bkp/Diseased.contigs.fasta .
    cp Diseased_bkp/Diseased.ec.fastq.gz .
    gunzip Diseased.ec.fastq.gz
else
    if [[ $READS != "Diseased.ec.fastq" ]]; then
        cp $READS Diseased.ec.fastq
    fi
    if [[ $READSFA != "Diseased.contigs.fasta" ]]; then
        cp $READSFA Diseased.contigs.fasta
    fi
fi

echo "1. BWA indexing and alignment..."
echo "
1. BWA indexing and alignment" >> $infotxt
datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt

# Index database sequences in the FASTA format.
bwa index Diseased.contigs.fasta 2>> $infotxt

# "align 70bp-1Mbp query sequences with the BWA-MEM algorithm"
if [[ $THREADS -eq 0 ]]; then
    bwa mem Diseased.contigs.fasta Diseased.ec.fastq > Diseased.match.sam 2>> $infotxt
else
    bwa mem -t $THREADS Diseased.contigs.fasta Diseased.ec.fastq > Diseased.match.sam 2>> $infotxt
fi
err=$? ; datenow=$(date +"%H:%M:%S") ; echo "END: $datenow with error status $err" >> $infotxt

echo "2. samtools sorting and indexing..."
echo "
2. samtools sorting and indexing" >> $infotxt
datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt

# Convert SAM into BAM
samtools view -bS Diseased.match.sam > Diseased.match.bam 2>> $infotxt

# "sort a BAM file by leftmost chromosomal coordinates"
# A5_pipeline.pl is bundled with an older version of samtools,
# not compatible with the 1.3.1+ syntax
samtools 2> stvers.txt
SAMVERS=$(grep Version stvers.txt)
SAMVERS=${SAMVERS:9:3}
let SAMVERS=${SAMVERS//.}
if [[ "$SAMVERS" < 13 ]]; then
    samtools sort Diseased.match.bam Diseased.sorted 2>> $infotxt
else
    samtools sort Diseased.match.bam -o Diseased.sorted.bam 2>> $infotxt
fi
rm stvers.txt

# "index a coordinate-sorted BAM file for fast random access"
samtools index Diseased.sorted.bam 2>> $infotxt

# Retrieve and print stats in the index file
# (ref seq name, seq length, # match reads, # unmatch reads)
samtools idxstats Diseased.sorted.bam > Diseased.sorted.csv 2>> $infotxt
err=$? ; datenow=$(date +"%H:%M:%S") ; echo "END: $datenow with error status $err" >> $infotxt

count=0
folderids[$count]=0

if [[ $CUTESTIM -eq 1 ]]; then
    echo ">  Starting cutoff value estimation..."
    extractPhyto 0

    let CUTMAXEST=$(countSequences.pl -r Diseased.ec.fastq)
    let NONMATCHDIM=$(countSequences.pl -r Stage2.0.nonmatch.fastq)
    if [[ $NONMATCHDIM -eq 0 ]]; then
        echo "There was a problem with the non-matching read size, exiting..."
        exit 1
    fi

    let CUTMAXEST=$NONMATCHDIM*10000/$CUTMAXEST
    let CUTMIN=0
    let CUTMAX=$CUTMAXEST
    let CUTSTEP=$(perl -w -e "use POSIX; print int($CUTMAX*30/10000+0.5), qq{\n}")

    if [[ $CUTSTEP -lt 1 ]]; then let CUTSTEP=1; fi
    let CUTMIN=$CUTSTEP
    let CUTMAX=$CUTSTEP
    echo "Estimated cutoff value: $CUTSTEP"
    let CUTESTIM=2
fi


#----------------MAIN LOOP-------------------------

for (( cutoffval=$CUTMIN; cutoffval<=$CUTMAX; cutoffval+=$CUTSTEP ))
do extractPhyto $cutoffval
done

#--------------------------------------------------


echo "8a. Calculating fasta/q sizes statistics..."
echo "
8a. Calculating fasta/q sizes statistics" >> $infotxt
datenow=$(date +"%H:%M:%S") ; echo "START: $datenow" >> $infotxt
#fd_results=Results_`date +"%y%m%d-%H%M"`
let cutoffval--
echo cutoff $cutoffval
fd_results=Results_$cutoffval

rm -r *.sam *.amb *.ann *.bwt *.pac *.sa *.bai *.s1 *.s2

mkdir $fd_results
cp Diseased.ec.fastq $fd_results/
cp Healthy.contigs.fasta $fd_results/

for (( cutoffval=$CUTMIN; cutoffval<=$CUTMAX; cutoffval+=$CUTSTEP )); do
    cp Stage1.$cutoffval.match.fastq $fd_results/
    cp Stage2.$cutoffval.nonmatch.fastq $fd_results/
    cp Stage3.$cutoffval.contigs.fasta $fd_results/
    cp Stage3.$cutoffval.contigs.phyto.fasta $fd_results/
done
mv $infotxt $fd_results/
cd $fd_results/

echo "SIZE    FILENAME" >> stats.txt
du -sh Diseased.ec.fastq >> stats.txt
du -sh Healthy.contigs.fasta >> stats.txt
for (( cutoffval=$CUTMIN; cutoffval<=$CUTMAX; cutoffval+=$CUTSTEP )); do
    du -sh Stage1.$cutoffval.match.fastq >> stats.txt
    du -sh Stage2.$cutoffval.nonmatch.fastq >> stats.txt
    du -sh Stage3.$cutoffval.contigs.fasta >> stats.txt
    du -sh Stage3.$cutoffval.contigs.phyto.fasta >> stats.txt
    echo >> stats.txt
done
echo "SEQUENCE COUNT" >> stats.txt
countSequences.pl Diseased.ec.fastq >> stats.txt
countSequences.pl Healthy.contigs.fasta >> stats.txt
for (( cutoffval=$CUTMIN; cutoffval<=$CUTMAX; cutoffval+=$CUTSTEP )); do
    countSequences.pl Stage1.$cutoffval.match.fastq >> stats.txt
    countSequences.pl Stage2.$cutoffval.nonmatch.fastq >> stats.txt
    countSequences.pl Stage3.$cutoffval.contigs.fasta >> stats.txt
    countSequences.pl Stage3.$cutoffval.contigs.phyto.fasta >> stats.txt
done

echo "8b. Compressing fasta/q files..."
echo "
8b. Compressing fasta/q files" >> $infotxt

for file in *.fasta; do gzip $file; done
for file in *.fastq; do gzip $file; done

mkdir fastq_match; mv *.match.* fastq_match
mkdir fastq_nonmatch; mv *.nonmatch.* fastq_nonmatch
mkdir fasta_contigs; mv Stage3.*.contigs.fasta.gz fasta_contigs

cd ..
mkdir Other_files
mv Stage* Diseased.* Healthy.* Other_files/

cd Other_files/
for file in *.fasta; do gzip $file; done
for file in *.fastq; do gzip $file; done
cd ..

err=$? ; datenow=$(date +"%H:%M:%S") ;
echo "END: $datenow with error status $err"

mv Other_files $fd_results/

echo "
*** Phytoassembly: end ***
"


