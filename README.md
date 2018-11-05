# Phytoassembly

Phytoplasmas are important plant pathogens that cannot however be cultured. To analyse their DNA, an alternative approach is to extract the sequences from those of a diseased plant. The phytoassembly procedure analyses a diseased plant sequencing to separate the plant reads, using an assembly from a healthy plant reference, from those of the phytoplasma and potential other pathogens present. Subsequent gene annotation can help identifying the pathogen species.

## Requirements
The procedure runs on 64-bit linux kernel 2.6.15 or later; requirements are working installations of [BioPerl](http://www.bioperl.org/wiki/Main_Page), [bwa](http://bio-bwa.sourceforge.net/bwa.shtml), [samtools](http://www.htslib.org/doc/), [blast](http://www.ncbi.nlm.nih.gov/books/NBK279690/) and [a5_pipeline](http://sourceforge.net/p/ngopt/wiki/A5PipelineREADME/).

## Installation
Once the requirements are met, the procedure can be launched from any location, preferably in a folder. The files to be analysed should be in the same folder as the procedure.

## Usage
```
phytoassembly.sh [-skipref] [-ref REFERENCE_GENOME.FASTA | -ref REFERENCE_GENOME.FASTQ ] (-ref2 REFERENCE_GENOME.FASTQ_2 ) [-skipreads] (-readsfa READS_CONTIGS.FASTA ) [-reads READS.FASTQ ] (-reads2 READS.FASTQ_2 ) [-threads THREADS ] [-min CUTOFFMIN ] [-max CUTOFFMAX ] [-step CUTOFFSTEP ]
```

## Methods
The procedure requires two files as input: a reference genome for the uninfected plant in fasta format, and the sequence reads to be analysed in fastq format. The reads can be either be stored in a single interleaved file or in two paired-end files, and the pipeline can also read files compressed in the gzip format.

0. The procedure calls the A5 pipeline to assembly the reference genome reads. The contigs list (`Healthy.contigs.fasta`) is kept for later, while the remaining files are archived in the bzip2 format.

1. The procedure calls the A5 pipeline to assembly the mixed genome reads. The contigs list (`Diseased.contigs.fasta`) is then indexed and aligned with the error corrected reads (`Diseased.ec.fastq`) using the bwa index and mem commands.
The resulting sam file (`Diseased.mapped.sam`) is converted to the bam format, and statistical data is produced using the samtools commands view, sort, index and idxstats, in order.

2. The resulting database (`Diseased.sorted.csv`) is passed to a perl procedure (`phytocount.pl`) that calculates the coverages (as a ratio between the number of mapped reads and the length, in %), and saves them in a csv file (`Diseased.sorted.cov.csv`). The contigs are then grouped according to coverage ratios (1%, 2%, 3%, etc.) and saved in a csv file (`Diseased.sorted.histo.csv`).

The following steps are run once with no cutoff, then using an optimized cutoff (`$cutoffval`) determined using the ratio between the sum of the lengths of the non-mapping reads with no cutoff and the sum of the lengths of the error corrected reads of the diseased plant.

3. Contigs with coverages higher than the cutoff value are exported to a fasta file (`Diseased.cutoff.$cutoffval.fasta`). This file is again indexed and aligned with the high quality reads (`Diseased.ec.fastq`) with bwa, again using the mem command, and saved in the sam format (as `Stage1.$cutoffval.match.sam`).

4. The sam file is passed to a second perl procedure (`phytofilter.pl`) that extracts and exports both mapping and non-mapping sequences to fasta and fastq files (`Stage1.$cutoffval.nonmatch.fasta/q`, `Stage1.$cutoffval.match.fasta/q`), according to the sam flag #4, (“the query sequence itself is unmapped”).

5. The mapping reads (`Stage1.$cutoffval.match.fastq`) are then aligned with bwa mem against the healthy assembly (`Healthy.contigs.fasta`) and saved in the sam format (`Stage2.$cutoffval.match.sam`).  
The sam file is passed to `phytofilter.pl` again, producing mapping and non-mapping fasta and fastq files (`Stage2.$cutoffval.match.fasta/q` and `Stage2.$cutoffval.nonmatch.fasta/q`).

6. The non-mapping sequences are assembled with the A5 pipeline, to obtain the non-mapping contigs list (`Stage3.$cutoffval.contigs.fasta`).

7. A third perl procedure (`phytoblast.pl`) creates a blast nucleotide database from the reference genome contigs file (`Healthy.contigs.fasta`) and queries it with tblastx (translated nucleotide query vs. translated nucleotide database blast) using the non-mapping sequences (`Stage3.$cutoffval.contigs.fasta`).  
The results are saved in a csv file (`Stage3.$cutoffval.contigs.csv`), which is then filtered according to the identity percentage (I.P.): entries with an I.P. greater than 95% are attributed to the plant (`Stage3.$cutoffval.contigs.plant.csv`), while those with a lower I.P. are attributed to the phytoplasma (`Stage3.$cutoffval.contigs.phyto.csv`).  
The IDs saved in the last csv file are used to extract the phytoplasma sequences from the query and saved in a fasta file (`Stage3.$cutoffval.contigs.phyto.fasta`).

8. Lastly, all the outputs (Diseased.*, Healthy.*, Stage3.*) are moved to folders with a timestamp (e.g. `PhytoFiles_<cutoff>_<date>`) and compressed in the tbz2 format.
