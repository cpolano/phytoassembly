#!/usr/bin/env perl

# blast-ing non-mapping sequences

# Version 0.4

# input:  Healthy.contigs.fasta (db), Stage3.<cutoffval>.contigs.fasta (query) <threads>
# output: Healthy.contigs.bdb (Blast DataBase), Stage3.<cutoffval>.contigs.csv,
#         Stage3.<cutoffval>.contigs.plant.csv, Stage3.<cutoffval>.contigs.phyto.csv,
#         Stage3.<cutoffval>.contigs.phyto.fasta

use strict;
#use warnings;
use Getopt::Std;

my %opts; getopts('h', \%opts);
my $argc = $#ARGV + 1;
if ( ($argc < 3) || ( $opts{h} ) ) {
    print "\nUSAGE:
    
$0 Healthy.contigs.fasta Stage3.<cutoffval>.contigs.fasta <threads>

Blasting and separating phytoplasma and plant genes.

Output: Stage3.<cutoffval>.plant.csv
        Stage3.<cutoffval>.phyto.csv
        Stage3.<cutoffval>.phyto.fasta

OPTIONS:

   -h    show this message\n\n";
    exit 1
}

my $in_database = "$ARGV[0]";
my $query = "$ARGV[1]";
my $threnum = "$ARGV[2]";

# removing spaces from the query file
# if there are spaces in the the IDs,
# blast only keeps the string before the first space
my $bkfile = $query;
$bkfile =~ s/\.fasta/\.bak\.fasta/;
if (-f $bkfile) {
    print "$query backup is present...\n";
} else {
    system("mv $query $bkfile");
}

open(my $INQUERY, "<", $bkfile) or die "$bkfile reading error: ".$!;
open(my $OUTQUERY, ">", $query) or die "$query writing error: ".$!;
while(1) {
    my $row = readline $INQUERY;
    last unless defined $row;
    $row =~ s/ /_/g;
    print $OUTQUERY $row;
}
close($INQUERY);
close($OUTQUERY);

# initializing database
my $database = $in_database.".bdb";
$database =~ s/\.fasta//;
my $dbcheck = $database.".nsq";
if (-f $dbcheck) {
    print "Database $database is present...\n";
} else {
    system("makeblastdb -in $in_database -dbtype nucl -out $database");
}

# blast-ing nucleotides
my $result = $query.".csv";
$result =~ s/\.fasta//;
print "Starting tblastx analysis... ";
if ($threnum == 0) {
    system("tblastx -db $database -query $query -outfmt \"6 qseqid sseqid pident length qseq sseq\" -max_target_seqs 1 -out $result");
} else {
    system("tblastx -db $database -query $query -num_threads $threnum -outfmt \"6 qseqid sseqid pident length qseq sseq\" -max_target_seqs 1 -out $result");
}

print "Done.\nRemoving duplicates in tblastx output... ";
my $bresult = $result.".blast.csv";
$bresult =~ s/\.csv\.blast/\.blast/;
system("mv $result $bresult");
system("sort -k1,1 -k2,2 -k3,3nr $bresult | sort -u -k1,1 --merge > $result");

my %physeqs; my $contIDs;
open(my $PHYTOSEQ, "<", $query) or die "$query reading error: ".$!;
print "Done.\nBuilding contigs database... ";
while(1) {
    my $row = readline $PHYTOSEQ;
    last unless defined $row;
    $row =~ s/\n//; $row =~ s/\r//;
    if ($row =~ m/^>/) {
        $contIDs = $row;
        $contIDs =~ s/\>//;
    } else {
        $physeqs{$contIDs} = $physeqs{$contIDs}.$row;
    }
}
close($PHYTOSEQ);

# initializing output
my $phyfasnam = $result.".phyto.fasta";
$phyfasnam =~ s/\.csv\.phyto/\.phyto/;
my $plant = $result.".plant.csv";
$plant =~ s/\.csv\.plant/\.plant/;
my $phyto = $result.".phyto.csv";
$phyto =~ s/\.csv\.phyto/\.phyto/;

my $OUTPLANT; my $OUTPHYTO; my $OUTFASTA;
open($OUTPLANT, ">", $plant) or die "$OUTPLANT writing error: ".$!;
open($OUTPHYTO, ">", $phyto) or die "$OUTPHYTO writing error: ".$!;
open($OUTFASTA, ">", $phyfasnam) or die "$OUTFASTA writing error: ".$!;
print $OUTPLANT "query id\tsubject id\t\% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\n";
print $OUTPHYTO "query id\tsubject id\t\% identity\talignment length\tmismatches\tgap opens\tq. start\tq. end\ts. start\ts. end\tevalue\tbit score\n";

# divides blast output in two files, one for the plant and one for the phytoplasma
# and exports the phytoplasma sequences in a fasta file

my $singleID = "non-matching tag";
my $discrim = 95;   # matches with more than this value
#                     % identity are assigned to the plant
print "Done.\nDividing plant (>$discrim% identity) and phytoplasma contigs... ";

open(my $BLASTDATA, "<", $result) or die "$result reading error: ".$!;
while(1) {
    my $row = readline $BLASTDATA;
    last unless defined $row;
    my @tabs = split(/\t/, $row);
    if ($tabs[0] =~ m/^#/) {next;}
    if ($singleID eq $tabs[0]) {next;}
    if ($tabs[2] > $discrim) {
        print $OUTPLANT $row;
    } else {
        print $OUTPHYTO $row;
        print $OUTFASTA ">".$tabs[0]."\n";
        print $OUTFASTA $physeqs{$tabs[0]}."\n";
    }
    $singleID = $tabs[0];
}
print "Done.\n";
close($BLASTDATA);
close($OUTFASTA);
close($OUTPHYTO);
close($OUTPLANT);
