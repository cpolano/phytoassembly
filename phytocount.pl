#!/usr/bin/env perl

use strict;
use Getopt::Std;

my %opts; getopts('h', \%opts);
my $argc = $#ARGV + 1;
if ( ($argc < 3) || ( $opts{h} ) ) {
    print "\nUSAGE:
    
$0 Diseased.sorted.csv Diseased.contigs.fasta <cutoffval>

Coverage calculation and histogram output.

Output: Diseased.sorted.cov.csv
        Diseased.sorted.histo.csv
        Diseased.cutoff.<cutoffval>.fasta

OPTIONS:

   -h    show this message\n\n";
    exit 1
}

# initializing variables and output
my $checkfile = $ARGV[0];
if ($checkfile !~ /\.csv$/) { print "\nThe first file should end with .csv\n\n"; exit 1; }

$checkfile = $ARGV[1];
if ($checkfile !~ /\.fasta$/) { print "\nThe second file should end with .fasta\n\n"; exit 1; }

my $cutoff = $ARGV[2];

# initializing file handles
my $contigsIN; 
my $INFILE; my $OUTFILE1; my $OUTFILE2; my $OUTFASTA;
my %contigs; my $contID;

# reading input
open($contigsIN, "<", $ARGV[1]) or die "$ARGV[1] reading error: ".$!;
while(1) {
    my $seq = readline $contigsIN;
    last unless defined $seq;
    $seq =~ s/\n//; $seq =~ s/\r//;
    if ($seq =~ m/^>/) {
        $contID = $seq;
        $contID =~ s/\>//;
    } else {
        $contigs{$contID} = $seq;
    }
}
close($contigsIN);

# initializing csv output
my $filename1 = $ARGV[0].".cov.csv";
$filename1 =~ s/\.csv\.cov/\.cov/;
my $filename2 = $ARGV[0].".histo.csv";
$filename2 =~ s/\.csv\.histo/\.histo/;

# from the samtools idxstats output calculate coverage %
my $i = 0; my $row; my $maxclass = 0;
my @cl_count; my @cl_names; my @lines; my @rd_count;

open($INFILE, "<", $ARGV[0]) or die "$ARGV[0] reading error: ".$!;
open($OUTFILE1, ">", $filename1) or die "$OUTFILE1 writing error: ".$!;
print $OUTFILE1 "NAME\tLENGTH\tALIGNED_READS\tCOVERAGE\n";
while(1) {
    $row = readline $INFILE;
    last unless defined $row;
    my @tabs = split(/\t/, $row);
    next if ($tabs[0] eq "*");
    my $coverage = $tabs[2] / $tabs[1] * 100;
    my $class = int($coverage + 0.5);
    print $OUTFILE1 $tabs[0]."\t".$tabs[1]."\t".$tabs[2]."\t".$class."\n";

    if ($cl_count[$class]) {
        $cl_count[$class]++;
        $cl_names[$class] = $cl_names[$class].",".$tabs[0];
        $rd_count[$class] = $rd_count[$class] + $tabs[2];
    }
    else {
        $cl_count[$class] = 1;
        $cl_names[$class] = $tabs[0];
        $rd_count[$class] = $tabs[2];
    }
    $maxclass = $class if ($maxclass < $class);
    $lines[$i] = $tabs[0]."\t".$tabs[1]."\t".$tabs[2]."\t".$class;
    $i++;
}
close($OUTFILE1);
close($INFILE);
print "Calculated coverages for $i contigs.\n";
print "Max coverage value: $maxclass\n";
print $filename1." generated.\n";


# saving in .csv
open($OUTFILE2, ">", $filename2) or die "$OUTFILE2 writing error: ".$!;
print $OUTFILE2 "CLASS\tCOUNT\tREADS\tCONTIGS\n";
my $j = $i = 0;
foreach my $row ( @cl_count ) {
    if ($row) {
        print $OUTFILE2 $i."\t".$row."\t".$rd_count[$i]."\t".$cl_names[$i]."\n";
        $j++;
    }
    $i++;
}
close($OUTFILE2);
print "$j classes determined.\n";
print $filename2." generated.\n";


# saving in .fasta
my $OUTFASTA;
my $outfastaname = $ARGV[1];
$outfastaname =~ s/contigs.fasta//;
$outfastaname=$outfastaname."cutoff.".$cutoff.".fasta";

open($OUTFASTA, ">", $outfastaname) or die "$outfastaname writing error: ".$!;

for ($i=$cutoff; $i <= $maxclass; $i++) {
    # generating the fasta file, keeping only the contigs above $cutoff % coverage
    my @IDs = split(/,/, $cl_names[$i]);
    foreach my $id ( @IDs ) {
        print $OUTFASTA ">".$id."\n".$contigs{$id}."\n";
    }
}
close($OUTFASTA);
print $outfastaname." generated.\n";


