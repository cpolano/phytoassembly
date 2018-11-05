#!/usr/bin/env perl

use strict;
use warnings;
use Getopt::Std;
use Bio::Perl;
use Bio::SeqIO;
use List::Util qw[min max];

my %opts; getopts('vrh', \%opts);
my $argc = $#ARGV + 1;
if ( ($argc < 1) || ( $opts{h} ) ) {
    print "\nUSAGE:

$0 [h|v|r] infile.fasta|infile.fastq

This program counts sequences and nucleotides
in a fastA or fastQ file.

OPTIONS:

   -v    more detailed statistics
   -r    write the reads total count in a temporary file
   -h    show this message\n\n";
    exit 1
}

my $inFile = $ARGV[0];
my $seqIN;


if (($inFile =~ /\.fa$/) ||
    ($inFile =~ /\.faa$/) ||
    ($inFile =~ /\.fna$/) ||
    ($inFile =~ /\.fas$/) ||
    ($inFile =~ /\.fasta$/)) {
    $seqIN = Bio::SeqIO->new('-format' => 'fasta', '-file' =>  "$inFile");
} elsif ($inFile =~ /\.fastq$/) {
    $seqIN = Bio::SeqIO->new('-format' => 'fastq', '-file' =>  "$inFile");
} else {
    print "Wrong format! Should be .fasta or .fastq\n\n"; exit 1;
}

my $count=0;
my $tot=0;
my $x;
my $g = my $c = 0;
my @contigL;
my @contigNames;
my $prog = 0;


print "File: $inFile" if ( !$opts{r} );
while(my $seq = $seqIN->next_seq()) {

    if    ($prog < 1000)     { $prog = int($count/100)*100 }
    elsif ($prog < 10000)    { $prog = int($count/1000)*1000 }
    elsif ($prog < 100000)   { $prog = int($count/10000)*10000 }
    elsif ($prog < 1000000)  { $prog = int($count/100000)*100000 }
    elsif ($prog < 10000000) { $prog = int($count/1000000)*1000000 }

    $count ++;
    $x = $seq->length();
    push @contigNames, $seq->display_id;
    $tot=$tot+$x;
    push(@contigL, $x);
    $c += $seq->seq =~ tr/C//;
    $c += $seq->seq =~ tr/c//;
    $g += $seq->seq =~ tr/G//;
    $g += $seq->seq =~ tr/g//;
}

my $gc = $c + $g;
$gc = $gc / $tot * 100;


if ( $opts{r} ) {
    print $tot; exit 1;
} elsif ( !$opts{v} ) {
    printf "\tsequences: %d\tnt|aa: %d\tmean: %.1f\tG+C: %.1f%%\n\n", $count, $tot, $tot/$count, $gc;
    exit 1;
} else {
    printf "\nThere are %d sequences for a total of %d nucleotides/aminoacids (mean = %.1f)\n", $count, $tot, $tot/$count;
    my @his = (0,0,0,0,0,0,0,0,0,0,0,0,0,0);
    my @contigA = sort { $b <=> $a } @contigL;
    my $a=0;
    my $j;
    my $N50 =0;
    my $N50c =1;
    my $N95 =0;
    my $newval =0;

    # log classes: 10-100 ; 100-1000; 1000-10000; 10000-100000;
    for ($j = 0; $j <= scalar(@contigA)-1; $j++) {

        if    ($prog < 1000)     { $prog = int($j/100)*100 }
        elsif ($prog < 10000)    { $prog = int($j/1000)*1000 }
        elsif ($prog < 100000)   { $prog = int($j/10000)*10000 }
        elsif ($prog < 1000000)  { $prog = int($j/100000)*100000 }
        elsif ($prog < 10000000) { $prog = int($j/1000000)*1000000 }

        if($contigA[$j] != 0) {
            $a = int(log($contigA[$j])/log(10));
            $his[$a]++; 
            $newval=$contigA[$j]+$newval;
            if($N50 == 0) {
                $N50c++;
                if($newval >= ($tot/2)) {
                    $N50 = $contigA[$j];
                }
            }
            if($N95 == 0) {
                if($newval >= ($tot * 0.95)) {
                    $N95 = $contigA[$j];
                }
            }
        } else {
            print "at position $j found 1 seq with 0 len : name $contigNames[$j] , prev name $contigNames[$j-1]\n";
        }
    }
    print "     10 -     100 -> $his[1] \n";
    print "    100 -   1,000 -> $his[2] \n";
    print "  1,000 -  10,000 -> $his[3] \n";
    print " 10,000 - 100,000 -> $his[4] \n";
    print "        > 100,000 -> $his[5] \n";

    printf "min: %d  max: %d  N50: %d  N50 contigs: %d  G+C: %.1f%%\n\n", min(@contigL), max(@contigL), $N50, $N50c, $gc;
}
