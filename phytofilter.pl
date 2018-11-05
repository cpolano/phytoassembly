#!/usr/bin/env perl

use strict;
use Getopt::Std;

my %opts; getopts('aqh', \%opts);
my $argc = $#ARGV + 1;
if ( ($argc < 1) || ( $opts{h} ) ) {
    print "\nUSAGE:
    
$0 StageX.<cutoffval>.match.sam

Separates non-mapping reads.

Output: StageX.<cutoffval>.nonmatch.fasta
        StageX.<cutoffval>.match.fasta
        StageX.<cutoffval>.nonmatch.fastq
        StageX.<cutoffval>.match.fastq

OPTIONS:

   -a    output fastA files
   -q    output fastQ files
   -h    show this message\n\n";
    exit 1
}

my $checkfile = $ARGV[0];
if ($checkfile !~ /\.sam$/) {
    print "\nInput file should end with .sam\n\n"; exit 1;
}

if ( !$opts{a} && !$opts{q} ) {
    print "\nYou need to choose at least one output type:\n -a (fasta) and/or -q (fastq)\n\n"; exit 1;
}

# initializing file handles
my $INFILE;
my $OUT_UNM_Q; my $OUT_MAP_Q;
my $OUT_UNM_A; my $OUT_MAP_A;

# initializing variables and output
my $unmap_q = $ARGV[0].".fastq";
$unmap_q =~ s/\.match\.sam/\.nonmatch/;
my $unmap_a = $ARGV[0].".fasta";
$unmap_a =~ s/\.match\.sam/\.nonmatch/;
my $map_q = $ARGV[0].".fastq";
$map_q =~ s/\.match\.sam/\.match/;
my $map_a = $ARGV[0].".fasta";
$map_a =~ s/\.match\.sam/\.match/;

my $row;
my $unmap_id; my $unmap_seq; my $unmap_qual;
my $map_id; my $map_seq; my $map_qual;
my $m = my $nm = my $mcount1 = my $mcount = my $nmcount1 = my $nmcount = 0;

# separating non-mapping sequences
open($INFILE, "<", $ARGV[0]) or die "$ARGV[0] reading error: ".$!;
if ( $opts{q} ) {
    open($OUT_UNM_Q, ">", $unmap_q) or die "$unmap_q writing error: ".$!;
    open($OUT_MAP_Q, ">", $map_q) or die "$map_q writing error: ".$!;
}
if ( $opts{a} ) {
    open($OUT_UNM_A, ">", $unmap_a) or die "$unmap_a writing error: ".$!;
    open($OUT_MAP_A, ">", $map_a) or die "$map_a writing error: ".$!;
}
while(1) {
    $row = readline $INFILE;
    last unless defined $row;
    $row =~ s/\n//; $row =~ s/\r//;
    my @tabs = split(/\t/, $row);
    if ($tabs[0] =~ m/^@/) {next;}
    my $flagcheck = $tabs[1];
    # sam file flag: 4 = "the query sequence itself is unmapped"
    if ( $flagcheck == 4 ) {
        $nm++;
        if ( $opts{a} ) {
            print $OUT_UNM_A "\>".$tabs[0]."\n";
            print $OUT_UNM_A $tabs[9]."\n";
        }
        if ( $opts{q} ) {
            print $OUT_UNM_Q "\@".$tabs[0]."\n";
            print $OUT_UNM_Q $tabs[9]."\n";
            print $OUT_UNM_Q "+".$tabs[0]."\n";
            print $OUT_UNM_Q $tabs[10]."\n";
        }
    } else {
        $m++;
        if ( $opts{a} ) {
            print $OUT_MAP_A "\>".$tabs[0]."\n";
            print $OUT_MAP_A $tabs[9]."\n";
        }
        if ( $opts{q} ) {
            print $OUT_MAP_Q "\@".$tabs[0]."\n";
            print $OUT_MAP_Q $tabs[9]."\n";
            print $OUT_MAP_Q "+".$tabs[0]."\n";
            print $OUT_MAP_Q $tabs[10]."\n";
        }
    }
    # progress meter
    #$mcount = $m;
    #$nmcount = $nm;
    #$mcount = int($mcount / 100000) * 100000;
    #$nmcount = int($nmcount / 100000) * 100000;
    #if ( ( $mcount != $mcount1 ) || ( $nmcount != $nmcount1 ) ) {
    #   print $mcount." matching reads, ".$nmcount." non-matching reads filtered...\n";
    #}
    #$mcount1 = $mcount;
    #$nmcount1 = $nmcount;
}
print $m." matching reads, ".$nm." non-matching reads filtered.\n";
close($INFILE);
if ( $opts{q} ) {
    close($OUT_MAP_Q);
    close($OUT_UNM_Q);
}
if ( $opts{a} ) {
    close($OUT_MAP_A);
    close($OUT_UNM_A);
}
