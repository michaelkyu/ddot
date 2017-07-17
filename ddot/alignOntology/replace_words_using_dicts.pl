#!/usr/bin/perl
use strict;

my %originalToNew = ();

for (my $i = 1; $i <= $#ARGV; $i++) {
    open IDS_TO_NAMES_FILE, $ARGV[$i] or die $!;
    while (<IDS_TO_NAMES_FILE>) {
	chomp($_);
	my ($original, $new);
	($original, $new) = split(/\t/, $_);
	#print "original: " . $original . " new: " . $new . "\n";
	$originalToNew{ $original } = $new;
    }
    close IDS_TO_NAMES_FILE;
}

open ONTOLOGY_FILE, $ARGV[0] or die $!;
while (<ONTOLOGY_FILE>) {
    chomp($_);
    my @vals = split(/\t/,$_);
    my $val;
    my $first = 0;
    foreach $val (@vals) {
	my $new_val;
	if (exists $originalToNew{$val} ) {
	    $new_val = $originalToNew{$val};
	} else {
	    $new_val = $val;
	}
	if ($first == 1) {
	    print "\t";
	} else {
	    $first = 1;
	}
	print "$new_val";
    }
    print "\n";
}
close ONTOLOGY_FILE;
