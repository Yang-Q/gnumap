#!/usr/bin/perl -w

# Author: cjhong
use strict;
use warnings;

my $inSam = $ARGV[0];
my $outSam = $ARGV[1];

my (@samFields);
my ($samLine,$outP);

my ($SAM_QUERY_STRAND2,$SAM_QUERY_STRAND)=(4,16);

open(INP,"<$inSam") || die "cannot read $inSam file!\n";
open($outP,">$outSam") || die "cannot read $outSam file!\n";

print STDERR "processing $inSam ...\n";

while ($samLine = <INP>){
	@samFields = split("\t",$samLine);
	
	if ($samFields[1] & $SAM_QUERY_STRAND) {
		
		$samFields[1] &= ~(1 << $SAM_QUERY_STRAND2);

		#$samFields[9] = reverse_complement($samFields[9]);
		$samLine = join("\t",@samFields);
	}
	printf { $outP } "%s", $samLine;
}
close(INP);
close($outP);

print STDERR "done.\n";

sub reverse_complement {
	my $dna = shift;
	# reverse the DNA sequence
	my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ACGTNacgtn/TGCANtgcan/;
	return $revcomp;
}
