#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd;

my ($readDir,$fastaF,$orgFq1,$orgFq2,$samDfile,$outDir) = @ARGV; #orgFq is _1.txt only and $samDfile is the one done by gnumap, which contains both pe12 and discordants.

#read samDfile and store readId to hash table
my %ReadIdMapped=();
my (@FqLine,@ORGFQ,@MAPFQ,@OrgRead,@FaDir)=();
my ($j,$seqId,$samLine,$pairLnk,$i,$i2,$pat,$fileBase);

my $fafqFlag = '@';
$fastaF = $fastaF+0;
if ($fastaF==1){
	$fafqFlag = '>';
}

print STDOUT "selecting unmapped fastq reads from an initial mapping...\n";
if ($orgFq2 eq "X") {
	$i2=1;
	@OrgRead=($orgFq1);
}
else {
	$i2=2;
	@OrgRead=($orgFq1,$orgFq2);
}

open(SAMFP,"<$samDfile") || die "cannot read $samDfile !\n";
#store all seqId sucessfully mapped into a hash
while($samLine = <SAMFP>) {
	if ($samLine =~ /^(\S+)\t/) { #capture common readId since we know this is a valid paired end sam!
		$ReadIdMapped{$1}=1;
	}
}
close(SAMFP);

for ($i=0;$i<$i2;$i++) { #open fq files
	open($ORGFQ[$i],"<$readDir/$OrgRead[$i]") || die "cannot read $readDir/$OrgRead[$i] !\n";
	$FaDir[$i]="$outDir";

	unless (-d $FaDir[$i]) {
		system("mkdir -p $FaDir[$i]");
	}

	$OrgRead[$i] =~ /^(\S+)\./;
	$fileBase = $1;
	
	open($MAPFQ[$i],">$FaDir[$i]/$fileBase.fq") || die "cannot read $FaDir[$i]/$OrgRead[$i] !\n";
}

if ($i2==2) { #for PER
	#explore qname format
	$samLine = readline($ORGFQ[0]);
	my @words=split('\t',$samLine);
	#print STDERR "$words[0]\n";
	my @pairLink=('\/','_');
	my @pairLink2=('/','_');
	
	for ($i=0;$i<2;$i++) {
		for ($j=1;$j<3;$j++) {
			$pat="$pairLink[$i]$j";
			#print STDERR "$pat\n";
			if ($words[0] =~ m/$pat$/) {
				$pairLnk=$pairLink2[$i];
				last;
			}
		}
	}

	my $pair1="${pairLnk}1";
	my $pair2="${pairLnk}2";
	seek($ORGFQ[0],0,0);

	#read orgFq Fasta and extract reads that were not mapped anywhere else,
	while (defined($FqLine[0] = readline($ORGFQ[0])) && defined($FqLine[1] = readline($ORGFQ[1]))) {
	
		if ($FqLine[0] =~ /^$fafqFlag(\S+)$pair1/) {
			
			if (exists($ReadIdMapped{$1})) {

				printf { $MAPFQ[0] } "\@%s%s\n", ($1,$pair1);
				printf { $MAPFQ[1] } "\@%s%s\n", ($1,$pair2);
				
				for ($i=0;$i<2;$i++) { # for each pair read, print seq
					for ($j=0;$j<3;$j++) {
						$FqLine[$j] = readline($ORGFQ[$i]);
						printf { $MAPFQ[$i] } "%s", $FqLine[$j];
					}
				}
			}
		}
	}
}
else {
	print STDERR "ser not supported yet...\n";
}

for ($i=0;$i<$i2;$i++) {
	close($ORGFQ[$i]);
	close($MAPFQ[$i]);
	#print STDERR "$FaDir[$i]/$OrgRead[$i]\n"; #debug
}
print STDOUT "done.\n"
#print STDERR "$samDfile\n"; #debug
