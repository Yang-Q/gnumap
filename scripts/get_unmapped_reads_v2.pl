#!/usr/bin/perl -w

use strict;
use warnings;
use Cwd;

sub check_fastq_pid_head_format {

	my ($fastq) = @_;
	my ($line, $pID_in_header,$sSize,$FQ);
	my @s;
	
	open($FQ,"<$fastq") || die "cannot open input fastq file, $fastq!\n";
	$line = readline($FQ);
	$pID_in_header=1;
	if ($line =~ m/^\S+\/(\d+)/) {
		$pID_in_header=1;
	}
	else {
		@s = split(' ',$line);
		$sSize = $#s+1;
		if ($sSize>1){
			$pID_in_header=0;
		}
	}
	close($FQ);
	return $pID_in_header;
}

my ($readDir,$inFmt, $outFmt,$orgFq1,$orgFq2,$filtReadDfn,$samDfile,$outDir) = @ARGV; #orgFq is _1.txt only and $samDfile is the one done by gnumap, which contains both pe12 and discordants.

#read samDfile and store readId to hash table
my %hMapped=();
my %hFiltReadId=();

my (@FqLine,@ORGFQ,@UNMAPFQ,@OrgRead,@FaDir)=();
my ($j,$seqId,$samLine,$fqLine,$pairLnk,$i,$i2,$pat,$fileBase,$pIdInHeader,$pair1,$pair2,$delimiterInHead,$filterF);

my $fafqFlag = '@';
if ($inFmt eq "fa"){
	$fafqFlag = '>';
}
else {
	$pIdInHeader = check_fastq_pid_head_format("$readDir/$orgFq1");
	$pair1="/1";
	$pair2="/2";
	$delimiterInHead = "";
	if ($pIdInHeader == 1) {
		$delimiterInHead = "\/1";
	}
}

my $fafqFlag2 = '@';
if ($outFmt eq "fa"){
	$fafqFlag2 = '>';
}

print STDOUT "selecting unmapped reads from an initial mapping in $outFmt format...\n";
if ($orgFq2 eq "X") {
	$i2=1;
	@OrgRead=($orgFq1);
}
else {
	$i2=2;
	@OrgRead=($orgFq1,$orgFq2);
}

my $prevMappedF = 0;
if ($samDfile eq "X") {
	$prevMappedF = 1;
}
else {
	open(SAMFP,"<$samDfile") || die "cannot read $samDfile !\n";
	#store all seqId sucessfully mapped into a hash
	while($samLine = <SAMFP>) {
		if ($samLine =~ /^(\S+)\t/) { #capture common readId since we know this is a valid paired end sam!
			$hMapped{$1}=1;
		}
	}
	close(SAMFP);
}

if ($filtReadDfn ne "X") {
	open(FILT,"<$filtReadDfn") || die "cannot read $filtReadDfn !\n";
	while($seqId = <FILT>) {
		$seqId = chomp($seqId);
		$hFiltReadId{$seqId}=1;
	}
	close(FILT);
}


for ($i=0;$i<$i2;$i++) { #open fq files
	open($ORGFQ[$i],"<$readDir/$OrgRead[$i]") || die "cannot read $readDir/$OrgRead[$i] !\n";
	$FaDir[$i]="$outDir";

	unless (-d $FaDir[$i]) {
		system("mkdir -p $FaDir[$i]");
	}

	$OrgRead[$i] =~ /^(\S+)\./;
	$fileBase = $1;
	
	open($UNMAPFQ[$i],">$FaDir[$i]/$fileBase.$outFmt") || die "cannot read $FaDir[$i]/$OrgRead[$i] !\n";
}

my $readUnitCnt;
my $fqCnt;

if ($i2==2) { #for PER
	seek($ORGFQ[0],0,0);
	seek($ORGFQ[1],0,0);
	
	$readUnitCnt = 4;
	$fqCnt = 4;
	if ($inFmt eq "fq") {
		$fqCnt = 8;
		$readUnitCnt = 8;
	}

	#read orgFq Fasta and extract reads that were not mapped anywhere else,
	while (defined($FqLine[0] = readline($ORGFQ[0])) && defined($FqLine[1] = readline($ORGFQ[1]))) {
	
		if (($fqCnt==$readUnitCnt) && ($FqLine[0] =~ /^$fafqFlag(\S+)$delimiterInHead/)) {
			$fqCnt=2;
			if ($prevMappedF || !(exists($hMapped{$1}))) {
			
				$filterF = 0;
				if (exists($hFiltReadId{$1})) {$filterF = 1;}
				
				unless ($filterF) {
					printf { $UNMAPFQ[0] } "${fafqFlag2}%s%s\n", ($1,$pair1);
					printf { $UNMAPFQ[1] } "${fafqFlag2}%s%s\n", ($1,$pair2);
				}
				
				for ($i=0;$i<2;$i++) { # for each pair read, print seq
					$FqLine[$i] = readline($ORGFQ[$i]);
					$fqCnt += 1;
					unless ($filterF) {
						printf { $UNMAPFQ[$i] } "%s", $FqLine[$i];
					}
				}

				if ($inFmt eq "fq") {
					#skip quality information
					for ($j=0;$j<2;$j++) {  #for each pair
						for ($i=0;$i<2;$i++) { #for each line
							$FqLine[$i] = readline($ORGFQ[$i]); #+
							$fqCnt += 1;
							if (($filterF==0) && ($outFmt eq "fq")) {
								printf { $UNMAPFQ[$i] } "%s", $FqLine[$i]; #qual
							}
						}
					}
				}
			}
		}
		else {
			$fqCnt+=2;
		}
	}
}
else {
	seek($ORGFQ[0],0,0);
	$readUnitCnt = 2;
	$fqCnt = 2;
	if ($inFmt eq "fq") {
		$fqCnt = 4;
		$readUnitCnt = 4;
	}
	#read orgFq Fasta and extract reads that were not mapped anywhere else,
	while (defined($FqLine[0] = readline($ORGFQ[0]))) {
		if (($fqCnt==$readUnitCnt) && ($FqLine[0] =~ /^$fafqFlag(\S+)/)) {
			$fqCnt=1;
			if ($prevMappedF || !(exists($hMapped{$1}))) {
			
				$filterF = 0;
				if (exists($hFiltReadId{$1})) {$filterF = 1;}
				
				$FqLine[0] = readline($ORGFQ[0]);
				$fqCnt += 1;
				unless ($filterF){
					printf { $UNMAPFQ[0] } "${fafqFlag2}%s\n", $1;
					printf { $UNMAPFQ[0] } "%s", $FqLine[0];
				}
				
				#skip quality information
				if ($inFmt eq "fq") {
					for ($j=0;$j<2;$j++) {
						$FqLine[0] = readline($ORGFQ[0]); #+
						$fqCnt += 1;
						if (($filterF==0) && ($outFmt eq "fq")) {
							printf { $UNMAPFQ[0] } "%s", $FqLine[0]; #qual
						}
					}
				}
			}
		}
		else {
			$fqCnt+=1;
		}
	}
}

for ($i=0;$i<$i2;$i++) {
	close($ORGFQ[$i]);
	printf { $UNMAPFQ[$i] } "\n"; #qual
	close($UNMAPFQ[$i]);
	#print STDERR "$FaDir[$i]/$OrgRead[$i]\n"; #debug
}
print STDOUT "done.\n"
#print STDERR "$samDfile\n"; #debug
