#!/usr/bin/perl -w

## par_blat.pl -h oligomer.bumc.bu.edu -p 1045 -d /media/sleepysilver/data/innocentive_id_org/db/taxon_fa2/5cats/twoBit -t U.2bit,V.2bit -r /home/hong/projects/epigenomics/experiments/reads/simulation/greyBsrSim/dwgsimBsr-hs19-wta_1.fq -o  /home/hong/projects/epigenomics/experiments/reads/simulation/greyBsrSim/dwgsimBsr-hs19-wta_1-UV.bl8 -p 8 -e blast8

use strict;
use warnings;
use Getopt::Long;
use threads;
use threads::shared;

our ($gDebug);

my ($numThreads,$minIdn,$tileSz,$stepSz,$indexDfn,$alnType,$readDfn,$repMatch,$maxIntron,$mcType,$alnDfn,$help,$ref,$pairId,$debug);

my %opts = ();

$pairId = 1;
$numThreads=1;
$minIdn=90;
$tileSz=12;
$stepSz=6;
$repMatch=256;
$mcType='dna';
$alnType='pslx';
$debug = 0;

my @usage;
push @usage, "  --help     Displays this information\n";
push @usage, "  --ref  a target reference genome file to be searched in 2 bit file format\n";
push @usage, "  --read_file  a read file to align in FASTA file format\n";
push @usage, "  --pair_id    is this first sequenced read or second? [1] for _1.fastq, 2 for _2.fastq\n";
push @usage, "  --out_file   specify an output file in full path\n";
push @usage, "  --num_threads  the number of threads to use\n";
push @usage, "  --out_fmt  bl7 or [pslx]\n";
push @usage, "  --min_identity  a minimum percent identity in in BLAT search\n";
push @usage, "  --tile_size   specify a tileSz in BLAT search\n";
push @usage, "  --step_size   specify a stepSz in BLAT search\n";
push @usage, "  --rep_match   specify the number of repMatch in BLAT search\n";
push @usage, "  --molecule    is this rna or dna\n";
push @usage, "  --max_intron  maximum intron to be allowed in BLAT search\n";
push @usage, "  --debug       run this program in debug mode?[0]\n";

GetOptions
(
 'help'        => \$help,
 'ref=s'    => \$ref, #template config file
 'read_file=s'  => \$readDfn,
 'pair_id=i' => \$pairId,
 'out_file=s'    => \$alnDfn,
 'num_threads=i' => \$numThreads,
 'out_fmt=s' => \$alnType,
 'min_identity=i' => \$minIdn,
 'tile_size=i' => \$tileSz,
 'step_size=i' => \$stepSz,
 'rep_match=i' => \$repMatch,
 'molecule=s' => \$mcType,
 'max_intron=i' => \$maxIntron,
 'debug=i' => \$debug
);

$gDebug = $debug;
not defined $help or die @usage;
defined $ref or die @usage;
defined $readDfn or die @usage;
defined $alnDfn or die @usage;

my $blatBin = "blat";

my ($alnD,$alnBase,$alnExt)=parse_filename($alnDfn);
unless (-d $alnD) { system("mkdir -p $alnD"); }

print STDERR "splitting reads...\n";
my ($splitInPrefix,$splitAlnPrefix)=split_fastq_reads($readDfn,$pairId,$alnD,$alnBase,$numThreads);


print STDERR "Performing BLAT search...\n";
my $options="-minIdentity=$minIdn -out=$alnType -noHead -q=$mcType -repMatch=$repMatch";

#TODO:
#blat_option($mcType,$oocFn,$tileSz,$stepSz,$oneOff,$minMatch,$minScore,$minIdn,$maxGap,$noHead,$repMatch,$maskType,$qMask,$repeats,$alnType,$maxIntron);

par_blat_in_bash($blatBin,$ref,$splitInPrefix,$splitAlnPrefix,$alnDfn,$options);
merge_split_aln($splitInPrefix, $splitAlnPrefix,$alnDfn);


sub merge_split_aln {
	my ($splitInPrefix, $splitAlnPrefix,$alnDfn) = @_;
	
	my $cmd = "cat $splitAlnPrefix\* > $alnDfn\n";
	if ($gDebug == 1){
		print STDERR $cmd;
	}
	system($cmd);
	
	if ($gDebug == 0){
		$cmd = "rm -rf $splitInPrefix\* \n";
		print STDERR $cmd;
		system($cmd); #debug
		$cmd = "rm -rf $splitAlnPrefix\* \n";
		system($cmd); #debug
		print STDERR $cmd;
	}
}

sub par_blat_in_bash {
	my ($blatBin,$ref,$splitInPrefix,$splitAlnPrefix,$alnDfn,$options) = @_;
	
	my $scriptDfn="$alnD/gnupar.sh";
	
	print "running the script $scriptDfn...\n";
	
	open (PARP,">$scriptDfn") || die "cannot create $scriptDfn\n";
	print PARP "#! /bin/bash -l\n";
	print PARP "$blatBin \"\$1\" $options \"\$2\" \"\$3\"\n";
	close(PARP);
	system("chmod +x $scriptDfn");
	
	#check if $ref has its 2bit format file
	my $ref2bit;
	my ($refD,$refBase,$refExt) = parse_filename($ref);
	
	if ($refExt eq "2bit") {
		$ref2bit = $ref;
	}
	else {
		$ref2bit = "$refD/$refBase.2bit";
		unless (-e $ref2bit) {
			my $cmd = "faToTwoBit $ref $ref2bit";
			system($cmd);
		}
	}
	#die "stopme";
	my $cmd = "parallel --xapply $scriptDfn $ref2bit ::: $splitInPrefix\* ::: $splitAlnPrefix\*";
	if ($gDebug==1) {
		print STDERR "$cmd\n";
	}
	system($cmd); #debug
	
	print "done.\n";
}

sub parse_filename {
	my ($filename) = @_; 
	my ( $fileD, $fileBase, $fileExt, $file ) =  ( 'x', 'x', '','x' );
	
	if ( $filename =~ /(\S*)\/(\S+)/ )
	{    #check if it contains directory information
		$fileD = $1;
		$file  = $2;
		if ( $file =~ /(\S+)\.(\S+)/ ) {    #check if it contains file extension
			$fileBase = $1;
			$fileExt  = $2;
		}
		else {
			$fileBase = $file;
		}
	}
	return ( $fileD, $fileBase, $fileExt );
}

sub check_fastq_pid_head_format {

	my ($fastq) = @_;
	my ($line, $pID_in_header,$sSize,$FQ,$read_fmt);
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
	
	#rewind
	seek $FQ, 0, 0;
	
	if ($line =~ m/^>/) {
		$read_fmt = 1;
	}
	elsif ($line =~ m/^@/) {
		$read_fmt = 2;
	}
	else {
		die "check input file [$fastq]";
	}
	
	close($FQ);
	return ($pID_in_header,$read_fmt);
}

sub readSlashPairIdFastq {
	my ($readFmt,$fqFile,$numSplits,$splitFaWps) = @_;
	open(INFQ,"<$fqFile") || die "cannot read $fqFile...";
	my ($j,$line);
	my $c=0;
	my $readMark=">";
	if ($readFmt == 2){$readMark="@";}
	while ($line = <INFQ>) {
		if ($line =~ m/^$readMark(\S+)/) {
			$j = $c % $numSplits;
			printf { $splitFaWps->[$j] } ">%s\n", $1; #head
			$line = <INFQ>;
			printf { $splitFaWps->[$j] } "%s", $line; #nt
			if ($readFmt == 2){
				$line = <INFQ>; #+
				$line = <INFQ>; #base_qual
			}
			$c=$c+1;
		}
	}
	close(INFQ);
}

sub readSpaceHeadFastq {
	my ($readFmt,$fqFile,$pairId,$numSplits,$splitFaWps) = @_;
	open(INFQ,"<$fqFile") || die "cannot read $fqFile...";
	my ($j,$line);
	my $c=0;
	my $readMark=">";
	if ($readFmt == 2){$readMark="@";}
	while ($line = <INFQ>) {
		if ($line =~ m/^$readMark(\S+)/) {
			$j = $c % $numSplits;
			printf { $splitFaWps->[$j] } ">%s/$pairId\n", $1; #head
			$line = <INFQ>;
			printf { $splitFaWps->[$j] } "%s", $line; #nt
			if ($readFmt == 2){
				$line = <INFQ>; #+
				$line = <INFQ>; #base_qual
			}
			$c=$c+1;
		}
	}
	close(INFQ);
}


#TODO: check if input read is fastq or fasta format
sub split_fastq_reads{
	my ($fqFile,$pairId,$alnD,$alnBase,$numThreads) = @_;
	
	my ($fp2,$c,$p,$readi,$alni,$line,$j,@readp,@alnp);
	
	my ($slashPairIdInHeadF,$readFmt)=check_fastq_pid_head_format($fqFile);
	my $splitInPrefix = "$fqFile.par_";
	my $splitAlnPrefix = "$alnD/$alnBase.par_";
	
	
	for ($p=0; $p<$numThreads; $p++){
		$readi="${splitInPrefix}${p}";
		open($readp[$p],">$readi") || die "cannot open $readi\n";
	}

	for ($p=0; $p<$numThreads; $p++){
		$alni="${splitAlnPrefix}${p}";
		open($alnp[$p],">$alni") || die "cannot open $alni\n";
	}
		
	print STDERR "splitting input reads...\n";
	if ($slashPairIdInHeadF == 1) {
		readSlashPairIdFastq($readFmt,$fqFile,$numThreads,\@readp);
	}
	else {
		readSpaceHeadFastq($readFmt,$fqFile,$pairId,$numThreads,\@readp);
	}
	print STDERR "done!\n";
		
	for ($p=0;$p<$numThreads;$p++){
		close($readp[$p]);
		close($alnp[$p]);
	}
	close(INFQ);
	
	return ($splitInPrefix,$splitAlnPrefix);
}
