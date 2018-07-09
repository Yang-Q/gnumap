#!/usr/bin/perl -w

# Author: cjhong

# psl2sam_v3.pl read_fasta.fa blat.psl(x) blat.sam
# v6: works with blatSrc35.zip
# v7: control scoring rule. For example, it is useful to filter out any read containing false positves of a/g RNA editing events; support unstranded/strand1/strand2 lib
# v8: 
#-no introns allowed when 1)the read type is DNA or 2) RNA reads are aligned to transcriptome sequence
#-append a potential strand that is consistent with nt_conversion
#v9: stranded support for nt-conversion

#Note that Qstart begins from 0 bp whereas Qend from 1bp
# This script calculates a score using the BLAST scoring
# system. However, I am not sure how to count gap opens and gap
# extensions. It seems that column 5-8 are not what I am
# after. This script counts gaps from the last three columns. It does
# not generate reference skip (N) in the CIGAR as it is not easy to
# directly tell which gaps correspond to introns.

use strict;
use warnings;
use Getopt::Std;

#my %opts = (a=>1, b=>3, q=>5, r=>2);

#getopts('a:b:q:r:', \%opts);
#die("Usage: psl2sam.pl [-a $opts{a}] [-b $opts{b}] [-q $opts{q}] [-r $opts{r}] <in.psl>\n") if (@ARGV == 0 && -t STDIN);
our (@gNtConv)=();

die("Usage: psl2sam.pl <orgRead.fa> <in.psl> <out.sam> <query_type>\n") if (@ARGV == 0 && -t STDIN);

my $orgRead = $ARGV[0];
my $inPsl = $ARGV[1];
my $outSam = $ARGV[2];
my $queryType = $ARGV[3];
my $ntConv = $ARGV[4];
my $numMatchCutoff = $ARGV[5]+0.0;
my $qualOn = $ARGV[6]+0.0;
my $perId = $ARGV[7]+0.0; #0:pair_1.fq, 1:pair_2.fq
my $libType = $ARGV[8]+0.0; #0:unstranded, 1:wt1, 2:cr1

my $readLen = 150.0;
my $mapQconst = 255.0/log($readLen*3);
my ($ntConvF,$ntConvMapStrands);

$ntConv =~ s/^\s+|\s+$//g ;
($ntConvF,$ntConvMapStrands) = check_nt_conv_matrix($ntConv,$perId);

my ($qSt0,$tSt0,$score,$line,$validF,$idx,$clipLen2,$cigar,$blockSz0,$lq,$lt,$mmScore,$intronLen,$queryLen,$cutoffScore,$i,$keepGoF,$inPslSz);
#my @stack;
#my $last = '';
#my ($a, $b, $q, $r) = ($opts{a}, $opts{b}, $opts{q}, $opts{r});
my ($INS,$DEL) = (0,1);
my (@s,@t,@blockSz,@qSt,@tSt,@qStr,@tStr,@As);

my $MAX_LEN_DEL = 4;

print STDOUT "converting psl($inPsl) to sam...\n";
if (-e $inPsl) {
	$inPslSz = -s $inPsl;
	if ($inPslSz < 1){exit(0);}
}
else{exit(0);}	

my ($hReadPsl) = read_seq_id_from_psl($inPsl);
my ($ReadQual,$ReadSeq);
#read fasta file that BLAT originally generated psl(x) file and store it into hash
if ($qualOn == 1.0) {
	($ReadSeq,$ReadQual) = read_fastq_in_hash($orgRead,$hReadPsl);
}
else {
	$ReadSeq = read_fasta_in_hash($orgRead,$hReadPsl);
}

my @gap_open = (0,0);
my @gap_ext = (0,0);

my ($a,$b) = (3,3); # for DNA, (match,-mismatch), respectively
my @q = (4,4); #gap open (ins,del)
my @r = (2,2); #gap ext (ins,del)


my $intronAllowedF = 0;
if ($queryType eq "rna") { # for RNA, lower penalty for deletion
	$intronAllowedF = 1;
	@q = (4,4);
	@r = (2,0);
}

open(PSL,"<$inPsl") || die "cannot read input file, $inPsl !\n" ;
open(SAM,">$outSam") || die "cannot write output file $outSam !\n";


while ($line = <PSL>) {
  next unless ($line =~ /^\d/);
	$validF = pslSanity($line);
	unless ($validF) {next;}
	
  @t = split('\t',$line);
  $cigar = '';
# 	if ($libType>0) {
# 		$keepGoF = check_lib_type_cond($t[8],$libType,$perId);
# 		if ($keepGoF==0) {next;}
# 	}
  if ($t[8] eq '+') {
		@s[0..4] = ($t[9], 0, $t[13], $t[15]+1, 0);
	}
	else {
		@s[0..4] = ($t[9], 16, $t[13], $t[15]+1, 0);
		$t[12] = $t[10] - $t[11];
	}
  @s[6..10] = ('*', 0, 0, '*', '*');

	#psl indicates both query and subj same block size but it may have a different distance between two consecutive start bp. In such case, it means ins or del!
  @blockSz = split(',', $t[18]); #each bloc size
  @qSt = split(',', $t[19]); #qStarts
  @tSt = split(',', $t[20]); #tStarts
  
  (@qStr,@tStr)=();
  if ($ntConvF && ($t[1]>0)) { #if this mapping is considered to catch nt conversion like C2T or A2G, let us focus on the one that there is mismatch between query and subj
		@qStr = split(',', $t[21]); #qStrings
		@tStr = split(',', $t[22]); #tStrings
	}
	
  ($blockSz0, $qSt0, $tSt0) = ($blockSz[0], $qSt[0], $tSt[0]);

	if ($qSt0) {
		$cigar .= $qSt[0] . 'S' ; # 5'-end clipping
		#$s[3] -= $qSt[0]; # pos - qStart
		#if ($s[3]<0) {next;}
	}
	
	@gap_open = (0,0);
	@gap_ext = (0,0);
	@As = (0,0); #(score in mapping to Watson,score in mapping to Crick)
	$intronLen = 0;
	
  for ($i=1;$i<$t[17];$i++) { #for each block (note that if it contains a single block, it will skip.)
		$lq = $qSt[$i] - $qSt[$i-1] - $blockSz[$i-1];
		$lt = $tSt[$i] - $tSt[$i-1] - $blockSz[$i-1];
		if ($lq < $lt) { # del: the reference gap is longer

			$gap_ext[$DEL] += $lt - $lq - 1;
			$cigar .= ($qSt[$i] - $qSt0) . 'M';
			
			if ($ntConvF && ($t[1]>0)) {
				comp2Strs($qStr[$i-1],$tStr[$i-1],$a,$b,\@As);
			}
			++$gap_open[$DEL];

			if ((($lt - $lq) > $MAX_LEN_DEL) && ($intronAllowedF==1)) {
				$cigar .= ($lt - $lq) . 'N';
				$intronLen += ($lt - $lq);
			}
			else {
				$cigar .= ($lt - $lq) . 'D';
			}
			
			($blockSz0, $qSt0, $tSt0) = ($blockSz[$i], $qSt[$i], $tSt[$i]);
		} elsif ($lt < $lq) { # ins: the query gap is longer
			++$gap_open[$INS];
			$gap_ext[$INS] += $lq - $lt - 1;
			$cigar .= ($tSt[$i] - $tSt0) . 'M';
			
			if ($ntConvF && ($t[1]>0)) {
				comp2Strs($qStr[$i-1],$tStr[$i-1],$a,$b,\@As);
			}
			
			$cigar .= ($lq - $lt) . 'I';
			($blockSz0, $qSt0, $tSt0) = ($blockSz[$i], $qSt[$i], $tSt[$i]);
		}
  }

  #the following section is for handling the last block
  $lq = $t[12] - $qSt0;
	$lt = $t[16] - $tSt0;
	
	if ($lq < $lt) { # del: the reference gap is longer
		$gap_ext[$DEL] += $lt - $lq - 1;
		$cigar .= ($lq) . 'M';
		
		if ($ntConvF && ($t[1]>0)) {
			comp2Strs($qStr[$t[17]-1],$tStr[$t[17]-1],$a,$b,\@As);
		}
		++$gap_open[$DEL];

		if (($lt - $lq) > $MAX_LEN_DEL) {
			$cigar .= ($lt - $lq) . 'N';
			$intronLen += ($lt - $lq);
		}
		else {
			$cigar .= ($lt - $lq) . 'D';
		}
  } elsif ($lt < $lq) { # ins: the query gap is longer
		++$gap_open[$INS];
		$gap_ext[$INS] += $lq - $lt - 1;
		$cigar .= $lt  . 'M' ;
		if ($ntConvF && ($t[1]>0)) {
			comp2Strs($qStr[$t[17]-1],$tStr[$t[17]-1],$a,$b,\@As);
		}
  } else {
		$cigar .= $lq . 'M' ;
		if ($ntConvF && ($t[1]>0)) {
			comp2Strs($qStr[$t[17]-1],$tStr[$t[17]-1],$a,$b,\@As);
		}
  }

	$clipLen2 = $t[10] - $t[12]; #qSize - qEnd
  if ( $clipLen2 > 0) {
		$s[5] = $cigar . $clipLen2 . 'S' ;
	}
	else {
		$s[5] = $cigar;
	}

	#$score = $a * $t[0] - $b * $t[1] - $q[$INS] * $gap_open[$INS] - $r[$INS] * $gap_ext[$INS] - $q[$DEL] * $gap_open[$DEL] - $r[$DEL] * $gap_ext[$DEL];
	$idx = $hReadPsl->{$s[0]};
	$s[9] = $ReadSeq->[$idx];
	#$cutoffScore = (length($ReadSeq->[$idx]) - $intronLen) * $a * $numMatchCutoff;
	$cutoffScore = $a * $numMatchCutoff;

	if ($ntConvF && ($t[1]>0)) {
		
		if ($As[0]<$As[1]) {
			$mmScore = $As[1];
			$s[13] = "XS:Z:$ntConvMapStrands[2]";
		}
		elsif ($As[0]>$As[1]) {
			$mmScore = $As[0];
			$s[13] = "XS:Z:$ntConvMapStrands[1]";
		}
		else {
			$mmScore = $As[0];
			$s[13] = "XS:Z:$ntConvMapStrands[0]";
		}
		$score = $mmScore - $q[$INS] * $gap_open[$INS] - $r[$INS] * $gap_ext[$INS] - $q[$DEL] * $gap_open[$DEL] - $r[$DEL] * $gap_ext[$DEL];
		
	}
	else {
		$score = $a * $t[0] - $b * $t[1] - $q[$INS] * $gap_open[$INS] - $r[$INS] * $gap_ext[$INS] - $q[$DEL] * $gap_open[$DEL] - $r[$DEL] * $gap_ext[$DEL];
		$s[13] = "XS:Z:$ntConvMapStrands[0]";
	}

  $score = 1 if ($score < 1);
	
	if ($score<$cutoffScore) {next;}

	$s[4] = int($mapQconst*log($score)+0.5);
	
  $s[11] = "AS:i:$score";
	$s[12] = "PG:Z:blat";

	unless (exists $hReadPsl->{$s[0]}) {
		die "check if $orgRead is consistent with $inPsl.\n";
	}
	
	if ($t[8] eq '-') {
		$s[9] = reverse_complement($s[9]);
	}

	if ($qualOn) {
		if ($t[8] eq '-') {
			$s[10] = reverse($ReadQual->[$idx])
		}
		else {
			$s[10] = $ReadQual->[$idx];
		}
	}
	else {
		$s[10] = 'I' x length($s[9]);
	}
  print SAM join("\t", @s), "\n";
}

close(PSL);
close(SAM);
print STDOUT "done.\n";

sub pslSanity {
	my ($pslLine) = @_;
	my $validF = 0;
	
	if ($pslLine =~ /^\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t[\+-]\t\S+\t\d+\t\d+\t\d+\t\S+\t\d+\t\d+\t\d+\t\d+\t(\d+,)+\t(\d+,)+\t(\d+,)/) {
		$validF = 1;
	}
	return ($validF);
}

sub read_seq_id_from_psl {
	my ($pslDfn) = @_;
	my $line;
	my $pslFP;
	open($pslFP,"<$pslDfn") || die "cannot read psl file, $pslDfn!\n";
	my %hRead = ();
	my $cnt = 0;
	while($line = <$pslFP>) {
		if ($line =~ /^\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t\d+\t[\+-]\t(\S+)\t/) {
			unless (exists $hRead{$1}) {
				$hRead{$1} = $cnt++;
			}
		}
	}
	close($pslFP);
	return (\%hRead);
}

sub read_fasta_in_hash {
	my ($fastaDfn,$hReadPsl) = @_;
	
	open(BLAT_IN_FA,"<$fastaDfn") || die "cannot read fasta file that blat used!\n";
	my ($line,$printF,$readIdx,$cnt); $cnt=0;
	my @ReadSeq = ();
	
	#print STDOUT "storing fasta input read that BLAT used previously!\n";
	while($line = <BLAT_IN_FA>) {
		chomp($line);
		if ($line =~ /^>(.*)/) {
			$printF = 0;
			if (exists $hReadPsl->{$1}) {
				$printF = 1;
				$readIdx = $hReadPsl->{$1};
				$cnt++;
			}
		}
		elsif ($printF) {
			$ReadSeq[$readIdx] = $line;
		}
	}
	close(BLAT_IN_FA);
	#print STDOUT "$cnt reads are stored.\n";

	return (\@ReadSeq);
}

sub read_fastq_in_hash {
	my ($fastqDfn,$hReadPsl) = @_;
	
	open(BLAT_IN_FQ,"<$fastqDfn") || die "cannot read fasta file that blat used!\n";
	my ($line,$printF,$readIdx,$cnt); $cnt=0;
	my @ReadSeq = ();
	my @ReadQual = ();
	
	#print STDOUT "storing fasta input read that BLAT used previously!\n"; #debug
	while($line = <BLAT_IN_FQ>) {
		chomp($line);
		if ($line =~ /^@(.*)/) {
			$printF = 0;
			if (exists $hReadPsl->{$1}) {

				$readIdx = $hReadPsl->{$1};
				$line = <BLAT_IN_FQ>; # seq
				chomp($line);
				$ReadSeq[$readIdx] = $line;
				$line = <BLAT_IN_FQ>; # +
				$line = <BLAT_IN_FQ>; # qual
				chomp($line);
				$ReadQual[$readIdx] = $line;
				$cnt++;
			}
		}
	}
	close(BLAT_IN_FQ);
	#print STDOUT "$cnt reads are stored.\n"; #debug

	return (\@ReadSeq,\@ReadQual);
}

sub reverse_complement {
	my $dna = shift;
	# reverse the DNA sequence
	my $revcomp = reverse($dna);

	# complement the reversed DNA sequence
	$revcomp =~ tr/ACGTacgt/TGCAtgca/;
	return $revcomp;
}

sub comp2Strs {
	my ($qStr,$tStr,$matchSc,$mismatchSc,$nt_offset,$As) = @_;
	
	my ($i,$q,$t);
	#check if two string length is same!
	
	#compute aligment score based on a/g or t/c
	for ($i=0; $i<length($qStr);$i++) {
		$q = substr($qStr,$i,1);
		$t = substr($tStr,$i,1);
		#As->(0,0)=(score in mapping to Watson, score in mapping to Crick)
		if ($q eq $t) {
			$As->[0] += $matchSc;
			$As->[1] += $matchSc;
		} elsif (($t eq $gNtConv[0][0]) && ($q eq $gNtConv[0][1])) { #a/g?
			$As->[0] += ($matchSc+$nt_offset);
			$As->[1] -= $mismatchSc;
		} elsif (($t eq $gNtConv[1][0]) && ($q eq $gNtConv[1][1])) { #t/c?
			$As->[1] += ($matchSc+$nt_offset);
			$As->[0] -= $mismatchSc;
		} else {
			$As->[0] -= $mismatchSc;
			$As->[1] -= $mismatchSc;
		}
	}
}


sub check_nt_conv_matrix {

	my ($ntConv, $perId) = @_;
	
	my @ntConvMapStrands = ("nn","X","X");
	my $ntConvF = 0;
	my $nt_offset = 0;
	
	if ($ntConv eq 'a2i') {
		$ntConvF=1;
		if (($libType eq 'wt1' && $perId == 0) || ($libType eq 'cr1' && $perId == 1)) {
			@gNtConv=(['a','g'],['t','c']);
			@ntConvMapStrands=("nn","ag","tc");
		}
		else {
			@gNtConv=(['t','c'],['a','g']);
			@ntConvMapStrands=("nn","tc","ag");
		}
		$nt_offset = -0.75;
	}
	elsif ($ntConv eq 'bs') { #TODO
		@gNtConv=(['c','t'],['g','a']);
		$ntConvF=1;
		@ntConvMapStrands=("nn","ct","ga");
		$nt_offset = 0.75;
	}
	return ($ntConvF,\@ntConvMapStrands);
}

#----------------
#only support strand Seq
sub ntConvUpdScore {
	my ($qStr,$tStr,$matchSc,$mismatch,$nt_offset,$strandi,$sensi,$As) = @_;
	my ($i,$q,$t);
	#check if two string length is same!
	
	#compute aligment score based on a/g or t/c
	for ($i=0; $i<length($qStr);$i++) {
		$q = substr($qStr,$i,1);
		$t = substr($tStr,$i,1);
		
		#note that As[0] is a score sum for sense
		if ($q eq $t) {
			$As->[0] += $matchSc;
			#$As->[1] += $matchSc;
		} elsif (($t eq $gNtConv[$strandi][0]) && ($q eq $gNtConv[$strandi][1])) { #for stranded
			$As->[0] += ($matchSc+$nt_offset);
			#$As->[1] -= $mismatchSc;
		} elseif ($unstranded && ($t eq $gNtConv[($strandi+1)%2][0]) && ($q eq $gNtConv[($strandi+1)%2][1])) {
			$As->[1] += ($matchSc+$nt_offset);
			$As->[0] -= $mismatchSc;
		} else {
			$As->[0] -= $mismatchSc;
			#$As->[1] -= $mismatchSc;
		}
	}
	
}

sub needleman_wunch {
	# Needleman-Wunsch  Algorithm 
	# usage statement
	#die "usage: $0 <sequence 1> <sequence 2>\n" unless @ARGV == 2;

	# get sequences from command line
	my ($seq1, $seq2) = @ARGV;

	# scoring scheme
	my $MATCH    =  1; # +1 for letters that match
	my $MISMATCH = -1; # -1 for letters that mismatch
	my $GAP      = -1; # -1 for any gap

	# initialization
	my @matrix;
	$matrix[0][0]{score}   = 0;
	$matrix[0][0]{pointer} = "none";
	for(my $j = 1; $j <= length($seq1); $j++) {
			$matrix[0][$j]{score}   = $GAP * $j;
			$matrix[0][$j]{pointer} = "left";
	}
	for (my $i = 1; $i <= length($seq2); $i++) {
			$matrix[$i][0]{score}   = $GAP * $i;
			$matrix[$i][0]{pointer} = "up";
	}

	# fill
	for(my $i = 1; $i <= length($seq2); $i++) {
			for(my $j = 1; $j <= length($seq1); $j++) {
					my ($diagonal_score, $left_score, $up_score);

					# calculate match score
					my $letter1 = substr($seq1, $j-1, 1);
					my $letter2 = substr($seq2, $i-1, 1);                            
					if ($letter1 eq $letter2) {
							$diagonal_score = $matrix[$i-1][$j-1]{score} + $MATCH;
					}
					else {
							$diagonal_score = $matrix[$i-1][$j-1]{score} + $MISMATCH;
					}

					# calculate gap scores
					$up_score   = $matrix[$i-1][$j]{score} + $GAP;
					$left_score = $matrix[$i][$j-1]{score} + $GAP;

					# choose best score
					if ($diagonal_score >= $up_score) {
							if ($diagonal_score >= $left_score) {
									$matrix[$i][$j]{score}   = $diagonal_score;
									$matrix[$i][$j]{pointer} = "diagonal";
							}
					else {
									$matrix[$i][$j]{score}   = $left_score;
									$matrix[$i][$j]{pointer} = "left";
							}
					} else {
							if ($up_score >= $left_score) {
									$matrix[$i][$j]{score}   = $up_score;
									$matrix[$i][$j]{pointer} = "up";
							}
							else {
									$matrix[$i][$j]{score}   = $left_score;
									$matrix[$i][$j]{pointer} = "left";
							}
					}
			}
	}

	# trace-back

	my $align1 = "";
	my $align2 = "";

	# start at last cell of matrix
	my $j = length($seq1);
	my $i = length($seq2);

	while (1) {
			last if $matrix[$i][$j]{pointer} eq "none"; # ends at first cell of matrix

			if ($matrix[$i][$j]{pointer} eq "diagonal") {
					$align1 .= substr($seq1, $j-1, 1);
					$align2 .= substr($seq2, $i-1, 1);
					$i--;
					$j--;
			}
			elsif ($matrix[$i][$j]{pointer} eq "left") {
					$align1 .= substr($seq1, $j-1, 1);
					$align2 .= "-";
					$j--;
			}
			elsif ($matrix[$i][$j]{pointer} eq "up") {
					$align1 .= "-";
					$align2 .= substr($seq2, $i-1, 1);
					$i--;
			}    
	}

	$align1 = reverse $align1;
	$align2 = reverse $align2;
	#print STDOUT "$align1\n"; #debug
	#print STDOUT "$align2\n"; #debug
}
# 
#   0    1     2		3 		4			5			6			7			8			9					10		11		12	13				14			15			16				17				18							19					20
# match	mis- 	rep. 	N's	Q gap	Q gap	T gap	T gap	strand	Q        	Q   	Q[0] 	Q[1] T       	T   	  T  	  	T 		 	block			blockSizes 			qStarts	 		tStarts
#      	match	match	   	count	bases	count	bases	      	name     	size	start	end	name     	size	  start		end			count
# ---------------------------------------------------------------------------------------------------------------------------------------------------------------
# 45	    1	    0	   0	  0	    0	     1	  22	  -	  YourSeq	    101	   0	  46 chrIV	 17493793	7256464	7256532			2					23,23,					55,78,		7256464,7256509,
# 24      1     0    0    0     0      0    0     -  HWI-ST141     36   11    36 chrIII  13783700   65345  65370      1          25,             0,       65345,
# 35      0	    0	   0	  0	    0	     0	  0	    +	 HWI-ST141 	   69	  34	  69 chrIII	 13783700	7595644	7595679	    1	         35,	           34,	    7595644
#	35	    0	    0	   0	  0	    0	     0	  0	    -	 HWI-ST141_    69	   0	  35 chrIII	 13783700	7595644	7595679	    1	         35,						 34,			7595644

# sub check_lib_type_cond {
# 	my ($samStrand,$libType,$perId) = @_;
# 	my $keepGo =0;
# 	if ($libType==1){ #stranded RNA; 1st read corresponds to sense gene
# 		if (($perId==0 && $samStrand eq '+') || ($perId==1 && $samStrand eq '-')) {
# 			$keepGo=1;
# 		}
# 	}
# 	elsif ($libType==2){
# 		if (($perId==0 && $samStrand eq '-') || ($perId==1 && $samStrand eq '+')) {
# 			$keepGo=1;
# 		}
# 	}
# 	else {
# 		$keepGo=1;
# 	}
# 	return $keepGo;
# }