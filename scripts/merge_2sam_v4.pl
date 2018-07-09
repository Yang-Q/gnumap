#!/usr/bin/perl -w
# postprocessing script to merge two sam files that aligned to +/- ref queryStrand, which considers the requirement of lib prep in paired-end sequencing
# TODO 
#1) collecting all possible range of PER, then select the pair whose the sum of the alignment scores is the best!
#2) Not only /1 and /2 format paired tag, to support all general per id tag.
#3) Provide two options for PER distance between two reads (fragment size method and mate inner distance)

# INCOMPELTE!
use File::Copy;
#use Statistics::Distributions;
use strict;
use warnings;
use Cwd;

#--------------#sam flag bit------------------------------
our ($SAM_SEQ_PAIRED,$SAM_VALID_PAIR,$SAM_QUERY_UNMAP,$SAM_MATE_UNMAPPED,$SAM_QUERY_STRAND,$SAM_STRAND_MATE,$SAM_1ST_IN_PAIR,$SAM_2ND_IN_PAIR,@SAM_IN_PAIR,$SAM_NOT_PRIMARY,@SAM_MAP2STRAND,@gMAP2STRAND,@gStrandOri,$VALID_PAIR,$INVALID_DIST,$MULTIPLE_PAIRS,$ONE_MATE_UNMAPPED,@gInsLib,$gPOS_STRAND,$gNEG_STRAND,$gMap1stR2up,$gMinTLen1,$gMapper,$gInclDiscordantF);

$SAM_SEQ_PAIRED=0;$SAM_VALID_PAIR=1;$SAM_QUERY_UNMAP=2;$SAM_MATE_UNMAPPED=3;$SAM_QUERY_STRAND=4;$SAM_STRAND_MATE=5;$SAM_1ST_IN_PAIR=6;$SAM_2ND_IN_PAIR=7;$SAM_NOT_PRIMARY=8;
$VALID_PAIR=1; $INVALID_DIST=0; $MULTIPLE_PAIRS=2; $ONE_MATE_UNMAPPED=3;
$gPOS_STRAND=0; $gNEG_STRAND=1;
#--------------#sam flag bit------------------------------

$|=1;
if(@ARGV < 10) {
	print STDERR 
	"
Usage: merge2sam.pl <input_sam_1> <input_sam_2> <paired_merge_out> <refType> <pair2ref> <scratchDirectory> <howToHandleIniputSamAfter>
	Where: <input_sam_?> is a sam file derived a pair input read file.
				<paired_merge> is the name of the file to be written that will form paired reads where stays within a valid fragment length.
				<refType> can be either 'dna' or 'rna'. Tell you this read was mapped to DNA ref sequence or Transcriptome sequence where splice junction was removed.
				<pair2ref> is 1 or 2. for example, 1 if <input_sam_1> aligns to +tive strand and <input_sam_2> aligns to -tive strand. Otherwise, set to 2.
				<scratchDirectory> is a full name including directory/path to use a space to sort
				<howToHandleIniputSamAfter> choose one among 'delete', 'zip', 'keep'.
				
				Note that depending on the lib prep, user may modify gInsLib[1] in line 51.
	INPUT:
	-----

	OUTPUT: paired_merge.txt and missing_merge.txt
	------
	";
	exit(1);
}
my $infile1 = $ARGV[0]; #_1.txt
my $infile2 = $ARGV[1]; #_2.txt
my $outfile1 = $ARGV[2]; #valid PE sam output
my $readType = $ARGV[3]; #dna or rna
my $map1stR2up = $ARGV[4]; #1 indicates _1 aligned to up and _2 aligned to dn. Otherwise, use 2!
#my $mapper = $ARGV[5];
my $per_dist = $ARGV[5]+0;
my $max_frag_len = $ARGV[6]+0;
my $includeDiscordant = $ARGV[7]+0;
my $scratchDir2sort = $ARGV[8];
my $handleInputFile = $ARGV[9];


my ($x,$GNU_SAM_X0,$GNU_SAM_XA,$GNU_SAM_XP,$GNU_SAM_XZ,$line,$pairID,$tmp_file,$outSamFp,$seqL,$seqId2upd);
my (@new_denom,@pos_denom,@neg_denom,@pos_denom_updF,@neg_denom_updF,@bufCnt,@queryStrand,@mapLib,@chrom,@startPos,@readLen,@readLenF,@words);

#$gMinTLen1 = 40; #it suggested from trim_galore (default 35)

# my ($GNUMAP,$BLAT) = (0,1);
# $gMapper = $GNUMAP;
# if ($mapper eq "blat") {
# 	$gMapper = $BLAT;
# }
# print STDOUT "validating PER ($infile1,$infile2) mappings...\n";

chomp($map1stR2up);
$gMap1stR2up = $map1stR2up+0;
$gInclDiscordantF = $includeDiscordant+0;
$GNU_SAM_XA = 11; #column index of raw_align_score
$GNU_SAM_XP = 12; #posterior prob.
#my $GNU_SAM_X0 = 13; #num of mappings
#my $GNU_SAM_XZ = 14; #gen or tran

#this tag will not be used for now. Todo
if ($readType eq "dna") {
	@gInsLib=($per_dist,$max_frag_len); #from hci rnaseq bioanalysis fragment length
}
elsif ($readType eq "rna") {
	@gInsLib=($per_dist,50000); #blat_pslx from gnumap_serial_v5.pl is consistent1
}
else {
	die "specify read type, dna or rna?\n";
}

my @noInput = (0,0);
unless (-e $infile1) {
	$noInput[0]=1;
}
# 
unless (-e $infile2) {
	$noInput[1]=1;
}

my $rndTag = int(rand(100000));

$tmp_file = $infile1 . "-temp1-$rndTag";
if ($noInput[0]==1 && $noInput[1]==1) {
	exit(0)
}

my $cmd = "cat ";
if ($noInput[0]==0) {
	$cmd = "$cmd $infile1 ";
}
if ($noInput[1]==0) {
	$cmd = "$cmd $infile2";
}

$x = system("$cmd | sort -T $scratchDir2sort -n -k1,3 -o $tmp_file");

if ( $? == -1 ){
  die "command failed: $!\n";
}

my $fileSz = -s $tmp_file;
if (!$fileSz){
	#print STDOUT "skipping forming PER...\n";
	exit(0)
}

open(INFILE, $tmp_file); #this is the merged file!
open($outSamFp, ">$outfile1") or die "\nERROR: in script merge2sam_v4.pl: Cannot open file '$outfile1' for writing\n";

@pos_denom = (0,0); @neg_denom = (0,0); @pos_denom_updF = (0,0); @neg_denom_updF = (0,0); 
$seqId2upd = 'X';
@bufCnt = (0,0); @readLenF=(0,0); @readLen=(0,0);

@SAM_IN_PAIR = ($SAM_1ST_IN_PAIR,$SAM_2ND_IN_PAIR);

if ($gMap1stR2up==1) {
	@SAM_MAP2STRAND = ($SAM_STRAND_MATE,$SAM_QUERY_STRAND);
	@gMAP2STRAND = ($gPOS_STRAND,$gNEG_STRAND);
	@gStrandOri = (1,-1);
}
else {
	@SAM_MAP2STRAND = ($SAM_QUERY_STRAND,$SAM_STRAND_MATE);
	@gMAP2STRAND = ($gNEG_STRAND,$gPOS_STRAND);
	@gStrandOri = (-1,1);
}

while ($line = <INFILE>) {
	#print "$line\n"; #debug
	@words=split('\t',$line);
	$seqL=length($words[9]);

	if ($words[0] =~ m/^$seqId2upd\/(\d+)/) { #if it stays in the same seqID block?, then, keep filling up
		$pairID = $1 - 1;
		unless ($readLenF[$pairID]){
			$readLen[$pairID]=$seqL;
			$readLenF[$pairID]=1;
		}

		fill_array2print($line,$pairID,\@mapLib,\@queryStrand,\@bufCnt
			,\@chrom,\@startPos,\@pos_denom,\@neg_denom,\@pos_denom_updF,\@neg_denom_updF);
	}
	else { #found a new entry of seqId2upd, so we need clean up buffer we filled up so far
		$new_denom[0] = $pos_denom[0] + $neg_denom[0];
		$new_denom[1] = $pos_denom[1] + $neg_denom[1];
		unless ($seqId2upd eq 'X') {
			update_sam_print(\@bufCnt,\@new_denom,$outSamFp,\@mapLib,\@queryStrand,\@chrom,\@startPos,\@readLen);
		}
		
		#okay, now let us start to fill a new block to a buffer
		$seqId2upd = 'X';
		if ($words[0] =~ m/^(\S+)\/1/) {
			$seqId2upd = $1;
			$pairID = 0;
			@readLenF=(1,0);
			$readLen[0]=$seqL;
		}
		elsif ($words[0] =~ m/^(\S+)\/2/) {
			$seqId2upd = $1;
			$pairID = 1;
			@readLenF=(0,1);
			$readLen[1]=$seqL;
		}

		@bufCnt = (0,0);
		@pos_denom = (0,0); @neg_denom = (0,0); @pos_denom_updF = (0,0); @neg_denom_updF = (0,0);

		fill_array2print($line,$pairID,\@mapLib,\@queryStrand,\@bufCnt
			,\@chrom,\@startPos,\@pos_denom,\@neg_denom,\@pos_denom_updF,\@neg_denom_updF);
	}
}

#finish up buffer in the last block
$new_denom[0] = $pos_denom[0] + $neg_denom[0];
$new_denom[1] = $pos_denom[1] + $neg_denom[1];
unless ($seqId2upd eq 'X') {
	update_sam_print(\@bufCnt,\@new_denom,$outSamFp,\@mapLib,\@queryStrand,\@chrom,\@startPos,\@readLen);
}

close($outSamFp);
close(INFILE);

$x = `rm -rf $tmp_file`;

if ($handleInputFile eq "delete") {
	$x = `rm -rf $infile1 $infile2`;
}
elsif ($handleInputFile eq "zip") {
	system("gzip -f $infile1 $infile2");
}

#print STDOUT "done.\n";

##===============================================================================================
sub fill_array2print {
	my ($line,$pairID,$mapLib,$queryStrand,$bufCnt,
		$chrom,$startPos,$pos_denom,$neg_denom,$pos_denom_updF,$neg_denom_updF) = @_;
	
	#$sam_buffer->[$pairID][$bufCnt->[$pairID]] = $line; #cmmon~, I didn't check the duplicated copy here and I assume that all input sam files doesn't so. If you are not sure, we recommend you use picard to remove all duplicated entries before running this script
# 	if ($gMapper == $GNUMAP) {
# 		if ($line =~ /^\S+\t\d+\t(\S+)\t(\d+)\t.*\tXA:f:(\S+)\t/) {
# 			$chrom->[$pairID][$bufCnt->[$pairID]] = $1;
# 			$startPos->[$pairID][$bufCnt->[$pairID]] = $2;
# 			#$alignScore->[$pairID][$bufCnt->[$pairID]] = $3 + 0;
# 		}
# 	}
# 	else {
# 		if ($line =~ /^\S+\t\d+\t(\S+)\t(\d+)\t.*\tAS:i:(\d+)\t/) {
# 			$chrom->[$pairID][$bufCnt->[$pairID]] = $1;
# 			$startPos->[$pairID][$bufCnt->[$pairID]] = $2;
# 			#$alignScore->[$pairID][$bufCnt->[$pairID]] = $3 + 0;
# 		}
# 	}

	if ($line =~ /^\S+\t\d+\t(\S+)\t(\d+)\t/) {
		$chrom->[$pairID][$bufCnt->[$pairID]] = $1;
		$startPos->[$pairID][$bufCnt->[$pairID]] = $2;
	}
	
	if ($line =~ /^\S+\t[0|256]+\t/) {
		if (!($pos_denom_updF->[$pairID])) {
			$pos_denom_updF->[$pairID] = 1;
		}
		$mapLib->[$pairID][$bufCnt->[$pairID]] = $line;
		$queryStrand->[$pairID][$bufCnt->[$pairID]] = $gPOS_STRAND;
	}
	elsif ($line =~ /^\S+\t[16|272]+\t/) {
		if (!($neg_denom_updF->[$pairID])) {
			$neg_denom_updF->[$pairID] = 1;
		}
		$mapLib->[$pairID][$bufCnt->[$pairID]] = $line;
		$queryStrand->[$pairID][$bufCnt->[$pairID]] = $gNEG_STRAND;
	}

	$bufCnt->[$pairID]++;
}


##=============================================================================================
sub check_pair_orient {
	my ($pos1,$polar1,$len1,$pos2,$polar2,$len2) = @_;
	
	my ($from1,$from2,$to1,$to2,$fromDist,$toDist);
	
	if ($polar1==1) {
		$fromDist = abs($pos1+$len1-1-$pos2);
		$toDist = abs($pos1-$pos2-$len2+1);
	}
	else {
		$fromDist = abs($pos1-$pos2-$len2+1);
		$toDist = abs($pos1+$len1-1-$pos2);
	}
	
	my $discordantF = 0;
	if ($fromDist<$toDist) {
		$discordantF = 1;
	}
	
	return $discordantF;
}
##=============================================================================================
sub build_paired_end_core2 {
	my ($mPos,$polarity,$readL,$mapLib) = @_;
	my (@mapLen,@peMapLoc_sorted,@peMapLoc,@samRep,@orgIdx,@a,@b,@dummy);
	#my ($cnt,$c1,$c1m,$r1,$absDistF,$which,$i,$j,$mapStatus,$diffAveInnerDist);
	my ($cnt,$absDistF,$which,$i,$j,$mapStatus,$diffAveInnerDist,$peProb,$discordantF);
	my $distTooLong = 0;
	my $peDist;	
	$mapLen[0] = $#{$mPos->[0]}+1;
	$mapLen[1] = $#{$mPos->[1]}+1;

	$cnt=0;
	#$r1=0;
	#$c1m = 0;
	#compute distance between two paired maps
	if ($mapLen[0]>0 && $mapLen[1]>0){
		for ($i=0;$i<$mapLen[0];$i++){ #for each entry of 1st pair mapping
			#$c1=0;
			for ($j=0;$j<$mapLen[1];$j++){ #for each entry of 2nd pair mapping
				$peDist = abs($mPos->[1][$j] - $mPos->[0][$i]); #here, $peDist = st1bp - st2bp inestead of ed1bp - st2bp
				if ((($polarity->[0][$i] + $polarity->[1][$j]) == 1) && (($peDist + $readL->[0]) <= $gInsLib[1])) {
					#check pair reads mapped inside
					$discordantF = check_pair_orient($mPos->[0][$i],$polarity->[0][$i],$readL->[0],$mPos->[1][$j],$polarity->[1][$j],$readL->[1]);
					
					if ($discordantF==1) {
 						last;
					}
					#mapId, distance-mdDist, loc1, loc2
					
					#make it close to 0 so that we can pick up the pairs whose distance is closed to $gInsLib[0] after sorting in descending order
					if (($gInsLib[0]<0) && ($readL->[0]<=(-1*$gInsLib[0]))) {
						$diffAveInnerDist = $peDist;
					}
					else {
						$diffAveInnerDist = abs($peDist + $readL->[0] - $gInsLib[0]);
					}
					$peMapLoc[$cnt]="$cnt\t$diffAveInnerDist\t$mPos->[0][$i]\t$mPos->[1][$j]\t$i\t$j";
					$cnt++;
				}
			}
		}

		if ($cnt<1) {
			$distTooLong = 1;
			for ($i=0;$i<$mapLen[0];$i++){ #for each entry of 1st pair mapping
				#$c1=0;
				for ($j=0;$j<$mapLen[1];$j++){ #for each entry of 2nd pair mapping
					$peDist = abs($mPos->[1][$j] - $mPos->[0][$i]);
					if (($polarity->[0][$i] + $polarity->[1][$j]) == 1) {
						#mapId, distance-mdDist, loc1, loc2
						$diffAveInnerDist=abs($peDist - $gInsLib[0]); #Note that $gInsLib[0] = averageFragLengthFromDist -primerL -adapterL
						#$sumAscore = $ascore->[0][$i] + $ascore->[1][$j];
						#$peMapLoc[$cnt]="$cnt\t$diffAveInnerDist\t$sumAscore\t$mPos->[0][$i]\t$mPos->[1][$j]";
						$peMapLoc[$cnt]="$cnt\t$diffAveInnerDist\t$mPos->[0][$i]\t$mPos->[1][$j]\t$i\t$j";
						$cnt++;
					}
				}
			}
		}
		
		#do one-to-one mapping following lib-quality-report
		#=> cnt, r1, c1, peMapLoc
		#make one dim array and - peak_dist and sorted in an ascending order
		@peMapLoc_sorted = sort { (split '\t', $a)[1] <=> (split '\t', $b)[1] } @peMapLoc; #mimic sortrows in matlab
		
		@peMapLoc_sorted = select_only_shortest_perDist($cnt,@peMapLoc_sorted);
		$cnt = scalar(@peMapLoc_sorted);
		
		@peMapLoc = str2mat('\t',$cnt,@peMapLoc_sorted);
		
		if ($distTooLong != 1) {
			#do one-to-one mapping following lib-quality-report
			($mapStatus,@samRep)= one2one_map($mapLen[0],$mapLen[1],$cnt,$readL,\@peMapLoc,$mapLib);
		}
		else  {
			set_sam_bit_flag($INVALID_DIST,$INVALID_DIST,-1,\@samRep,\@peMapLoc,$cnt,$mapLib,$mPos,$readL);
		}
	}
	else {
		$which=1;
		$cnt=$mapLen[1];
		if ($mapLen[0]>0){
			$which=0;
			$cnt=$mapLen[0];
		}
		@dummy = ((0,0),(0,0));
		set_sam_bit_flag($ONE_MATE_UNMAPPED,$which,0,\@samRep,\@dummy,0,$mapLib,$mPos,$readL);
	}

#debug===============================
# 	$j = $#samRep+1;
# 	for ($i=0;$i<$j;$i++){
# 		print STDOUT join(", ",@{$samRep[$i]});
# 		print STDOUT "\n";
# 	}
# 	print STDOUT "done.\n";
#debug===============================
	return @samRep
}

##===============================================================================================
sub update_sam_print {
	my ($bufCnt,$new_denom,$outSamFp,$mapLib,$queryStrand,$chrom,$startPos,$readLen)= @_;
	my (@samRep,@mPos,@polarity,@readL,@entr)=();
	my ($j, $i);
	
	#create two arrays, polarity and pn, as inputs to build_paired_end_core();
	for ($j=0; $j<2; $j++) { #per pairID
		for ($i=0;$i<$bufCnt->[$j];$i++) { #per mappedLocation
			$mPos[$j][$i]=$startPos->[$j][$i];
			$polarity[$j][$i]= $queryStrand->[$j][$i]#info from queryStrand at line 185!
			#$ascore[$j][$i] = $alignScore->[$j][$i];
		}
		$readL[$j]=$readLen->[$j];
	}
	@samRep = build_paired_end_core2(\@mPos,\@polarity,\@readL,\@mapLib); #todo: make sure that samRep contains rnext, pnext, tlen(|--------->
	#																																																						 <----------|)
	#																																																		 --->|              |<---

	#combine samRep and sam_buffer to finalize
	$j = $#samRep+1;

	for ($i=0; $i<$j; $i++) {
		@entr=split('\t',$samRep[$i][2]);
		$entr[0] =~ m/^(\S+)\/\d+/;
		$entr[0] = $1;
		$entr[1]=$samRep[$i][1];
		$entr[3]=$samRep[$i][0];
		$entr[7]=$samRep[$i][3];
		
		if ($entr[7]>0) {
			$entr[6]='=';
		}
		else {
			$entr[6]='*';
		}
		$entr[8]=$samRep[$i][4];
		
		#for two possiblity to finally report as potential mapping
		#1. concordant and the template length > minTempL2
		#2. mate unmapped
		if ((($entr[1]&(1<<$SAM_VALID_PAIR))!=0) || (($gInclDiscordantF) && (($entr[1]&(1<<$SAM_MATE_UNMAPPED))!=0))) {
			print $outSamFp join("\t",@entr);
		}
	}
}

sub set_sam_bit_flag {
	my ($mapStatus,$which,$dimS,$samRep,$peMapLoc,$N,$mapLib,$mPos,$readL) = @_;
	my ($i,$j,$I,$k,$cnt);
	my @I2=(0,0);

	if ($mapStatus==$ONE_MATE_UNMAPPED) {
		
		$I = $#{$mPos->[$which]}+1;
		
		for ($i=0;$i<$I;$i++) {
			$samRep->[$i][5] = $readL->[$which];
			$samRep->[$i][3] = 0;
			$samRep->[$i][4] = 0;
			$samRep->[$i][0] = $mPos->[$which][$i]; 
			$samRep->[$i][1] |= (1 << $SAM_SEQ_PAIRED);
			$samRep->[$i][1] |= (1 << $SAM_MATE_UNMAPPED);
			
			if ($gMAP2STRAND[$which]==$gNEG_STRAND){
				$samRep->[$i][1] |= (1 << $SAM_QUERY_STRAND);
			}
			$samRep->[$i][1] |= (1 << $SAM_IN_PAIR[$which]);
			$samRep->[$i][2] = $mapLib->[$which][$i];

			if ($I>1) {
				$samRep->[$i][1] |= (1 << $SAM_NOT_PRIMARY);
			}
		}
	}
	elsif ($mapStatus == $INVALID_DIST) {
		$cnt=0;
		for ($i=0;$i<$N;$i++){
			for ($j=0;$j<2;$j++) {
				$k=$peMapLoc->[$i][$j+4];
				$samRep->[$cnt][3] = 0;
				$samRep->[$cnt][4] = 0;
				$samRep->[$cnt][0] = $mPos->[$j][$k]; #query name
				$samRep->[$cnt][1] |= (1 << $SAM_SEQ_PAIRED);
				$samRep->[$cnt][1] |= (1 << $SAM_MAP2STRAND[$j]);
				$samRep->[$cnt][1] |= (1 << $SAM_IN_PAIR[$j]);
				$samRep->[$cnt][2] = $mapLib->[$j][$k]; 
				$samRep->[$cnt][5] = $readL->[$j];
				
				if ($N>1) {
					$samRep->[$cnt][1] |= (1 << $SAM_NOT_PRIMARY);
				}
				$cnt++;
			}
		}
	}
	elsif ($mapStatus == $MULTIPLE_PAIRS){
		for ($i=0;$i<$dimS;$i++) {
			$samRep->[$i][1] |= (1 << $SAM_NOT_PRIMARY);
			##$samRep->[$i][4] != #pnext;
		}
	}
}

sub checkValidInsert{
	my ($dist,$maxIns) = @_;
	my $validF=1;
	if ($dist>$maxIns) {
		$validF=0;
	}
	return $validF;
}

sub checkOne2oneAvail2{
	my ($r,$c,$PairList) = @_;
	my $vacancy = 1;
	if (($PairList->[0][$r]>0) || ($PairList->[1][$c]>0)){
		$vacancy=0;
	}
	return ($vacancy,$r,$c);
}

sub ind2sub_ele { #column-wise subindex
	my ($idx1,$numR) = @_;
	my ($rem,$quo);
	#todo
	if ($numR>0) {
		$rem = $idx1 % $numR;
		$quo = ($idx1-$rem) / $numR;
	}
	else {
		die "divided by zero?\n";
	}
	return($rem,$quo);
}

sub str2mat {
	my ($delim,$R,@narrowStripStr) = @_;
	my ($r,$c);
	my @mat=();
	
	for ($r=0;$r<$R;$r++) {
		@{$mat[$r]} = split($delim,$narrowStripStr[$r]);
	}
	return @mat;
}

sub select_only_shortest_perDist {
	my ($N,@peMapLoc)=@_;
	my $i;
	my @words = ();
	
	@words=split(' ',$peMapLoc[0]);
	my $shortest = $words[1]+0;
	
	for ($i=0;$i<$N;$i++) {
		@words=split(' ',$peMapLoc[$i]);
		$words[1] = $words[1]+0;
		if ($words[1]>$shortest){
			last;
		}
	}
	my @peMapLocTrunc = splice(@peMapLoc,0,$i);
	return @peMapLocTrunc;
}


sub one2one_map {
	my ($R,$C,$N,$readL,$peMapLoc,$mapLib)=@_; #R and C indicates a dim of the original pairs
	
	##my $R = $_[0]; my $C = $_[1]; my @peMapLoc = @{$_[2]}; my @gInsLib = @{$_[3]};
	#my $N = $R * $C;
	my (@pairList,@samRep)=();
	my ($cnt,$i,$j,$r,$c,$mapStatus,$vacancy,$validF,$fragLen,$key,$prevFlip);

	for ($i=0;$i<$R;$i++) {$pairList[0][$i]=0;}
	for ($j=0;$j<$C;$j++) {$pairList[1][$j]=0;}

	#this time looking for only good-looking couples
	$cnt = 0;
	my $multiPairF = 0;
	my %hPairLoc=();
	
	for ($i=0;$i<$N;$i++) {
		$validF=checkValidInsert($peMapLoc->[$i][1]+$gInsLib[0],$gInsLib[1]);
		if ($validF) { # 1 indicates valid distance
			($vacancy,$r,$c)=checkOne2oneAvail2($peMapLoc->[$i][4],$peMapLoc->[$i][5],\@pairList); # 0 indicates mapId
			if ($vacancy) {
				$pairList[0][$r]=1; #mark 1st pair as found!
				$pairList[1][$c]=1; #mark 2nd pair as found!

				$samRep[$cnt][0] = $peMapLoc->[$i][2];
				$samRep[$cnt][1] |= (1 << $SAM_SEQ_PAIRED);
				$samRep[$cnt][1] |= (1 << $SAM_VALID_PAIR);
				$samRep[$cnt][1] |= (1 << $SAM_MAP2STRAND[0]);
				$samRep[$cnt][1] |= (1 << $SAM_1ST_IN_PAIR);
				$samRep[$cnt][2] = $mapLib->[0][$peMapLoc->[$i][4]];
				$samRep[$cnt][3] = $peMapLoc->[$i][3]; # 3 indicates mapping position of 2nd read

				$fragLen = abs($peMapLoc->[$i][3] - $peMapLoc->[$i][2]) + $readL->[$gMAP2STRAND[1]];
				$samRep[$cnt++][4]= $gStrandOri[0] * $fragLen;

				$samRep[$cnt][0] = $peMapLoc->[$i][3];
				$samRep[$cnt][1] |= (1 << $SAM_SEQ_PAIRED);
				$samRep[$cnt][1] |= (1 << $SAM_VALID_PAIR);
				$samRep[$cnt][1] |= (1 << $SAM_MAP2STRAND[1]);
				$samRep[$cnt][1] |= (1 << $SAM_2ND_IN_PAIR);
				$samRep[$cnt][2] = $mapLib->[1][$peMapLoc->[$i][5]];
				$samRep[$cnt][3] = $peMapLoc->[$i][2]; # 2 indicates mapping position of 1st read
				$samRep[$cnt++][4]= $gStrandOri[1] * $fragLen;
			}
		}
	}

	$mapStatus = $VALID_PAIR;
	if ($cnt==0) {#woops, too far from each other
		$mapStatus = $INVALID_DIST;
	}
	elsif ($cnt>2){ #multiple pairs
		$mapStatus = $MULTIPLE_PAIRS;
		set_sam_bit_flag($mapStatus,-1,$cnt,\@samRep,$peMapLoc,0,$mapLib,-1);
	}
	return ($mapStatus,@samRep);
}
