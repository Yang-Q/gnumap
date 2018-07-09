#!/usr/bin/perl -w
#Note that read_id in an input sam file should not contain any pair_id tag but it should be retrieved from SAM FLAG bit!!!
use strict;
use warnings;
use Cwd;
use File::Copy;

our ($gMDUST_MAX_REP,$GNUMAP,$BLAT,$gDISC2OP, $SINGLE,$SEP_DISC2,$ALLOW_DISC2,$NO_DISC2,$gUpdXPF,%gHreadGroup) = ();

($GNUMAP,$BLAT) = (1,2);
($SINGLE,$SEP_DISC2,$ALLOW_DISC2,$NO_DISC2) = (0,1,2,3);

my $tmp = $ARGV[0]; #mapper_id
my $mapper_id = $GNUMAP;
if ($tmp eq "blat") {$mapper_id = $BLAT;}

$tmp = $ARGV[1]; chomp($tmp); #paired_end?
my $paired_end = 0;
if ($tmp eq "pe") {$paired_end = 1;}
my $gnuNtConvF = $ARGV[2]; 
$gnuNtConvF = $gnuNtConvF + 0;

my $inAmbSam = $ARGV[3]; #gnumap_sam_to_resolve
my $outBestSam = $ARGV[4]; #out_sam
my $scratchDir2sort = $ARGV[5]; #scratch dir
my $maxNumBest = $ARGV[6];
$gDISC2OP = $ARGV[7] + 0;


my $discordantOutD = $ARGV[8];
my @sepDiscPtrW=();
if ($gDISC2OP == $SEP_DISC2) {
	unless (-e $discordantOutD) {
		die "check if $discordantOutD exist!";
	}
	
	for (my $i=1;$i<3;$i++) {
		print STDERR "$discordantOutD/$i.sam\n";
		open($sepDiscPtrW[$i-1],">$discordantOutD/$i.sam")
	}
}
else {
	@sepDiscPtrW=(-1,-1);
}

$gUpdXPF = $ARGV[9] + 0;
my $multiMapF = $ARGV[10]+0;
my $filtRheadFn = $ARGV[11];

my 	$gFiltReadF = 1;
if ($filtRheadFn eq "X") {
	$gFiltReadF = 0;
}

$gMDUST_MAX_REP = $maxNumBest + 0;
my $x;

my $inAmbSamSz;
if (-e $inAmbSam) {
	$inAmbSamSz = -s $inAmbSam;
	if ($inAmbSamSz < 1){exit(0);}
}
else{exit(0);}

my $rndTag = int(rand(10000));
my $tmp_file = "${inAmbSam}-temp1-${rndTag}";

my $rndTag2 = int(rand(10000000));
my $scratchDir2sort2 = "$scratchDir2sort/$rndTag2";
unless (-d $scratchDir2sort2) {system("mkdir -p $scratchDir2sort2");}
print STDOUT "selecting the best entries among ambiguous gnumap sam ($inAmbSam)...\n"; #debug

#print STDOUT "sorting multiple mapped reads in sam file....\n"; #debug
my $sortDoneMsgFn = "${inAmbSam}-sortDone";
if (-e $sortDoneMsgFn) {
	system("mv $inAmbSam $tmp_file");
}
else {
	$x = system("cat $inAmbSam | sort -T $scratchDir2sort2 -n -k1 -o $tmp_file");
	if ( $? == -1 ){
		die "command failed: $!\n";
	}
	#print STDOUT "done!\n"; #debug
}

print STDERR "mapperID=$mapper_id\n";

my ($SAM_QUERY_STRAND,$FIRST_PAIR,$SECOND_PAIR,$FIRST_PAIR2,$SECOND_PAIR2) = (16,64,128,0,1);
my (@recSamLine,@recAlScore,@recFragLen,@recMapLoc)=();
my ($read_id,$flag,$ref_id,$strand,$ntConv,$map_loc,$mapq,$mate_map_loc,$frag_len,$ascore,$read_id_group,$entryIdx,$samRegExp,$pair_id,$parseError,$filtReadF,$filtReadFp);

open(SAMFP,"<$tmp_file") || die "cannot read $tmp_file after sorting $inAmbSam!\n";
open(OPTSAM,">$outBestSam") || die "cannot create $outBestSam!\n";

if ($gFiltReadF==1) {
	open($filtReadFp,">$filtRheadFn") || die "cannot create $filtRheadFn!\n";
}

my $samLine = <SAMFP>;
$samLine =~ /^(\S+)\t.*/;
$read_id_group = $1;
seek(SAMFP,0,0);
$entryIdx = 0;

$ntConv = "nn"; #it can be ct, ga, ag, or tc

while ($samLine = <SAMFP>) {
	$parseError = 1;
	if ($mapper_id == $GNUMAP) {
		if (($gnuNtConvF==1) && ($samLine =~ /^(\S+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t\S+\t(\d+)\t(-{0,1}\d+)\t\S+\t\S+\tAS:i:(\d+)\t.*\tXS:Z:(\S+)\t/)) {
			$ascore = $8 +0.0;
			($read_id,$flag,$ref_id,$map_loc,$mapq,$mate_map_loc,$frag_len,$ntConv)=($1,$2,$3,$4,$5,$6,$7,$9);
			$parseError = 0;
		}
		elsif ($samLine =~ /^(\S+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t\S+\t(\d+)\t(-{0,1}\d+)\t\S+\t\S+\tAS:i:(\d+)\t.*/) {
			$ascore = $8 +0.0;
			($read_id,$flag,$ref_id,$map_loc,$mapq,$mate_map_loc,$frag_len)=($1,$2,$3,$4,$5,$6,$7);
			$parseError = 0;
		}
	}
	elsif ($mapper_id == $BLAT) {
		if ($samLine =~ /^(\S+)\t(\d+)\t(\S+)\t(\d+)\t(\d+)\t\S+\t\S+\t(\d+)\t(-{0,1}\d+)\t\S+\t\S+\tAS:i:(\d+)\tPG:Z:blat/) {
			($read_id,$flag,$ref_id,$map_loc,$mapq,$mate_map_loc,$frag_len,$ascore)=($1,$2,$3,$4,$5,$6,$7,$8);
			$parseError = 0;
			
		}
	}
	else {
		die "check this sam file to see if it is corrupted?";
	}
	
	if ($parseError == 1) {die "check $samLine !";}
	
	$tmp = $flag & $SECOND_PAIR;
	$pair_id = $FIRST_PAIR2;
	if ($tmp == $SECOND_PAIR) {$pair_id = $SECOND_PAIR2;}

	$tmp = $flag & $SAM_QUERY_STRAND;
	$strand = 0;
	if ($tmp == $SAM_QUERY_STRAND) {$strand = 1;}
		
	if ($read_id ne $read_id_group) { #come to new entry and let's find the best optimals!
		$filtReadF = find_max_entry($mapper_id,\@recSamLine,\@recAlScore,\@recMapLoc,$entryIdx,$paired_end,$multiMapF,\*OPTSAM,\@sepDiscPtrW);
		if (($gFiltReadF==1) && $filtReadF) {
			printf { $filtReadFp } "%s\n", $read_id_group;
		}
		
		%gHreadGroup = (); 
		$entryIdx = 0;
		
		@recSamLine = ();
		@recAlScore = ();
		@recFragLen = ();
		@recMapLoc = ();
		
		$entryIdx = fill_read_group_buffer($paired_end,$read_id,$pair_id,$ref_id,$strand,$map_loc,$mate_map_loc,$samLine,$ascore,$frag_len,$entryIdx,$ntConv,\@recSamLine,\@recAlScore,\@recFragLen,\@recMapLoc);
		$read_id_group = $read_id;
	}
	else { #keep store them into a collection
		$entryIdx = fill_read_group_buffer($paired_end,$read_id,$pair_id,$ref_id,$strand,$map_loc,$mate_map_loc,$samLine,$ascore,$frag_len,$entryIdx,$ntConv,\@recSamLine,\@recAlScore,\@recFragLen,\@recMapLoc);
	}
}

$filtReadF = find_max_entry($mapper_id,\@recSamLine,\@recAlScore,\@recMapLoc,$entryIdx,$paired_end,$multiMapF,\*OPTSAM,\@sepDiscPtrW);
if (($gFiltReadF==1) && $filtReadF) {
	printf { $filtReadFp } "%s\n", $read_id_group;
}
		
close(SAMFP);
close(OPTSAM);
if ($gFiltReadF==1) {close($filtReadFp);}

system("mv $tmp_file $inAmbSam");
system("echo 'sorted sam files' > ${inAmbSam}-sortDone");
if (-d $scratchDir2sort2) {system("rm -rf $scratchDir2sort2");}

print STDOUT "done!\n";

#==========================================
sub find_max_entry {
	my ($mapper_id,$recSamLine, $recAlScore, $recMapLoc, $numList,$paired_end, $multiMapF,$filePtrW,$sepDiscPtrW) = @_;

	my ($i,$maxScore,$max_cnt,$idx,$samLine2,$filtReadF,$discF,$noScorePairi); $maxScore = 0.0;
	my @scoreArray = ();
	my @EmptySlot = ();
	my %hMapLoc = ();
	
	for ($i=0; $i<$numList; $i++) { #sum alignment score
		if ($paired_end) {
			unless ($recAlScore->[$i][0]) {$recAlScore->[$i][0]=0.0;}
			unless ($recAlScore->[$i][1]) {$recAlScore->[$i][1]=0.0;}
			$scoreArray[$i] = $recAlScore->[$i][0] + $recAlScore->[$i][1];
		}
		else {
			$scoreArray[$i] = $recAlScore->[$i][0];
		}

		if ($scoreArray[$i] > $maxScore) {
			$maxScore = $scoreArray[$i];
		}
	}
	
	#count total num of entries having maxScore
	$max_cnt = 0;
	for ($i=0; $i<$numList; $i++) {
		if ($scoreArray[$i] == $maxScore) {
			$noScorePairi = 0;
			if (($gDISC2OP != $SINGLE) && (($recAlScore->[$i][0]==0.0) || ($recAlScore->[$i][1]==0.0))) {
				$noScorePairi = 1;
			}
			$EmptySlot[$i]=$noScorePairi;
			$idx = $hMapLoc{$recMapLoc->[$i]};
			if (defined($idx)) {next;}
			$max_cnt++;
			$hMapLoc{$recMapLoc->[$i]} = 1;
		}
	}

	%hMapLoc = ();
	#my $maxCutoff=$maxScore*1.0;
	
	if ($max_cnt <= $gMDUST_MAX_REP) {
		$filtReadF = 0;
		for ($i=0; $i<$numList; $i++) {
			if ($scoreArray[$i] >= $maxScore) {
			
				#use hash to avoid mapped reads on the same location w/ the maxScore
				$idx = $hMapLoc{$recMapLoc->[$i]};

				if (defined($idx)) {next;}

				if (($EmptySlot[$i] == 0) || ($gDISC2OP != $NO_DISC2)) {
					for (my $j=0;$j<2;$j++) { #for each pair
						if ($recAlScore[$i][$j] > 0.0){
							$samLine2 = $recSamLine->[$i][$j];
							if ($gUpdXPF==1) {
								$samLine2 = update_XP_field_in_sam($recSamLine->[$i][$j],$mapper_id,$multiMapF,$max_cnt);
							}
							
							if (($EmptySlot[$i] == 0) || ($gDISC2OP == $ALLOW_DISC2)) {
								printf { $filePtrW } "%s", $samLine2;
							}
							elsif ($gDISC2OP == $SEP_DISC2) {
								$samLine2 = revert_sam_flag($samLine2,($j+1));
								printf { $sepDiscPtrW->[$j] } "%s", $samLine2; #qual
							}
						}
					}
				}
				
				$hMapLoc{$recMapLoc->[$i]} = 1;
			}
		}
	}
	else {
		$filtReadF = 1;
	}
	return $filtReadF;
}

#===================
sub revert_sam_flag {
	my ($samLine,$pair_id) = @_;
	my $SAM_QUERY_STRAND = 16;
	my ($tmp,$flag);
	
	my @entr=split('\t',$samLine);
	$flag = $entr[1]+0;
	$tmp = $flag & $SAM_QUERY_STRAND;
	$entr[1] = "0";
	if ($tmp == $SAM_QUERY_STRAND) {$entr[1] = "16";}
	
	$entr[0] = "$entr[0]/$pair_id";
	$samLine = join("\t",@entr);

	return $samLine;
}
#==========================================
sub update_XP_field_in_sam {
	my ($samLine,$mapper_id,$multiMapF,$dupCnt) = @_;
	my $pProb = 1.0/$dupCnt;
	my $probField = "XP:f:$pProb";
	
	if ($mapper_id == $GNUMAP) {
		#replace XP:f:x
		$samLine =~ s/XP:f:\S+/$probField/g;
	}
	elsif ($mapper_id == $BLAT) {
		#add XP:f:x
		chomp($samLine);
		
		if ($multiMapF && $dupCnt>1) {
			if ($samLine =~ m/^\S+\s+(\d+)\s+/) {
				$flag = $1+256;
				$samLine =~ s/\b$1\b/$flag/;
			}
		}

		$samLine = "$samLine\t$probField\n";
	}
	
	return ($samLine);
}

#==========================================
sub fill_read_group_buffer {

	my ($paired_end,$read_id,$pair_id,$ref_id,$strand,$map_loc,$mate_map_loc,$samLine,$ascore,$frag_len,$entryIdx,$ntConv,$recSamLine,$recAlScore,$recFragLen,$recMapLoc) = @_;
	
	my ($SAM_QUERY_STRAND,$FIRST_PAIR,$SECOND_PAIR,$FIRST_PAIR2,$SECOND_PAIR2) = (16,64,128,0,1);
	my ($hash_key,$entry);
	
	if ($paired_end==1) { #if a user wants to use the paired_end information, use a first pair as hash_key
		$hash_key = "$read_id<$ref_id<$map_loc<$mate_map_loc<$ntConv";
		if ($pair_id == $SECOND_PAIR2) {
			$hash_key = "$read_id<$ref_id<$mate_map_loc<$map_loc<$ntConv";
		}
	}
	else {
		$hash_key = "$read_id<$ref_id<$strand<$map_loc<$ntConv";
		$pair_id = 0;
	}
	
	$entry = $gHreadGroup{$hash_key};
	unless (defined($entry)) { #first time to see this entry (so, need initialize)
		$gHreadGroup{$hash_key} = $entryIdx;
		$entry = $entryIdx;
		$entryIdx += 1;
		$recAlScore->[$entry][0] = 0;
		$recAlScore->[$entry][1] = 0;
	}

	$recSamLine->[$entry][$pair_id] = $samLine;
	$recFragLen->[$entry][$pair_id] = $frag_len;
	$recAlScore->[$entry][$pair_id] = $ascore;
	
	if ($paired_end) {
		if ($pair_id == $FIRST_PAIR2) {
			$recMapLoc->[$entry] = "$ref_id:$map_loc:$mate_map_loc";
		}
		else {
			$recMapLoc->[$entry] = "$ref_id:$mate_map_loc:$map_loc";
		}
	}
	else {
		$recMapLoc->[$entry] = "$ref_id:$map_loc";
	}

	return $entryIdx;
}
