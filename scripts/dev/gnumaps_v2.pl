#!/usr/bin/perl -w
## TODO: 
#1) to perform a comprehensive test on the quality control module
#3) to combine two sams of blat before forming PER [to debug]
#4) reuse gnumap index

use strict;
use warnings FATAL => 'all';
use Getopt::Std;
use Getopt::Long;
use File::Basename;
use Cwd qw[abs_path];
use File::Spec;
use Time::HiRes qw( time );

######################################################################
# define global variables
######################################################################
our (@gRead,@gRef,@gSrchDim,@gBIN,@gOut,@gSens2map,$gREF_DELIM,$gSTRAND_DELIM,$gGNREF,$gGNUIDX,$gTNREF,$gGNREF_BASE,$gGNREF_FEXT,$gTNREF_BASE,$gTNREF_FEXT,$rFILE_BASE,$rDIR,$rFEXT,$rDIR_QC,$rREAD_TYPE,$rLIB_TYPE,$rPER,$rPER_AVE_DIST,$rMAX_FRAG_LEN,$rCONV,$rPHRED,$oWORKD,$oALIGND,$oGNU_IDX_FEXT,$oREFTYPE,$oTEMPD,$oPI,$oMAPPED,$oSTRAND,$oPCR,$oPER,$oMEXT,$oZIP,$oBAM,$ALL,$UP,$DN,$UNSTRAND,$STRAND1,$STRAND2,$bQC,$bGNU,$bPER,$bBEST,$bGNUSGR,$bBLAT,$bSTP,$bSAMT,$bPSL2SAM,$bUNMAP2FA,$b2BIT,$bqADAP,$bqADAP2,$bqCUTADPT,$bqPRTSEQ,$bgBIN,$bgOPRFX,$bgGNU_ACC,$bgGNU_KMER,$bgGNU_MAX_HASH,$bgBIN_LIB,$bgGNU_SLIDE,$bgGNU_TOPK_HASH,$bgGNU_SW,$bgiIDX_FEXT,$bgiOPRFX,$bgiNULLREAD,$bgiREVCOMP_SAM,$bpBIN,$bpPRFX,$b1BIN,$b1PRFX,$b1NUM_AMB,$bgsBIN,$bgsOPRFX,$bgsPOSTPROC,$bbBIN,$bbOPRFX,$bbPERIDN,$bbTILESZ,$bbSTEPSZ,$bbREPMAT,$bbMAX_INTRON,$bbREAD_FEXT,$bstBIN,$bstOPRFX,$bstBIN_LIB,$bpsBIN,$bpsOPRFX,$buBIN,$buOPRFX,$b2BIN,$b2OPRFX,$b2IDX_FEXT,$gTHREADS,$DNA,$RNA,$gDummy,$NAs,$VALID,$EMPTY,$NOT_EXIST,$gDebug);

($gREF_DELIM,$gSTRAND_DELIM,$gGNREF,$gGNUIDX,$gTNREF,$gGNREF_BASE,$gGNREF_FEXT,$gTNREF_BASE,$gTNREF_FEXT)=(0,1,2,3,4,5,6,7,8,9);
($rFILE_BASE,$rDIR,$rFEXT,$rDIR_QC,$rREAD_TYPE,$rLIB_TYPE,$rPER,$rPER_AVE_DIST,$rMAX_FRAG_LEN,$rCONV,$rPHRED) = (0,1,2,3,4,5,6,7,8,9,10);
($oWORKD,$oALIGND,$oGNU_IDX_FEXT,$oREFTYPE,$oTEMPD,$oPI,$oMAPPED,$oSTRAND,$oPCR,$oPER,$oMEXT,$oZIP,$oBAM)=(0,1,2,3,4,5,6,7,8,9,10,11,12); #gOut
($ALL,$UP,$DN) = (0,1,2); #gSens2map
($UNSTRAND,$STRAND1,$STRAND2) = (0,1,2);
($bQC,$bGNU,$bPER,$bBEST,$bGNUSGR,$bBLAT,$bSTP,$bSAMT,$bPSL2SAM,$bUNMAP2FA,$b2BIT) = (0,1,2,3,4,5,6,7,8,9,10);
($bqADAP, $bqADAP2, $bqCUTADPT, $bqPRTSEQ) = (0,1,2,3); #bQC subidx
($bgBIN,$bgOPRFX,$bgGNU_ACC,$bgGNU_KMER,$bgGNU_MAX_HASH,$bgBIN_LIB,$bgGNU_SLIDE,$bgGNU_TOPK_HASH,$bgGNU_SW,$bgiIDX_FEXT,$bgiOPRFX,$bgiNULLREAD,$bgiREVCOMP_SAM) = (0,1,2,3,4,5,6,7,8,9,10,11,12); #bGNU subidx
($bpBIN,$bpPRFX) = (0,1); #bPER
($b1BIN,$b1PRFX,$b1NUM_AMB) = (0,1,2); #bBEST
($bgsBIN,$bgsOPRFX,$bgsPOSTPROC) = (0,1,2); #bGNUSGR
($bbBIN,$bbOPRFX,$bbPERIDN,$bbTILESZ,$bbSTEPSZ,$bbREPMAT,$bbMAX_INTRON,$bbREAD_FEXT) = (0,1,2,3,4,5,6,7); #bBLAT
($bstBIN,$bstOPRFX,$bstBIN_LIB) = (0,1,2); #bSTP or bSAMT
($bpsBIN,$bpsOPRFX) = (0,1); #bPSL2SAM
($buBIN,$buOPRFX) = (0,1); #bUNMAP2FA
($b2BIN,$b2OPRFX,$b2IDX_FEXT) = (0,1,2); #b2BIT
$gTHREADS=0; $gDummy='X'; $NAs = 'X';
($DNA,$RNA)=("dna","rna");
($VALID,$EMPTY,$NOT_EXIST) = (0,1,2);
######################################################################
# parse arguments
######################################################################

my @usage;
push @usage, "Usage: " . basename($0) . " [options]\n";
push @usage, " --help     			Displays this information\n";
push @usage, " --genome			FASTA format reference genome file[required]\n";
push @usage, " --transc			a directory that contains FASTA format reference transcriptome files (Refer to README)\n";
push @usage, " --pair_1			if pair-end read specify _1.txt\n";
push @usage, " --pair_2			if pair-end read specify _2.txt\n";
push @usage, " --single			if single-end read specify the read file here[read1,2 or single required]\n";
push @usage, " --qcontrol			quality control enabled and make sure that all necessary programs are ready[0]\n";
push @usage, " --adaptor			adaptor sequence\n";
push @usage, " --adaptor2			adaptor sequence if the one for 2nd is different\n";
push @usage, " --lib_type			it can be [unstrand], cr1 or wt1\n";
push @usage, " --read_type			it can be [dna] or rna\n";
push @usage, " --phred_offset			speicify read quality offset. note that set to 0 if reads are in fasta file format (0, [33], or 64)\n";
push @usage, " --per_dist			if paired read, specify a distance between two reads[300]\n";
push @usage, " --num_amb			speicify the number of mapping locations to be reported [500]\n";
push @usage, " --map_quality			specify the mapping option, e.g., sensitive,balanced,fast\n";
push @usage, " --nt_conv			if the read is bisulfited or you want to explore rna edits, specify bs(bisulfited reads) or a2i (to discover RNA edit events (A-to-I))\n";
push @usage, " --skip_blat			if you want to skip blat mapping, set to 1 [0]\n";
push @usage, " --pileup			if you want to convert sam to bam file, set to 1(normal) or 2(flip 2nd read to the watson strand for a better representation in IGB) [0].\n";
push @usage, " --mpi				set to 1 if you want to enable mpi[0]\n";
push @usage, " --num_threads			the number of cores to enable[1]\n";
push @usage, " --outdir			an output directory[required]\n";
push @usage, " --debug			set to 1 if you want to perform in a debug mode\n";


my ( $help, $genome, $transc, $pair_1, $pair_2, $single, $qcontrol, $adap, $adap2, $lib_type, $read_type, $phred_offset, $per_dist, $nt_conv, $pileup1, $map_quality, $num_amb, $skip_blat, $mpi, $num_threads, $outdir, $debug_mode);

GetOptions(
	'help' => \$help,
	'genome=s' => \$genome,
	'transc=s' => \$transc,
	'pair_1=s' => \$pair_1,
	'pair_2=s' => \$pair_2,
	'single=s' => \$single,
	'qcontrol=i' => \$qcontrol,
	'adaptor=s'  => \$adap,
	'adaptor2=s' => \$adap2,
	'lib_type=s' => \$lib_type,
	'read_type=s' => \$read_type,
	'phred_offset=i' => \$phred_offset,
	'per_dist=i' => \$per_dist,
	'nt_conv=s' => \$nt_conv,
	'map_quality=s' => \$map_quality,
	'num_amb=i' => \$num_amb,
	'skip_blat=i' => \$skip_blat,
	'mpi=i' => \$mpi,
	'num_threads=i' => \$num_threads,
	'outdir=s' => \$outdir,
	'pileup=i' => \$pileup1,
	'debug=i' => \$debug_mode
);
not defined $help or die @usage;
defined $genome or die @usage;
defined $outdir or die @usage;

if (not defined $qcontrol){$qcontrol=0;}
if (not defined $lib_type){$lib_type="unstrand";}
if (not defined $read_type){$read_type="dna";}
if (not defined $phred_offset){$phred_offset=33;}
if (not defined $nt_conv){$nt_conv=$NAs;}
if (not defined $skip_blat){$skip_blat=0;}
if (not defined $pileup1){$pileup1=0;}
if (not defined $mpi){$mpi=0;}
if (not defined $num_threads){$num_threads=1;}
if (not defined $map_quality){$map_quality="balanced";}
if (not defined $num_amb){$num_amb=500;}
if (not defined $debug_mode){$debug_mode=0;}

########################################
# check input condition and requirements
########################################
my ( $read1, $read2, $installD ) = parse_arguments($genome,$transc,$pair_1,$pair_2,$single,$qcontrol,$adap,$adap2,$lib_type,$read_type,$phred_offset, $per_dist, $nt_conv, $map_quality, $num_amb, $mpi, $num_threads,$outdir,$debug_mode);

my $start = gnumap_pipeline_logo(0.0,0,"",$installD);

my @Reads = ($read1,$read2);

################################################################################################
# [optional but strongly recommened for BLAT to run on Roche454] trim adaptors and low quality.
# (in order to enable this feature, fastqc,cutadpt,prinseq should be installed!)
################################################################################################
if (  $qcontrol == 1 ) {
	gnumap_qcontrol( \@Reads ); #Reads array will be updated to point out reads under quality control
}

########################################
# check if ref_genome has index
########################################
my $gnuRefIdxD = build_ref_index_gnumap( $gRef[$gGNREF] ); 
my $blatGenRefIdx  = "X";
my $blatTrnRefIdx = "X";

#build blat index on ref genome sequence!
$blatGenRefIdx = build_ref_index_blat( $gRef[$gGNREF], 0 );
if ( $gRef[$gREF_DELIM] == 2 ) {
	$blatTrnRefIdx = build_ref_index_blat( $gRef[$gTNREF], 1 ) ;
}


##########################################
# gnumapping reads to genome
##########################################
my $gnuOutPrefix = run_gnumap( $gRef[$gGNREF], \@Reads );

##########################################
# postprocessing on gnumap mapping results
##########################################
my @Chrs = cmd_grep_awk_get_list_with_pat( $gRef[$gGNREF], '>' );
my $gnuAln = postprocessing( \@Chrs, $gBIN[$bGNU][$bgOPRFX], $gnuOutPrefix, $gOut[$oSTRAND][0], $gOut[$oSTRAND][1], $gOut[$oPI][$gRead[$rPER]+1] );

##########################################
# if it is rna, obtain leftover reads,
#   otherwise, go to sam2bam to pileup process
##########################################

my ($isEmpty, $unmappedD);
my $blatSamMappedFn = "";

if ($skip_blat==0) {
	($isEmpty,$unmappedD) = get_unampped_reads( $gnuOutPrefix, $gnuAln );

	##########################################
	# BLAT aligns unmapped reads to both genome
	#  and transcriptome
	##########################################
	if ($isEmpty==0) {
		$blatSamMappedFn =run_pblat_gen_tran3( $blatGenRefIdx, $blatTrnRefIdx, $unmappedD,$gRef[$gGNREF] );
	}
}

##########################################
# merge gnumap and blat sam file,
# then, it finally converts to bam file
##########################################
my ($outD,$samFn) = cmd_merge_gnuBlat2( $gRef[$gGNREF], $gnuAln, $blatSamMappedFn, $pileup1);

################################################
# clean up all temp or intermediately-gen files
################################################
if ($gDebug==0){
	cleanup_files($gOut[$oALIGND],$gBIN[$bGNU][$bgOPRFX],$gBIN[$bBLAT][$bbOPRFX]); #!debug
}

my $logMsg = "check output files in $outD\n$samFn\n";
$gDummy = gnumap_pipeline_logo($start,1,$logMsg,$installD);

	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
	
#===========================
sub cmd_merge_gnuBlat2 {

	my ( $refDfn, $sam1, $sam2, $pileup1) = @_;
	my ( $cmdLine, $sam2D, $samD, $bamFn, $s, $ntConvSuffix, $samStFn, $mergedSamDfn, $emptyFile, $prevD );

	
	my $pIdx = $gRead[$rPER]+1;
	$sam2D = "$gOut[$oALIGND]/$gBIN[$bSAMT][$bstOPRFX]";
	unless ( -d $sam2D ) { system("mkdir -p $sam2D"); }
	$mergedSamDfn = "$sam2D/$gOut[$oPI][$pIdx].$gOut[$oMEXT][0]";
		
	unless ( -e $sam1 ) {
		die "check file $sam1!\n";
	}
	
	if ($sam2 eq "") {
		rename $sam1, $mergedSamDfn;
	}
	else { #merge two files
		$cmdLine = "\ncat $sam1 $sam2 > $mergedSamDfn\n";
		if ($gDebug==0){
			$cmdLine = "$cmdLine\nrm -rf $sam1 $sam2\n"; #!debug
		}
		system($cmdLine);
	}
	
	#for nt_conv flag, we apply gnumap sam2sgr and get some conversion stats across genome
	my $gnuSamStrand = "XS:Z";
	if ($pileup1>0){
		if ($gRead[$rCONV] ne $NAs) {
			for ($s=0;$s<$gRef[$gSTRAND_DELIM];$s++){
				($ntConvSuffix,$gDummy) = select_gnumap_option($gRead[$rCONV], $s, 0);
				$samStFn = "$sam2D/$gOut[$oPI][$pIdx]-$ntConvSuffix.$gOut[$oMEXT][0]"; #to prepare sam mapped to nt_conv idx
				cmd_grep($mergedSamDfn,"$gnuSamStrand:$ntConvSuffix",$samStFn); #grep only selected sam
				($prevD,$bamFn) = cmd_sam2bam_full($refDfn,$samStFn,$pileup1,"keep");
				($prevD,$gDummy) = cmd_sam2gmp($refDfn,$ntConvSuffix,$samStFn,"zip");
			}
		}
		else {
			($prevD,$bamFn) = cmd_sam2bam_full($refDfn,$mergedSamDfn,$pileup1,"keep");
		}
	}
	return ($sam2D,$mergedSamDfn);
}

#=====================
sub gnumap_qcontrol {
	my ($Reads) = @_;

	my @at_tmp = ();
	my ( $at_good_tmp,$at_bad_tmp, $p, $pi, $cmd );

	#trim adaptors if needs it
	for ( $p = 0 ; $p < $gRead[$rPER] ; $p++ ) {
		$pi = $p + 1;
		$at_tmp[$p] = "$gRead[$rDIR_QC]/${pi}-at.fq";
		if (defined $gBIN[$bQC][$bqADAP+$p]) {
			$cmd = "$gBIN[$bQC][$bqCUTADPT] -f fastq -O 2 -a $gBIN[$bQC][$bqCUTADPT] $Reads->[$p] -o $at_tmp[$p]";
		}
		else {
			$cmd = "ln -s $Reads->[$p] $at_tmp[$p]";
		}
		system($cmd);
	}

	#trim low quality bases or homopolymer
	$cmd         = "";
	$at_good_tmp = "$gRead[$rDIR_QC]/at-good";
	$at_bad_tmp  = "$gRead[$rDIR_QC]/at-bad";

	if ( $gRead[$rPER] == 2 ) {    #output to ./qc directory
		$cmd = "$gBIN[$bQC][$bqPRTSEQ] -fastq $at_tmp[0] -fastq2 $at_tmp[1] -out_fomrat 3 -min_len 1 -no_qual_header -lc_method entropy -lc_threshold 50 -trim_qual_type min -trim_qual_rule lt -trim_qual_left 3 -trim_qual_right 3 -out_good $at_good_tmp -out_bad $at_bad_tmp";

		system($cmd);

		#change file name
		$Reads->[0] = "$gRead[$rDIR_QC]/1.fq";
		rename "$gRead[$rDIR_QC]/at-good_1_singletons.fastq", $Reads->[0];

		$Reads->[1] = "$gRead[$rDIR_QC]/2.fq";
		rename "$gRead[$rDIR_QC]/at-good_2_singletons.fastq", $Reads->[1];

	}
	else {
		$cmd = "$gBIN[$bQC][$bqPRTSEQ] -fastq $at_tmp[0] -out_fomrat 3 -min_len 1 -no_qual_header -lc_method entropy -lc_threshold 50 -trim_qual_type min -trim_qual_rule lt -trim_qual_left 3 -trim_qual_right 3 -out_good $at_good_tmp -out_bad $at_bad_tmp";

		system($cmd);

		#change file name
		$Reads->[0] = "$gRead[$rDIR_QC]/1.fq";
		rename "$gRead[$rDIR_QC]/at-good.fastq", $Reads->[0];
	}
	
	for ( $p = 0 ; $p < $gRead[$rPER] ; $p++ ) { unlink( $at_tmp[$p] ); }

}

#=========================================================
# merge two sam files mapped to both genome and transcriptome and then form a PER
sub run_pblat_gen_tran3 {
	my ( $blatGenRefIdx, $blatTrnRefIdx, $unmappedD, $ref) = @_;

	  my (
		$blatTGenAln, $tmp,	$blatMappedD, $blatSamFn, $stat,
		$blatSamMappedFn, $alignBestD,   $cmd, $p,
		$blatGenSamD, $samExt, $blatTranSamD, $subD2opt
	  );
	my (@Chrs,@blatGenAln,@blatTranAln);
	
	#mapping to reference genome ------------------------------
	($blatGenSamD,$subD2opt) = pBLAT_wrapper( $gOut[$oREFTYPE][0], $blatGenRefIdx, $unmappedD, 1 );
	@Chrs = cmd_grep_awk_get_list_with_pat( $ref, '>' );
	my @pId    = flatten( $gOut[$oPI] );
	$samExt = ".$gOut[$oMEXT][0]";      #only take sam file format

	#mapping to transcriptome sequence-------------------------
	if ( $gRef[$gREF_DELIM] == 2 ) { 
		($blatTranSamD,$subD2opt) = pBLAT_wrapper( $gOut[$oREFTYPE][1], $blatTrnRefIdx, $unmappedD, 1 );

		$alignBestD = "$blatTranSamD/$gOut[$oMAPPED][1]";    #mappped
		unless ( -d $alignBestD ) { system("mkdir -p $alignBestD"); }

		for ($p=0;$p<$gRead[$rPER];$p++){
			$blatGenAln[$p] = "$blatGenSamD/$gOut[$oSTRAND][2]/$pId[$p]${samExt}";
			$blatTranAln[$p] = "$blatTranSamD/$gOut[$oSTRAND][2]/$pId[$p]${samExt}";
			($stat,$tmp) = conv_transc2genome_coordinate($blatTranAln[$p],$p+1);
			if ($stat>0){
				$cmd = "cat $tmp >> $blatGenAln[$p]\n"; #merge to blat.gen.sam
				if ($gDebug==0){
					$cmd = "$cmd\nrm -rf $tmp\n"; #!debug
				}
				system($cmd);
			}
		}
		if ($gDebug==0){
			system("rm -rf $blatTranSamD"); #!debug
		}
	}
	
	$blatTGenAln = postprocessing( \@Chrs, $gBIN[$bBLAT][$bbOPRFX], $blatGenSamD, $subD2opt, 'X', $gOut[$oPI][$gRead[$rPER]+1]);
	
	#to prepare some final output file name -------------------------
	$blatMappedD     = "$gOut[$oALIGND]/$gBIN[$bBLAT][$bbOPRFX]/$gOut[$oMAPPED][1]";
	unless (-d $blatMappedD) {system("mkdir -p $blatMappedD");}
	$blatSamFn       = "$gOut[$oPI][$gRead[$rPER]+1].$gOut[$oMEXT][0]";
	$blatSamMappedFn = "$blatMappedD/$blatSamFn";
	rename $blatTGenAln, $blatSamMappedFn;
	return $blatSamMappedFn;
}

sub build_ref_index_gnumap {
	my ($refDfn) = @_;

  $gRef[$gGNUIDX] = [ ( 'x', 'x' ) ];

	my ( $s, $cmd, $gnuNtConvOpt, $ntConvSuffix );

	#gnumap : ref_genome : nt_conv
	my $gnuRefIdxD = "$gOut[$oWORKD]/$gOut[$oREFTYPE][0]/$gBIN[$bGNU][$bgiOPRFX]";

	unless ( -d $gnuRefIdxD ) { system("mkdir -p $gnuRefIdxD"); }

	my $dummy_out = "$gOut[$oTEMPD]/nullReadOutTmp";

	#if there is any nt conversion, we prepare two refIdx
	for ( $s = 0 ; $s < $gRef[$gSTRAND_DELIM] ; $s++ ) { #if nt_conv is defined, gnumap should prepare two (nt-converted) indexes

		( $ntConvSuffix, $gnuNtConvOpt ) = select_gnumap_option( $gRead[$rCONV], $s, 1 );

		$gRef[$gGNUIDX][$s] = "$gnuRefIdxD/$gRef[$gGNREF_BASE]-$ntConvSuffix-m$gBIN[$bGNU][$bgGNU_KMER]-s$gBIN[$bGNU][$bgGNU_SLIDE]-h$gBIN[$bGNU][$bgGNU_MAX_HASH].$gBIN[$bGNU][$bgOPRFX]";

		unless ( -e $gRef[$gGNUIDX][$s] ) {
			$cmd = "export LD_LIBRARY_PATH=$gBIN[$bGNU][$bgBIN_LIB]:\$LD_LIBRARY_PATH\n";

			$cmd = "$cmd\n$gBIN[$bGNU][$bgBIN] -g $refDfn --save=$gRef[$gGNUIDX][$s] -o $dummy_out -m $gBIN[$bGNU][$bgGNU_KMER] -s $gBIN[$bGNU][$bgGNU_SLIDE] -c 1 -h $gBIN[$bGNU][$bgGNU_MAX_HASH] -t 1 -a 0.95 $gnuNtConvOpt -y 2 $gBIN[$bGNU][$bgiNULLREAD]\n";

			#print "$cmd\n";    #!debug
			system($cmd);
		}
	}
	unlink("${dummy_out}.sam");
	return $gnuRefIdxD;
}

#====================

sub build_ref_index_blat {
	my ( $refDfn, $refType_i ) = @_;    #$refDfn can be file or directory

	my ( $blatIdxDfn, $cmdLine, $chr, $out, $tRef);
	my (@Chrs);
	my $blatIdxD = "$gOut[$oWORKD]/$gOut[$oREFTYPE][$refType_i]/$gBIN[$b2BIT][$b2OPRFX]";

	unless ( -d $blatIdxD ) { system("mkdir -p $blatIdxD"); }

#if ( $refType_i == 0 ) {
	$blatIdxDfn = "$blatIdxD/$gRef[$gGNREF_BASE].$gBIN[$b2BIT][$b2IDX_FEXT]";
	unless ( -e $blatIdxDfn ) {
		$cmdLine = "$gBIN[$b2BIT][$b2BIN] $refDfn $blatIdxDfn\n";
		system($cmdLine);
	}
	$out = $blatIdxDfn;
#}
# 	else {
# 		@Chrs = listFilesInDir( $refDfn, 0, ".$gBIN[$bBLAT][$bbREAD_FEXT]", 1 );    #TODO support all other file extension
# 		foreach $chr (@Chrs) {
# 			$blatIdxDfn = "$blatIdxD/$chr.$gBIN[$b2BIT][$b2IDX_FEXT]";
# 			$tRef       = "$refDfn/$chr.$gBIN[$bBLAT][$bbREAD_FEXT]";
# 			$cmdLine    = "$gBIN[$b2BIT][$b2BIN] $tRef $blatIdxDfn\n";
# 			system($cmdLine);
# 		}
# 		$out = $blatIdxD;
# 	}
	return $out;
}

#====================
sub getFileName {
	my ( $fn, $fileExt ) = @_;
	my $fn0 = "";

	if ( $fn =~ /(\S+)\.$fileExt/ ) { $fn0 = $1; }
	else { die "check the input file extension, $fn \n"; }
	return $fn0;
}

#====================
sub cmd_grep {
	my ( $inputDfn, $pat, $outDfn ) = @_;
	my $cmdLine = "grep -iP \'$pat\' $inputDfn > $outDfn\n";
	system($cmdLine);
}

#====================

sub cmd_sam2bam_full {
	my ($refFn, $samFn, $flp2F, $handleInAfter) = @_;
	my $samFn2flp;
	
	#check if this is valid sam file (TODO)
	my $samFnBase = getFileName( $samFn, $gOut[$oMEXT][0] );
	unless ($samFnBase) { die "check samFn (${samFn}) format!\n"; }
	unless ( -e $refFn ) { die "the ref file, ${refFn}, does not exist!\n"; }
	
	#flip 2nd pair read for PER to display properly in IGB (asked by David Nix)
	if ($flp2F == 2) {
		$samFnBase = "$samFnBase-$gOut[$oPI][4]";
		$samFn2flp = "$samFnBase.$gOut[$oMEXT][0]";
		($gDummy, $gDummy) = cmd_samTransPar_java(0, $samFn, $samFn2flp, $flp2F);
		$samFn = $samFn2flp;
	}

	my $tmpFn = "$samFnBase-tmp";
	my $prevD = 0;
	my $progCmd = "$gBIN[$bSAMT][$bstBIN] view -bT ${refFn} ${samFn} > ${tmpFn}\n";
	$progCmd = "${progCmd}\n$gBIN[$bSAMT][$bstBIN] sort ${tmpFn} ${samFnBase}\n";
	my $bamFn = "$samFnBase.$gOut[$oBAM][0]\n";
	
	$progCmd = "${progCmd}\n$gBIN[$bSAMT][$bstBIN] index ${bamFn}\n";
	( $gDummy, $progCmd, $gDummy ) = cmd_zip_util( $progCmd, $tmpFn, "delete" );
	( $gDummy, $progCmd, $gDummy ) = cmd_zip_util( $progCmd, $samFn, $handleInAfter );
	
	system($progCmd);
	
	return ( $prevD, $bamFn );
}

#===========================================================================
sub pBLAT_wrapper {
	my ( $refType, $idxFn, $faD, $onlyBest ) = @_;

	my ( $fa, $j, $cmd2qsub, $faDfn_j, $tmp, $blatSamDfn_j, $pslx_j, $samD, $psl2samOpt,
		$faDfn, $DelimSTs, $DelimEnds, $rBase, $pslD, $p, $alignBestDfn, $blatSamD,$pi,$blatSamDfn,$cmd,$pslDfn,$aliBaseD,$aliSubD,$handleInAfter);
	my $samExt = $gOut[$oMEXT][0];

	$handleInAfter="keep";
	if ($gDebug==0){
		$handleInAfter="delete";
	}
	my $blatOption = "";
	
	if ($refType eq $gOut[$oREFTYPE][0]) { #genome
		$blatOption = "-q=$gRead[$rREAD_TYPE]";
		$psl2samOpt = $gRead[$rREAD_TYPE];
	}
	else { #transcriptome
		$blatOption = "-q=dna"; #no intron allowed for transcriptome ref 
		$psl2samOpt = "dna";
	}
	$aliBaseD = "$gOut[$oALIGND]/$gBIN[$bBLAT][$bbOPRFX]/$refType";
	$aliSubD=$gOut[$oSTRAND][2];
	$blatSamD = "$aliBaseD/$aliSubD"; #temp directory for blat psl file
	my $blatAccCutoff=$gBIN[$bGNU][$bgGNU_ACC]-0.2;

	unless ( -d $blatSamD ) { system("mkdir -p $blatSamD"); }

	for ( $p = 0 ; $p < $gRead[$rPER] ; $p++ ) {# for each paired read (for single end, only _1.txt will be processed)
		$pi = $p + 1;
		$blatSamDfn = "$blatSamD/$pi.${samExt}";
		unless (-e $blatSamDfn) {
			$pslDfn = cmd_blat_pslx( $idxFn, $faD, $pi, $blatOption );
			$faDfn = "$faD/$pi.$gBIN[$bBLAT][$bbREAD_FEXT]";
			$cmd = "$gBIN[$bPSL2SAM][$bpsBIN] $faDfn $pslDfn $blatSamDfn $psl2samOpt $gRead[$rCONV] $blatAccCutoff 0 $p $gRead[$rLIB_TYPE]\n";
			#print "$cmd\n" #debug
			system($cmd);
		}
		
		$cmd="";
		if ( $onlyBest == 1 ) {
			$alignBestDfn = "$blatSamD/$pi.${samExt}.tmp";
			cmd_resolve_ambig_sam2( $gBIN[$bBLAT][$bbOPRFX], $gOut[$oPI][2], $blatSamDfn, $alignBestDfn, $handleInAfter); #!debug
			if (-e $alignBestDfn){
				rename $alignBestDfn, $blatSamDfn;
			}
		}
		
		if (defined($pslDfn) && (-e $pslDfn)) {
			if ($gDebug==0){
				$cmd = "rm -rf $pslDfn\n";
				system($cmd); #!debug
			}
		}
	}
	return ($aliBaseD,$aliSubD);
}

#================================
sub conv_transc2genome_coordinate {

	#for transcriptome, bring it back to genome coordinate.
	my ($blatTranAln,$pairId) = @_;
	my $filesize = 0;
	
	if (-e $blatTranAln){
		$filesize = -s $blatTranAln;
		
		if ($filesize<1){
			return (0,$blatTranAln);
		}
	}
	else {
		print "$blatTranAln does not exist...\n";
		return (-1,$blatTranAln);
	}
	
	my $blatTranAln2 = "$blatTranAln.stpTmp.sam";
	my $prevD;
	( $prevD, $blatTranAln ) = cmd_samTransPar_java( $pairId, $blatTranAln, $blatTranAln2, 0 );
	return (1,$blatTranAln);
}

#==================
sub append_read_pid_to_sam {
	my ($sam,$pairId) = @_;
	my (@s)=();
	my ($sam2,$line);
	
	$sam2 = "$sam.pairId.tmp";
	open(SAM,"<$sam") || die "cannot open input sam file, $sam!\n";
	open(SAM2,">$sam2") || die "cannot open input sam file, $sam2!\n";
	
	while ($line = <SAM>){
		@s = split('\t',$line);
		$s[0]="$s[0]/$pairId";
		print SAM2 join("\t", @s);
	}
	close(SAM);
	close(SAM2);
	rename $sam2, $sam;
}
#===========================================================
sub cmd_samTransPar_java {
	my ($pairId, $tranSam, $outSamDfn, $flp2F ) = @_;

	my ( $prevD, $samT2Gen, $tGenSam, $progCmd );
	$progCmd = "export CLASSPATH=$gBIN[$bSTP][$bstBIN_LIB]:.:\$CLASSPATH";

	my $stpOpt = "";
	if ( $flp2F == 1 ) { $stpOpt = "-r"; }

	$progCmd = "$progCmd\njava -Xmx1500M -jar $gBIN[$bSTP][$bstBIN] -f $tranSam -a 400 -n 25 -d $stpOpt -s $outSamDfn";
	print "$progCmd\n"; #debug
	system($progCmd);
	
	$tGenSam = "$outSamDfn.gz";
	( $prevD, $progCmd, $samT2Gen ) = cmd_zip_util( "", $tGenSam, "unzip-n-delete" );

	#remove header in sam file
	$progCmd = "$progCmd\nsed -i '1,4d' $samT2Gen\n";    #todo: this should be updated!!!
	print "$progCmd\n"; #debug
	system($progCmd);
	
	#todo: if $pairId defined, append its original id_tag
	if ($pairId>0) {
		append_read_pid_to_sam($samT2Gen,$pairId);
	}
	
	$progCmd="";
	if ( $flp2F == 0 ) {
		( $prevD, $progCmd, $gDummy ) = cmd_zip_util( $progCmd, $tranSam, "delete" );
		$progCmd = "$progCmd\nmv $samT2Gen $tranSam\n";
	}
	else {
		( $prevD, $progCmd, $gDummy ) = cmd_zip_util( $progCmd, $tranSam, "zip" );
		$tranSam = $samT2Gen;
	}
	print "$progCmd\n"; #debug
	system($progCmd);
	return ( $prevD, $tranSam );
}

sub cmd_blat_pslx {
	my ( $refIdx, $unmappedD, $readBase, $blatOption ) = @_;

	my $cmdLine = "cd $unmappedD\n";
	my $outPsl  = "$readBase.psl";

	$cmdLine = "$cmdLine\n$gBIN[$bBLAT][$bbBIN] $refIdx -threads=$gTHREADS -minIdentity=$gBIN[$bBLAT][$bbPERIDN] -tileSize=$gBIN[$bBLAT][$bbTILESZ] -stepSize=$gBIN[$bBLAT][$bbSTEPSZ] -repMatch=$gBIN[$bBLAT][$bbREPMAT] $blatOption -out=pslx $unmappedD/$readBase.$gBIN[$bBLAT][$bbREAD_FEXT] $outPsl\n"; #running blat

	system($cmdLine);
	$outPsl = "$unmappedD/$outPsl";
	return ($outPsl);

}

#====================
#Chrs : list of all chromosomes
#mapper: which aligner do you use?
#samDprefix: alignment base directory
#subD2opt1: subdirectory under samDprefix to indicate where nt-converted strand (watson or crick or both) is used.
sub postprocessing {
	my ( $Chrs, $mapper, $samDprefix, $subD2opt1, $subD2opt2, $perTag) = @_;

	my ($p,$cmdLine,$alignBestD,$alignBestDfn,$alignMergedSamDfn,$samDfn,$pairIdx);

	my $samExt = ".$gOut[$oMEXT][0]";      #only take sam file format
	my @pId    = flatten( $gOut[$oPI] );

	$alignBestD = "$samDprefix/$gOut[$oMAPPED][1]";    #mappped
	unless ( -d $alignBestD ) { system("mkdir -p $alignBestD"); }
	
	$alignBestDfn = "$alignBestD/$perTag${samExt}";

	if ( $gRead[$rPER] == 2 ) {            #merge in per mode
		demuxPeSam( $Chrs, $samDprefix, $subD2opt1, $subD2opt2 );
		perSam_core( $Chrs, $samDprefix, $subD2opt1, $subD2opt2, $mapper );
	}
	
	#--------------------
	#simply merge
	$alignMergedSamDfn = muxSam( $samDprefix );
	
	
	#--------------------
	#pickup best ones
	my $handleInAfter="keep";
	if ($gDebug==0){
		$handleInAfter="delete";
	}
	cmd_resolve_ambig_sam2( $mapper, $perTag, $alignMergedSamDfn, $alignBestDfn, $handleInAfter ); #!debug

	#now, clean up the original gnumap sam files
	$samDfn  = "$samDprefix/*/?${samExt}";
	if ($gDebug==0){
		$cmdLine = "rm -rf $samDfn\n";
		system($cmdLine); #debug
	}
	return $alignBestDfn;
}

#================================================================
sub cmd_demux_strand {

	my ( $cmdLine, $inSamD, $inSamfn ) = @_;
	my @Strand = ( "watson", "crick" );
	my ( $k, $outSamD, $gnuSamStrand, $ntConv );
	my (@SamFlag);

	$gnuSamStrand = "XS:Z";

	for ( $k = 0 ; $k < 2 ; $k++ ) {

		$outSamD = "$inSamD/$Strand[$k]";
		unless ( -d $outSamD ) { system("mkdir -p $outSamD"); }

		( $ntConv, $gDummy ) = select_gnumap_option( $gRead[$rCONV], $k, 0 );
		$cmdLine = "$cmdLine\ngrep -P '\\t$gnuSamStrand:$ntConv\' $inSamD/$inSamfn > $outSamD/$inSamfn\n";
	}

	$cmdLine = "$cmdLine\nrm -rf $inSamD/$inSamfn\n";
	return ($cmdLine);
}

sub selectTxtBetween2Lines {
	my ( $splitIdx, $cmd2qsub, $line1, $line2, $dfn ) = @_;

	my $dfn_j = "${dfn}_${splitIdx}";

	unless ( -e $dfn_j ) {
		$cmd2qsub = "$cmd2qsub\nawk \'BEGIN {min=$line1;max=$line2}{if (NR>=min){if(NR<=max) print}}\' $dfn > $dfn_j\n";
	}
	return ( $cmd2qsub, $dfn_j );
}

sub get_unampped_reads {

	my ( $samDprefix, $alnDfn ) = @_;    #get a hit sam file

	#prepare the location to store unmapped reads
	my $gnuUnmapD = "$samDprefix/$gOut[$oMAPPED][2]";
	unless ( -d $gnuUnmapD ) {
		system("mkdir -p $gnuUnmapD");
	}

	my @OrgRead = ( "$gOut[$oPI][0].$gRead[$rFEXT]", "X" );
	if ( $gRead[$rPER] == 2 ) {
		@OrgRead = ( "$gOut[$oPI][0].$gRead[$rFEXT]", "$gOut[$oPI][1].$gRead[$rFEXT]" );
	}

	my @UnmappedPe = cmd_extractUnmappedPeRead(\@OrgRead, $alnDfn, $gnuUnmapD);
	
	my $isEmpty=0;
	if ($UnmappedPe[0] eq "X"){
		$isEmpty=1;
	}

	return ($isEmpty,$gnuUnmapD);
}

sub cmd_extractUnmappedPeRead {
	my ($OrgRead, $gnuSamConcoDf, $gnuUnmapD ) = @_;

	my $readD = $gRead[$rDIR];
	if (defined $gRead[$rDIR_QC]){
		$readD = $gRead[$rDIR_QC];
	}
	my $fastaF = 0;
	if ($gRead[$rPHRED] == 0) {
		$fastaF = 1;
	}
	
	my $cmdLine = "$gBIN[$bUNMAP2FA][$buBIN] $readD $fastaF $OrgRead->[0] $OrgRead->[1] $gnuSamConcoDf $gnuUnmapD\n";
	
	system($cmdLine);

	my @UnmappedPe = ( "X", "X" );
	$UnmappedPe[0] = "$gnuUnmapD/$gOut[$oPI][0].$gBIN[$bBLAT][$bbREAD_FEXT]";
	
	my $filesize = -s $UnmappedPe[0];
	
	my $emptyUnmappedF = 0;
	if ($filesize < 1){
		$UnmappedPe[0]="X";
	}

	if ( $OrgRead->[1] ne "X" ) {
		$UnmappedPe[1] = "$gnuUnmapD/$gOut[$oPI][1].$gBIN[$bBLAT][$bbREAD_FEXT]";
	}
		
	return @UnmappedPe;
}

#==================================
sub cmd_trim_read {
	my ( $progCmd, $inD, $InFn, $outD, $handleInAfter ) = @_;

	my ( @valOut, @UnpOut ) = ();
	my ( $i, $prevD );

	my @adaptor = ( "TGGAATTCTCGGGTGCCAAGG", "GATCGTCGGACTGTAGAACTCTGAAC" );

	my @InDfn  = ( "$inD/$InFn->[0]",  "$inD/$InFn->[1]" );
	my @OutDfn = ( "$outD/$InFn->[0]", "$outD/$InFn->[1]" );

	$progCmd = "$progCmd\n$gBIN[$bQC][$bqCUTADPT] --phred64 -a $adaptor[0] -a2 $adaptor[1] --paired --retain_unpaired $InDfn[0] $InDfn[1]\n"; #default overlap = 1

	#this is a hard-coded output from this trim tool
	$valOut[0] = $InDfn[0] . "_val_1.fq";
	$valOut[1] = $InDfn[1] . "_val_2.fq";
	$UnpOut[0] = $InDfn[0] . "_unpaired_1.fq";
	$UnpOut[1] = $InDfn[1] . "_unpaired_2.fq";

	my $unpOutBookMark = "$UnpOut[0]\.bookmark";

	for ( $i = 0 ; $i < 2 ; $i++ ) {
		( $gDummy, $progCmd, $gDummy ) =
		  cmd_zip_util( $progCmd, $InDfn[$i], $handleInAfter );

		$progCmd = "$progCmd\ncat $UnpOut[$i] >> $valOut[$i]";

		#system("tail -4 $UnpOut[$i] >> $unpOutBookMark");
		$progCmd = "$progCmd\ntail -4 $UnpOut[$i] >> $unpOutBookMark";
		$progCmd = "$progCmd\nrm -f $UnpOut[$i]";
		$progCmd = "$progCmd\nmv $valOut[$i] $OutDfn[$i]\n";
	}
	$prevD = 0;
	return ( $prevD, $progCmd, @OutDfn );
}

#====================
sub select_gnumap_option {
	my ( $readSpecial, $readSpecial_strand, $sense ) = @_;
	my ($WATSON,$CRICK) = (0,1);
	my $gnuNtConvOpt = "";
	my $idxSuffix    = "na";

	#to prepare gnumap (nt converted if any) index
	if ( $readSpecial eq "bs" ) { #bisulfited watson
		if ($readSpecial_strand == $WATSON) {
			$gnuNtConvOpt = "-b";
			$idxSuffix    = "ct";
		}
		else {
			$gnuNtConvOpt = "-B";
			$idxSuffix    = "ga";
		}
	}
	elsif ( $readSpecial eq "a2i") { #a to i in watson
		if ($readSpecial_strand == $WATSON) {
			$gnuNtConvOpt = "-d";
			$idxSuffix    = "ag";
		}
		else {
			$gnuNtConvOpt = "-D";
			$idxSuffix    = "tc";
		}
	}
		
	if ( $sense == $UP ) {
		$gnuNtConvOpt = "$gnuNtConvOpt --up_strand";
	}
	elsif ( $sense == $DN ) {
		$gnuNtConvOpt = "$gnuNtConvOpt --down_strand";
	}
	else {    #both
		$gnuNtConvOpt = "$gnuNtConvOpt";
	}

	if ( $gRead[$rPHRED] == 64 ) {
		$gnuNtConvOpt = "$gnuNtConvOpt --illumina";
	}
	elsif ( $gRead[$rPHRED] == 0) {
		$gnuNtConvOpt = "$gnuNtConvOpt --noqual";
	}

	if ( $gBIN[$bGNU][$bgGNU_SW] == 1 ) {
		$gnuNtConvOpt = "$gnuNtConvOpt --sw";
	}

	return ( $idxSuffix, $gnuNtConvOpt );
}

#=============================================
sub run_gnumap {
	my ( $refDfn, $Reads ) = @_;

	#consider all cases
	#for example, (normal,nt_conv)x(per,ser)x(unstrand,strand)

	my (
		$rBase,  $p,        $pi,    $s,     $samD,
		$k,      $cmd1,     $jobId, $prevD, $targetReadD,
		$samDfn, $aStrandD, $validF
	);

	#gnumap only maps to genome sequence only for now.
	my $gnuOutPrefix =  "$gOut[$oALIGND]/$gBIN[$bGNU][$bgOPRFX]/$gOut[$oREFTYPE][0]";

	for ( $p = 0 ; $p < $gRead[$rPER] ; $p++ )
	{    # for each paired read (for single end, only _1.fq will be processed)
		$pi = $p + 1;

		for ( $s = 0 ; $s < $gRef[$gSTRAND_DELIM] ; $s++ )
		{    #for each strand, wt/cr in case nt_conv. O.w., only use wt
			
			$samD = "$gnuOutPrefix/$gOut[$oSTRAND][$s]";
			unless ( -d $samD ) { system("mkdir -p $samD"); }
			#print "$gSens2map[$p][$s]\n"; #!debug
			$samDfn = cmd_gnumap( $refDfn, $gRef[$gGNUIDX][$s],$samD, $Reads->[$p], $pi, $s, $gSens2map[$p][$s] ); #note that $gRef[$gGNUIDX][$s] from build_ref_index_gnumap() and $gSens2map[$p][$s] guides which strand(sense,antisense) gnumap should align
			
		}
	}
	return ($gnuOutPrefix);
}

#===========================================================================
sub cmd_gnumap {

	my ( $refDfn, $refIdx, $samD, $read, $pairIdx, $strandIdx, $senseIdx ) = @_;

	my ( $validF, $cmdLine_resume, $phredStr, $aliDf, $aliDfx, $readRe,
		$samDf, $samDfn );

	
	my ( $ntConv, $gnuNtConvOpt ) = select_gnumap_option( $gRead[$rCONV], $strandIdx, $senseIdx );
	$aliDfx = "$samD/$pairIdx";
	my $cmdLine = "export LD_LIBRARY_PATH=$gBIN[$bGNU][$bgBIN_LIB]:\$LD_LIBRARY_PATH\n";
	$cmdLine = "$cmdLine\n$gBIN[$bGNU][$bgBIN] -g $refDfn --read=$refIdx -o $aliDfx -m $gBIN[$bGNU][$bgGNU_KMER] -s $gBIN[$bGNU][$bgGNU_SLIDE] -c $gTHREADS -h $gBIN[$bGNU][$bgGNU_MAX_HASH] -t 25 -a $gBIN[$bGNU][$bgGNU_ACC] $gnuNtConvOpt -y $gBIN[$bGNU][$bgGNU_TOPK_HASH] $read\n";

	$samDfn = "$aliDfx.$gOut[$oMEXT][0]";
	if ($gDebug==1 && -e $samDfn){
		print STDOUT "use prev results ($samDfn)...\n";
	}
	else{
		#print "$cmdLine\n";    #!debug
		system($cmdLine);
	}
	return $samDfn;
}

#=================
sub resume_sam {
	my ( $samDfn, $readDfn ) = @_;
	my $rndTag        = int( rand(100000) );
	my $readDfnResume = "$readDfn-$rndTag-resumed";
	my ( $cmd, $x, $lastRead, $lineN );
	my $filesize = -s $samDfn;

	if ( $filesize > 0 ) {

		#first fix prev samDfn
		$cmd = "sed -i \'\$d\' $samDfn";

		#print STDERR "d:$cmd\n";
		system($cmd);
		sleep(1);

		#identify the last read in samDfn
		$cmd = "tail -n 1 $samDfn | awk -F\"\\t\" \'\{ print \$1 \}\'";

		#print STDERR "d:$cmd\n";
		$x = `$cmd`;

		chomp($x);

		#sleep(1);
		$lastRead = $x;

		#print STDERR "d:$x\n";
		#extract $readDfn from the last read
		$cmd = "grep -P -n -m 1 \'$lastRead\' $readDfn | awk -F: \'\{print \$1\}\'";

		#print STDERR "d:$cmd\n";
		$x = `$cmd`;

		chomp($x);
		$lineN = $x + 4;

		#sleep(1);
	}
	else {
		$lineN = 1;
	}

	my $progCmd = "tail -n \+$lineN $readDfn > $readDfnResume\n";

	#print STDERR "d:$progCmd\n";

	return ( $progCmd, $readDfnResume );
}

#==========================================
sub muxSam {
	my ($samDprefix) = @_;
	my ( $prevD, $catPat, $outCatDfn);

	my $outCatD = "$samDprefix/$gOut[$oMAPPED][0]";    #simply merged
	unless ( -d $outCatD ) { system("mkdir -p $outCatD"); }

	if ( $gRead[$rPER] == 2 ) {
		$catPat = "$samDprefix/*/*/$gOut[$oPI][3]/\*/$gOut[$oPI][3]\.$gOut[$oMEXT][0]";
		$outCatDfn = "$outCatD/$gOut[$oPI][3]\.$gOut[$oMEXT][0]";    #TODO: chr* pat limitation
	}
	else {     #simply collect sam files under strand directories
		$catPat    = "$samDprefix/*/$gOut[$oPI][0]\.$gOut[$oMEXT][0]";
		$outCatDfn = "$outCatD/$gOut[$oPI][2]\.$gOut[$oMEXT][0]";
	}
	
	my $handleInAfter="keep";
	if ($gDebug==0){
		$handleInAfter="delete";
	}
	( $prevD, $outCatDfn ) = cmd_bigCat_inPat( $catPat, $outCatDfn, $handleInAfter ); #!debug
	
	return $outCatDfn;
}

#=============================
sub demuxPeSam {
	my ($Chrs, $samDprefix, $subD2opt1, $subD2opt2) = @_;

	my ( $i, $j, $sam, $outD, $samD, $samFn, $chrD, $demuxSam, $demuxSamD, $s, $s1, $s2, $t, $samFnSz, $cmdLine);
	my (@SamFiles) = ();
	my (@SubD)=($subD2opt1,$subD2opt2);

	#$s1 = $gSrchDim[$sSTRAND][$ST];
	#$s2 = $gRead[$rPER];

	my $aliExtPat = "\.$gOut[$oMEXT][0]";
	my $aliExt    = ".$gOut[$oMEXT][0]";
	my @SamFlag   = ( [ 0, 256 ], [ 16, 272 ] );

	#my @Strand = ("up","dn");

	my $pairId = 1;

	#my $patRef = "\.$gRef[$gCHR_FEXT]";

	#print STDERR "demuxing sam...\n"; #!debug
	for ( $s = 0 ; $s < 2 ; $s++ ) {  #for each strand to filter
		if ($SubD[$s] eq 'X') {next;}
		
		$samD = "$samDprefix/$SubD[$s]";

		for ( $i = 0 ; $i < 2 ; $i++ ) {                  #for each pair
			$pairId = $i + 1;
			$samFn  = "$samD/${pairId}${aliExt}";         #input file
			if (-e $samFn) {
				$samFnSz = -s $samFn;
				if ($samFnSz<1) {next;}
			}
			else {next;}
			
			foreach $chrD ( @{$Chrs} ) {                  #for each chr
				for ( $j = 0 ; $j < 2 ; $j++ )
				{    #for each sense ($sense,antisense|up,dn)
					$demuxSamD = "$samD/$chrD/$gOut[$oPCR][$j]";
					unless ( -d $demuxSamD ) {
						system("mkdir -p $demuxSamD");
					}
					$demuxSam = "$demuxSamD/${pairId}${aliExt}";
					$cmdLine = "grep -P \'^\\S+\\t($SamFlag[$j][0]|$SamFlag[$j][1])\\t\\b$chrD\\b\' $samFn > $demuxSam\n";
					system($cmdLine);
				}
			}
		}
	}
	#print "done\n"; #!debug
}

#==================================
sub cleanup_files {
	my ($expBaseD,$subD1,$subD2) = @_;
	
	my $dir2unlink;
	if ($subD1 ne $NAs){
		$dir2unlink = "$expBaseD/$subD1";
		if (-e $dir2unlink) {
			system("rm -rf $dir2unlink");
		}
	}
	if (defined($subD2) && ($subD2 ne $NAs)){
		$dir2unlink = "$expBaseD/$subD2";
		if (-e $dir2unlink) {
			system("rm -rf $dir2unlink");
		}
	}
}
#==========================================
sub cmd_bigCat_inPat {
	my ( $smallFishes, $bigFish, $handleInAfter ) = @_;

	my $progCmd = "cat $smallFishes > $bigFish\n";
	system($progCmd);
	$progCmd="";
	( $gDummy, $progCmd, $gDummy ) = cmd_zip_util( $progCmd, $smallFishes, $handleInAfter );
	$progCmd = "$progCmd\nsleep 1\n";
	system($progCmd);
	my $prevD = 0;
	if ( -e $bigFish ) { $prevD = 1; }

	return ( $prevD, $bigFish );
}

#===========================================================================
sub perSam_core {
	my ( $Chrs, $samDprefix,  $subD2opt1, $subD2opt2, $mapper ) = @_;

	my @AstrandDfn=[('X','X'),('X','X')];
	my (@SubD)=($subD2opt1,$subD2opt2);
	
	my ($s,$j,$k,$peAliD,$peAliDfn,$cmd1,$chrBase,$cmdLine,$prevD,$pairId,$chrSamD,$dummy);    #todo:cleanup

	# 	my $patRead = "_$gSrchDim[$oPI][$ST].$gRead[$rFILE_EXT]";
	# 	my $patRef = ".$gRef[$gCHR_FEXT]";

	my $fopInAfter = "keep";
	if ($gDebug==0){
		$fopInAfter = "delete"; #!debug
	}
	
	for ($s = 0 ; $s < 2 ; $s++) {    #for each strand
		if ($SubD[$s] eq 'X') {next;}
		
		foreach $chrBase ( @{$Chrs} ) {                     #for each chr

			#$cmdLine="\n";
			$chrSamD = "$samDprefix/$SubD[$s]/${chrBase}";

			for ($j = 0 ; $j < 2 ; $j++) {                #for each pair
				$pairId = $j + 1;
				for ($k = 0 ; $k < 2 ; $k++)#for each sense i.e., up/dn; it should exist already!
				{    
					$AstrandDfn[$j][$k] = "$chrSamD/$gOut[$oPCR][$k]/${pairId}.$gOut[$oMEXT][0]";
				}
			}

			#merge 2 sam
			$peAliD = "$chrSamD/$gOut[$oPI][3]/$gOut[$oPER][0]";
			unless ( -d $peAliD ) { system("mkdir -p $peAliD"); }

			$peAliDfn = "$peAliD/$gOut[$oPI][3].$gOut[$oMEXT][0]";

			cmd_merge2peSam($AstrandDfn[0][0],$AstrandDfn[1][1],$peAliDfn,1,$gRead[$rREAD_TYPE],$gRead[$rPER_AVE_DIST],$gRead[$rMAX_FRAG_LEN], $mapper, $fopInAfter);


			$peAliD = "$chrSamD/$gOut[$oPI][3]/$gOut[$oPER][1]";
			unless ( -d $peAliD ) { system("mkdir -p $peAliD"); }

			$peAliDfn = "$peAliD/$gOut[$oPI][3].$gOut[$oMEXT][0]";

			cmd_merge2peSam($AstrandDfn[0][1],$AstrandDfn[1][0],$peAliDfn,2,$gRead[$rREAD_TYPE],$gRead[$rPER_AVE_DIST],$gRead[$rMAX_FRAG_LEN], $mapper, $fopInAfter);

		}
	}
}

#==========
sub cmd_grep_awk_get_list_with_pat {

	my ( $file, $pat ) = @_;
	my $rndTag  = int( rand(100000) );
	my $tmpFile = "$gOut[$oTEMPD]/chrList-$rndTag.txt";

	#	print STDERR "$tmpFile\n";
	my $cmd = "grep -iP \'$pat\' $file \| sed \'s/$pat//g\' > $tmpFile";

	#	print STDERR "$cmd\n";

	system($cmd);

	my $x = `cat $tmpFile`;
	my @List4Pat = split( '\n', $x );

	system("rm -rf $tmpFile");
	return (@List4Pat);
}

#====================
sub check_file_sanity {
	my ($file2exam) = @_;
	my ($fileSize,$status);

	if (-e $file2exam) {
		$fileSize = -s $file2exam;
		$status   = $VALID;
		if ( $fileSize == 0 ) {
			$status = $EMPTY;
		}
	}
	else {
		$status = $NOT_EXIST;
	}
	return ($status);
}

#================================================================
sub cmd_merge2peSam {
	my ( $samA, $samB, $peDf, $pairAlignMode, $refType, $avePerDist,
		$maxFragLen, $mapper, $handleInAfter )
	  = @_;
	
	my $progCmd = "";
	
	my $allowDiscordant = 0;
	if ($mapper eq "blat"){
		$allowDiscordant = 1;
	}
	
	$progCmd ="$gBIN[$bPER][$bpBIN] $samA $samB $peDf $refType $pairAlignMode $avePerDist $maxFragLen $allowDiscordant $gOut[$oTEMPD] $handleInAfter";
	#print STDOUT "$progCmd\n"; #!debug
	system($progCmd);
	
}

#==========================================================================
sub cmd_resolve_ambig_sam2 {
	my ( $mapper, $per, $ambSam, $bestSam, $handleInAfter ) = @_;

	my $mapperId = 1;
	if ( $mapper eq "blat" ) { $mapperId = 2; }

	my $gnuNtConvF = 0;
	if ( $gRead[$rCONV] ne 'X' ) {
		$gnuNtConvF = 1;
	}

	my $cmdLine ="$gBIN[$bBEST][$b1BIN] $mapperId $per $gnuNtConvF $ambSam $bestSam $gOut[$oTEMPD] $gBIN[$bBEST][$b1NUM_AMB]\n";
	
	if ($gDebug==1) {
		print STDOUT "$cmdLine\n";
	}
	system($cmdLine);
	$cmdLine ="";
	( $gDummy, $cmdLine, $gDummy ) = cmd_zip_util( $cmdLine, $ambSam, $handleInAfter );
	system($cmdLine);
}

#=========================================================================
sub getFaSplitDelim {

	my ( $faDfn, $nSplits ) = @_;
	my $FASTA_LINE_UNIT = 2;
	my ( $i, $M );

	#my $tmp = "wc -l $faDfn \| awk -F\" \" \'{ print \$1 }\'";
	#print $tmp;

	my $x = `wc -l $faDfn \| awk -F\" \" \'{ print \$1 }\'`;
	chomp($x);
	my $N        = $x + 0;                 #num of lines
	my $lineUnit = $FASTA_LINE_UNIT;
	my $M1       = int( $N / $nSplits );

	$M = $M1;
	if ( $M1 % 2 == 1 ) {
		$M = $M1 + 1;
	}

	my $S1 = $nSplits - 1;
	my ( @DelimSTs, @DelimEnds ) = ();

	for ( $i = 0 ; $i < $S1 ; $i++ ) {
		$DelimSTs[$i] = $M * $i + 1;
		$DelimEnds[$i] = $M * ( $i + 1 );
	}
	$DelimSTs[$S1]  = $M * $S1 + 1;
	$DelimEnds[$S1] = $N;

	return ( \@DelimSTs, \@DelimEnds );
}


#===========================================================================
sub listFilesInDir {
	my ( $dr, $prefixDirF, $readUp2pat, $prefixFileF ) = @_;

	my @fList = ();
	my ( $fn, $pDr );
	#print STDERR "reading $dr...\n";
	opendir( $pDr, $dr );
	while ( $fn = readdir($pDr) ) {
		if ( $fn =~ /^\.+/ ) { next; }
		if (   ( $fn =~ /(\S+)$readUp2pat$/ )
			|| ( $fn =~ /(\S+)$readUp2pat\.gz$/ ) )
		{
			if ( $prefixDirF == 1 ) {
				if ( $prefixFileF == 1 ) {
					push( @fList, "$dr/$1" );
				}
				else {
					push( @fList, "$dr/$fn" );
				}
			}
			else {
				if ( $prefixFileF == 1 ) {
					push( @fList, $1 );
				}
				else {
					push( @fList, $fn );
				}
			}
		}
	}
	closedir($pDr);
	my $listSz = $#fList + 1;

	if ( $listSz < 1 ) {
		print "hey, no valid files in this directory, $dr !\n";
	}

	return @fList;
}

sub cmd_zip_util {
	##---------------------------
	## Be aware of what you are doing This is related with delete job
	## make sure that the progCmd associated with file2remove does not use any background process
	##---------------------------
	my ( $progCmd, $file2handle, $opSel ) = @_;

	my ($outFn,$file2handleSz);
# 	if (-e $file2handle) {
# 		$file2handleSz = -s $file2handle;
# 		if ($file2handleSz < 1){return ("",$progCmd,"");}
# 	}
# 	else{return ("",$progCmd,"");}
	
	if ( $opSel eq "delete" ) {
		$progCmd = "$progCmd\nsleep 1\nrm -rf $file2handle\n";
		$outFn   = "X";
	}
	elsif ( $opSel eq "zip" ) {
		$progCmd = "$progCmd\nsleep 1\ngzip -f $file2handle\n";
		$outFn   = "$file2handle.gz";
	}
	elsif ( $opSel eq "unzip-n-keep" ) {
		$outFn = substr( $file2handle, 0, -3 );
		$progCmd = "$progCmd\nsleep 1\ngunzip -f -c $file2handle > $outFn\n";
	}
	elsif ( $opSel eq "unzip-n-delete" ) {
		$outFn = substr( $file2handle, 0, -3 );
		$progCmd = "$progCmd\nsleep 1\ngunzip -f $file2handle\n";
	}
	my $prevD = 0;

	return ( $prevD, $progCmd, $outFn );
}

sub flatten {
	my @work = @_;
	my @result;
	while (@work) {
		my $next = shift @work;
		if ( ref($next) eq 'ARRAY' ) {
			unshift @work, @$next;
		}
		else {
			push @result, $next;
		}
	}
	return @result;
}

#==============
sub messageBox {
	my ( $msg, $notice ) = @_;
	print STDERR "===========================================\n";
	print STDERR "$msg submitted\n";
	print STDERR "follow $notice after the submitted jobs are all completed\n";
	print STDERR "===========================================\n";
}

#==============
sub gnumap_pipeline_logo {
	my ($tick,$runPoint,$msg2display,$installD) = @_;
	system "clear";
	
	my ($tock,$elapsed,$elapsedStr);
	$tock = time;
	if ($runPoint == 0) {
		$msg2display = "In progress ...";
		$elapsedStr = "-elapsed time: 0.0";
	}
	elsif ($runPoint == 1) {
		$msg2display = "$msg2display\ncompleted!";
		$elapsed = $tock - $tick;
		$elapsedStr = "-elapsed time: $elapsed";
	}
	print STDERR "\n\n===========================================\n";
	print STDERR "GNUMAP_PIPELINE(v0.1)\n";
	print STDERR "-http://jlab.bu.edu\n";
	print STDERR "-cjhong\@bu.edu\n";
	print STDERR "-running at $installD\n";
	system "date";
	print STDERR "$elapsedStr\n";
	print STDERR "===========================================\n";
	print STDERR "$msg2display\n";
	return $tock;
}

#===================
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

#=======================
sub parse_arguments {
	my (
		$ref_genome, $ref_transc, $pair_1,       $pair_2,
		$single,     $qcontrol,       $adap,         $adap2,
		$lib_type,   $read_type,      $phred_offset, $per_ave_dist,
		$nt_conv,    $map_quality, $num_amb, $mpi, $num_threads,  $outD, $debug_mode
	) = @_;
	
	my ($read1,$read2,$read0,$dummy);
	my $installD = $ENV{'GNUMAPS'};
	unless (defined $installD) {
		die "please set environment varialble, export GNUMAPS=gnumap_wrapper_install_directory_path !\n";
	}

	$gDebug = $debug_mode;
	 #1. check if ref genomes exist
	unless (defined $ref_genome) {
		die "specify reference genome.";
	}

	unless (-e $ref_genome) {
		die "$ref_genome does not exist!";
	}
	
	$ref_genome = abs_path($ref_genome);
	
	
	$gRef[$gGNREF] = $ref_genome;
	( $dummy, $gRef[$gGNREF_BASE], $gRef[$gGNREF_FEXT] ) = parse_filename($ref_genome);


	#3. check reads conditions
	if ( ( !defined($pair_1) || !defined($pair_2) ) && !defined($single) ) {
		die "check input reads";
	}

	#3.1 check if read is per or ser
	#it will be used a delimiter of mapping index
	$gRead[$rPER] = 2;
	my $frag_len = 1500; #TODO
	if (defined $single) {
		$gRead[$rPER] = 1;
		$pair_1 = abs_path($single);
	  $read0 = $pair_1;    #it will be used to get a read base
	  unless (-e $read0) {
			die "there is no input file, $read0!";
		}
	} 
	else {
		$pair_1 = abs_path($pair_1);
		$pair_2 = abs_path($pair_2);
		unless (-e $pair_2) {
			die "there is no input file, $pair_2!";
		}
		
		$read0 = $pair_1;
		$gRead[$rPER_AVE_DIST] = 0;
		if (defined $per_ave_dist) {
			$gRead[$rPER_AVE_DIST] = int($per_ave_dist) + 0;
		} 
		$gRead[$rMAX_FRAG_LEN] = 1500;
		if (defined $frag_len) {
			$gRead[$rMAX_FRAG_LEN] = int($frag_len);
		}
	}

	#3.1b fill read file, base, and file extension etc
	($dummy,$gRead[$rFILE_BASE],$gRead[$rFEXT]) = parse_filename($read0);

	#5. prepare output directories
	$outD = abs_path($outD);
	$gOut[$oWORKD] = $outD;
	

	$gRead[$rDIR] = "$outD/$gRead[$rFILE_BASE]"; #we create soft link copy of the original input to exp output directory so that the original dataset won't be intact!
	
	unless ( -d $gRead[$rDIR] ) {
		system("mkdir -p $gRead[$rDIR]");
	}

	$read1 = "$gRead[$rDIR]/1.$gRead[$rFEXT]";
	unless (-l $read1) {
		system("ln -s $pair_1 $read1");
	}
	
	if ( $gRead[$rPER] == 2 ) {
		$read2 = "$gRead[$rDIR]/2.$gRead[$rFEXT]";
		unless (-l $read2) {
			system("ln -s $pair_2 $read2");
		}
	}

	$gOut[$oALIGND] = "$gRead[$rDIR]/aligned";
	if ( -d $gOut[$oALIGND] ) {
		if ($gDebug==0){
			#system("rm -rf $gOut[$oALIGND]"); #!debug
		}
	}
	unless ( -d $gOut[$oALIGND] ) {
		#system("mkdir -p $gOut[$oALIGND]");
	}

	#3.2 check if read is dna or rna
	$gRef[$gREF_DELIM]  = 1;
	$gRead[$rREAD_TYPE] = $DNA;
	my $tRefD = 'X';
	if ((defined $read_type) && ($read_type eq "rna" )) {
		$gRead[$rREAD_TYPE] = $RNA;
		#check if transcriptome ref exists
		if (defined $ref_transc) {
			$ref_transc = abs_path($ref_transc);
			$gRef[$gTNREF] = $ref_transc;
			( $dummy, $gRef[$gTNREF_BASE], $gRef[$gTNREF_FEXT] ) = parse_filename($ref_transc);
			$gRef[$gREF_DELIM]  = 2;
		}
	}

	if ($nt_conv eq $NAs) {
		$gRef[$gSTRAND_DELIM] = 1;
		$gRead[$rCONV] = $NAs;
	}
	elsif (($nt_conv eq "bs") || ($nt_conv eq "a2i")) {
		$gRead[$rCONV]  = $nt_conv;
		$gRef[$gSTRAND_DELIM] = 2;
	}
	else {
		die "check nt_conv in the input argument!\n";
	}
	
	#3.3 check if it is stranded or not
	@gSens2map = ( [ $ALL, $ALL ], [ $ALL, $ALL ]);
	
	$gRead[$rLIB_TYPE] = $UNSTRAND;
	if (defined $lib_type && ($lib_type ne "unstrand")) {
		$gRead[$rLIB_TYPE] = $STRAND1;
		@gSens2map = ( [ $UP, $DN ], [ $DN, $UP ]);
	  if ($lib_type eq "cr1") {
			@gSens2map = ( [ $DN, $UP ], [ $UP, $DN ]);
		  $gRead[$rLIB_TYPE] = $STRAND2;
		}
	} 
	else { $lib_type = "unstrand"; }

	my $depD = "$installD/dep";
	my $scriptD = "$installD/scripts";

	#4. check if quality control is enabled (TODO)
	if ( (defined $qcontrol) && $qcontrol == 1 ) {
		$gRead[$rDIR_QC] = "$gRead[$rDIR]/qc";
		unless ( -d $gRead[$rDIR_QC] ) {
			system("mkdir -p $gRead[$rDIR_QC]");
		}
		# about trimmer
		if (defined $adap) {
			$gBIN[$bQC][$bqADAP]         = $adap;
			  $gBIN[$bQC][$bqCUTADPT] = "$depD/cutadapt/bin/cutadpt";
		  }

		  if (defined $adap2) {
			$gBIN[$bQC][$bqADAP2] = $adap2;
		  } $gBIN[$bQC][$bqPRTSEQ] = "$depD/prinseq-lite/prinseq-lite.pl";
	}

#$gBIN[$bQC][$oOPRFX] = "qc"; # once you perform quality control, gRead will be updated too!

	#check phredoffset
	$gRead[$rPHRED] = 33;
	if (defined $phred_offset) {
		$gRead[$rPHRED] = $phred_offset;
		if ($gRead[$rFEXT] eq "") {
			$gRead[$rFEXT] = "fa";
		}
	}
	elsif ($gRead[$rFEXT] eq "") {
		$gRead[$rFEXT] = "fq";
	}


	$gOut[$oGNU_IDX_FEXT] = [ ( "gnu", ) ];
	$gOut[$oREFTYPE] = [ ( "gen", "tgen", ) ];

	$gOut[$oTEMPD] = "$gOut[$oWORKD]/tmp";
	unless ( -d $gOut[$oTEMPD] ) {
		system("mkdir -p $gOut[$oTEMPD]");
	}

	#----------------------------------
	# about gnumap options
	my $appendMpiStr="";
	if ( $mpi == 1 ) {
		$appendMpiStr="mpirun -np 1 ";
	}
	
	$gBIN[$bGNU][$bgBIN]      = "${appendMpiStr}${depD}/gnumap_MPI_opt/bin/gnumap-stl"; 
	print STDOUT "$gBIN[$bGNU][$bgBIN]\n"; #debug
	
	$gBIN[$bGNU][$bgOPRFX]    = "gnu";
	
	$gBIN[$bGNU][$bgGNU_ACC]  = 0.95;
	$gBIN[$bGNU][$bgGNU_SLIDE] = 1;
	$gBIN[$bGNU][$bgGNU_TOPK_HASH] = 3;
	if ($map_quality eq "sensitive") {
		$gBIN[$bGNU][$bgGNU_ACC]  = 0.93;
		#$gBIN[$bGNU][$bgGNU_SLIDE] = 1; 
		$gBIN[$bGNU][$bgGNU_TOPK_HASH] = 4;
	}
	elsif ($map_quality eq "fast") {
		$gBIN[$bGNU][$bgGNU_ACC]  = 0.97;
		#$gBIN[$bGNU][$bgGNU_SLIDE] = 3; 
		$gBIN[$bGNU][$bgGNU_TOPK_HASH] = 2;
	}
	
	
	$gBIN[$bGNU][$bgGNU_KMER] = 17;

	#$gBIN[$bGNU][$oGNU_KMER_TRAN] = $gBIN[$bGNU][$oGNU_KMER] + 2;
	$gBIN[$bGNU][$bgGNU_MAX_HASH] = 500;

	#$gBIN[$bGNU][$oGNU_MAX_HASH_TRAN] = 600;
	$gBIN[$bGNU][$bgBIN_LIB]   = "$depD/gnumap_MPI_opt/lib/lib";
	
	$gBIN[$bGNU][$bgGNU_SW]    = 0;

	#-----------------------------------
	# about merge_2sam
	$gBIN[$bPER][$bpBIN]  = "$scriptD/merge_2sam_v4.pl";
	$gBIN[$bPER][$bpPRFX] = "pe";

	#-----------------------------------
	# about pickup_best
	$gBIN[$bBEST][$b1BIN]  = "$scriptD/resolve_amb_sam.pl";
	$gBIN[$bBEST][$b1PRFX] = "best";
	$gBIN[$bBEST][$b1NUM_AMB] = $num_amb; #only best unique allowed

	#-----------------------------------

	# about gnumap_index
	$gBIN[$bGNU][$bgBIN]      = "${appendMpiStr}${depD}/gnumap_MPI_opt/bin/gnumap-stl";
	$gBIN[$bGNU][$bgiIDX_FEXT] = "gnu";
	$gBIN[$bGNU][$bgiOPRFX]    = "idxGnu";
	$gBIN[$bGNU][$bgiNULLREAD] = "$scriptD/gnuIdxRead.fq";
	$gBIN[$bGNU][$bgiREVCOMP_SAM]  = "$scriptD/revCompSam.pl";

	#-----------------------------------
	# about sam2sgr
	$gBIN[$bGNUSGR][$bgsBIN]     = "${appendMpiStr}${depD}/gnumap_MPI_opt/bin/sam2sgr";
	$gBIN[$bGNUSGR][$bgsOPRFX]   = "gmp";
	$gBIN[$bGNUSGR][$bgsPOSTPROC] = "$scriptD/filterSam2GmpReport.pl";

	#----------------------------------
	# about blat options

	$gBIN[$bBLAT][$bbBIN]        = "$depD/blat/blat";
	$gBIN[$bBLAT][$bbOPRFX]      = "blat";
	$gBIN[$bBLAT][$bbPERIDN]     = 93;
	$gBIN[$bBLAT][$bbTILESZ]     = 11;
	$gBIN[$bBLAT][$bbSTEPSZ]     = 11;
	$gBIN[$bBLAT][$bbREPMAT]     = 256;
	$gBIN[$bBLAT][$bbMAX_INTRON] = 50000;
	$gBIN[$bBLAT][$bbREAD_FEXT]  = "fa";
	

	#-----------------------------------
	# about UseqSTP
	$gBIN[$bSTP][$bstBIN]     = "$depD/USeq/Apps/SamTranscriptomeParser";
	$gBIN[$bSTP][$bstOPRFX]   = "tgen";
	$gBIN[$bSTP][$bstBIN_LIB] = "$depD/USeq/Apps";

	#-----------------------------------
	$gBIN[$bSAMT][$bstBIN]   = "$depD/samtools/samtools";
	$gBIN[$bSAMT][$bstOPRFX] = "output";

	#-----------------------------------
	# about psl2sam
	$gBIN[$bPSL2SAM][$bpsBIN]   = "$scriptD/psl2sam_v8.pl";
	$gBIN[$bPSL2SAM][$bpsOPRFX] = "sam";

	#-----------------------------------
	# about get_unmapped_read
	$gBIN[$bUNMAP2FA][$buBIN]   = "$scriptD/get_unmapped_reads_v2.pl";
	$gBIN[$bUNMAP2FA][$buOPRFX] = "unmapped";

	#-----------------------------------
	# about fato2bit
	$gBIN[$b2BIT][$b2BIN]      = "$depD/blat/faToTwoBit";
	$gBIN[$b2BIT][$b2OPRFX]    = "idxBlat";
	$gBIN[$b2BIT][$b2IDX_FEXT] = "2bit";

	#-----------------------------------
	$gOut[$oPI]     = [ ( "1",      "2",      "se", "pe", "2flp" ) ];
	$gOut[$oMAPPED] = [ ( "merged", "mapped", "unmapped" ) ];
	$gOut[$oSTRAND] = [ ( "watson", "crick", "both") ]; #this indicates nt conversion strand if user enables nt_convF

	$gOut[$oPCR]  = [ ( "up",     "dn", ) ];
	$gOut[$oPER]  = [ ( "up1dn2", "dn1up2", ) ];
	$gOut[$oMEXT] = [ ( "sam",    "psl", ) ];
	$gOut[$oZIP]  = [ ( "gz", ) ];
	$gOut[$oBAM]  = [ ( "bam", ) ];

	$gTHREADS = 1;
	if (defined $num_threads) {
		$gTHREADS = $num_threads;
	} 
	
	return ( $read1, $read2, $installD);
}


#====================
# sub sam2bam {
# 
# 	my ($wStep) = @_;
# 
# 	my (
# 		$s,     $samBase, $samD,    $samDfn, $bamFn, $samStFn,
# 		$jobId, $prevD,   $cmdLine, $p,      $rBase
# 	);
# 
# 	my @pId             = flatten( $gOut[$oPI] );
# 	my $gnuSamStrand    = "XS:Z";
# 	my @gnuBsrIdxSuffix = ( "ct", "ga" );
# 
# 	$samBase = $pId[3];
# 	if ( $gRead[$rPER] == 1 ) {
# 		$samBase = $pId[2];
# 	}
# 
# 	foreach $rBase (@gReadList) {
# 		$samD    = "$gAliD/$rBase/$gBIN[$bGNU][$oOPRFX]/$gOut[$pMAPPED][1]";
# 		$samDfn  = "$samD/$samBase\.$gOut[$oMEXT][0]";
# 		$cmdLine = "\n";
# 
# 		for ( $s = $gSrchDim[$sSTRAND][0] ; $s < $gSrchDim[$sSTRAND][1] ; $s++ )
# 		{    #watson or crick
# 			$samStFn = "$samD/$samBase-$gnuBsrIdxSuffix[$s]\.$gOut[$oMEXT][0]";
# 			$cmdLine = cmd_grep( $cmdLine, $samDfn,
# 				"$gnuSamStrand:$gnuBsrIdxSuffix[$s]", $samStFn );
# 			( $prevD, $cmdLine, $bamFn ) =
# 			  cmd_sam2bam_full( $cmdLine, $gRefDfn, $samStFn, "keep" );
# 			( $prevD, $cmdLine, $gDummy ) =
# 			  cmd_sam2gmp( $cmdLine, $gRefDfn, $gnuBsrIdxSuffix[$s], $samStFn,
# 				"zip" );
# 		}
# 		$cmdLine = "$cmdLine\nrm -rf $samDfn\n";
# 		system($cmdLine);
# 	}
# }

#========================
sub cmd_sam2gmp {
	my ( $refDfn, $gnuConvNt, $inSam, $fopInAfter ) = @_;

	my $samBase = getFileName( $inSam, $gOut[$oMEXT][0] );

	#flipSamAligned2Crick
	my $flpCrTag     = "flpCr";
	my $flpCrSamBase = "$samBase-$flpCrTag";
	my $flpCrSam     = "$flpCrSamBase.$gOut[$oMEXT][0]"; # temp file
	
	my $cmdLine = "";
	my $gnuNtConvOpt;
	$cmdLine = cmd_flipCrickSam($cmdLine, $inSam, $flpCrSam);
	($gDummy, $gnuNtConvOpt) = select_gnumap_option($gRead[$rCONV], 0, 0); #since we flip crick sam, we only care about watson mapping and sense=all
	$cmdLine = cmd_sam2gmp_core($cmdLine, $refDfn, $flpCrSamBase, $gnuNtConvOpt, $flpCrSam);
	my $flpCrGmp = "$flpCrSamBase.$gBIN[$bGNUSGR][$bgsOPRFX]"; #this is output file and also temp file
	
	#filtering gmp
	my $gmpFn = "$samBase.$gBIN[$bGNUSGR][$bgsOPRFX]"; #to prepare final output gmp file
	$cmdLine = cmd_filterSam2GmpReport( $cmdLine, $gnuConvNt, $flpCrGmp, $gmpFn );

	#remove all tmp files
	$cmdLine = "$cmdLine\nrm -rf $flpCrSam $flpCrGmp\n";
	( $gDummy, $cmdLine, $gDummy ) = cmd_zip_util( $cmdLine, $inSam, $fopInAfter );
	system($cmdLine);
	return ( 0, $cmdLine, $gmpFn );
}

#============
sub cmd_sam2gmp_core {
	my ( $cmdLine, $refDfn, $flpCrSamBase, $gnuNtConvOpt, $flpCrSam ) = @_;

	$cmdLine = "$cmdLine\nexport LD_LIBRARY_PATH=$gBIN[$bGNU][$bgBIN_LIB]:\$LD_LIBRARY_PATH\n";
	$cmdLine = "$cmdLine\n$gBIN[$bGNUSGR][$bgsBIN] -g $refDfn -o $flpCrSamBase -c 1 $gnuNtConvOpt $flpCrSam\n";

	return ($cmdLine);
}

#====================
sub cmd_flipCrickSam {
	my ( $cmdLine, $inSam, $flpCrSam ) = @_;

	$cmdLine = "$cmdLine\n$gBIN[$bGNU][$bgiREVCOMP_SAM] $inSam $flpCrSam\n";
	return ($cmdLine);
}

#================================================================
sub cmd_filterSam2GmpReport {
	my ( $cmdLine, $gnuConvNt, $rawGmpFn, $gmpFn ) = @_;

	$cmdLine =
	  "$cmdLine\n$gBIN[$bGNUSGR][$bgsPOSTPROC] $rawGmpFn $gnuConvNt $gmpFn\n";

	return ($cmdLine);
}

# #===========================================================================
# #don't work for PER
# sub convertSam2fastq {
# 	my ( $cmdLine, $sam, $readDfn ) = @_;
# 
# 	$cmdLine =
# 	  "$cmdLine\nsleep 5\n$gBIN[$bSAM2FQ][$oBIN] -s $sam -f $readDfn\n";
# 	return ($cmdLine);
# }


# #===================
# sub qsubmitter1 {
# 	my ( $jobDsc, $cmdLine, $nodes, $ppn, $mem, $walltime ) = @_;
# 
# # Open the .qsub file for writing and then start dumping in commands until 'EOF'
# 	$gQsubCnt++;
# 
# 	chomp($cmdLine);
# 	my $jobId = "X";
# 	my ( $jname, $mem_option );
# 
# 	if ( $cmdLine eq '' ) { return $jobId; }
# 
# 	if ($gDebug) {
# 		print "qsub1:$jobDsc\n$cmdLine\n";
# 	}
# 
# 	my $qsubf = "$gQsubDir/${jobDsc}_${gQsubCnt}.qsub";
# 
# 	if ( $gHpcId eq "byu" ) {
# 		open( QSUBFILE, ">$qsubf" ) || die "cannot create $qsubf\n";
# 		print QSUBFILE "#! /bin/bash\n";
# 		print QSUBFILE "#PBS -N ${jobDsc}_${gQsubCnt}\n";
# 
# #print QSUBFILE "#PBS -l nodes=1:ppn=${ppn},mem=${mem},walltime=${walltime}:00:00\n";
# 		print QSUBFILE
# "#PBS -l nodes=1:ppn=${ppn}:beta,mem=${mem},walltime=${walltime}:00:00\n";
# 
# 		#print QSUBFILE "#PBS -k oe\n";
# 		#print QSUBFILE "#PBS -m e\n";
# 		#print QSUBFILE "#PBS -M cjhong\@bu.edu\n";
# 		print QSUBFILE "cd $gQsubDir\n\n";
# 
# 		#print QSUBFILE "echo \"$gQsubCnt:$cmdLine\"\n\n";
# 		print QSUBFILE "${cmdLine}\n\n";
# 		close(QSUBFILE);
# 
# 		print "qsub1:qsub ${qsubf}\n";
# 
# 		if ( $gOut[$pQSUB][0] == 1 ) {
# 			$jobId = `qsub ${qsubf}`;
# 		}
# 		else {
# 			system("$cmdLine");
# 		}
# 		sleep(5);
# 		chomp($jobId);
# 	}
# 	elsif ( $gHpcId eq "rubysky" ) {
# 
# 		#refer to Adam's email to complete this section
# 		open( QSUBFILE, ">$qsubf" ) || die "cannot create $qsubf\n";
# 		$jname = "${jobDsc}_${gQsubCnt}";
# 		print QSUBFILE "#!/bin/bash\n\n";
# 		print QSUBFILE "#\$ -cwd\n";
# 		print QSUBFILE "#\$ -N $jname\n";
# 		print QSUBFILE "#\$ -j y\n\n";
# 		print QSUBFILE "#\$ -m be\n";
# 		print QSUBFILE "#\$ -M cjhong\@bu.edu\n";
# 
# 		print QSUBFILE "#\$ -o $gQsubDir/$jname.log\n";
# 		print QSUBFILE "#\$ -pe single_node 1-${ppn}\n";
# 
# 		#print QSUBFILE "#\$ -l mem_free=${mem}\n";
# 		print QSUBFILE "#\$ -P johnsonlab\n";
# 
# #print QSUBFILE "#PBS -l nodes=1:ppn=${ppn}:bigmem,mem=${mem},walltime=${walltime}:00:00\n";
# 		print QSUBFILE
# "echo \"==========================================================\"\n";
# 		print QSUBFILE "echo \"Starting on       : \$(date)\"\n";
# 		print QSUBFILE "echo \"Running on node   : \$(hostname)\"\n";
# 
# 		#print QSUBFILE "echo \"Current job ID    : \$(JOB_ID)\"\n";
# 		#print QSUBFILE "echo \"Current job name  : \$(JOB_NAME)\"\n";
# 		#print QSUBFILE "echo \"Task index number : \$(TASK_ID)\"\n";
# 		print QSUBFILE
# "echo \"==========================================================\"\n";
# 
# 		print QSUBFILE "module load openmpi-x86_64\n";
# 		print QSUBFILE "cd $gQsubDir\n";
# 		print QSUBFILE "echo \"$gQsubCnt\"\n";
# 		print QSUBFILE "${cmdLine}\n";
# 
# 		print QSUBFILE
# "echo \"==========================================================\"\n";
# 		print QSUBFILE "echo \"Finished on       : \$(date)\"\n";
# 		print QSUBFILE
# "echo \"==========================================================\"\n";
# 		close(QSUBFILE);
# 		$mem_option = "-l mem_free=${mem}";
# 		print "$mem_option\n";
# 		system("qsub ${mem_option} ${qsubf}");
# 		sleep(5);
# 		$jobId = "1";
# 	}
# 	elsif ( $gHpcId eq "oligomer" ) {
# 
# 		open( QSUBFILE, ">$qsubf" ) || die "cannot create $qsubf\n";
# 		print QSUBFILE "#! /bin/bash\n";
# 		print QSUBFILE "#PBS -N ${jobDsc}_${gQsubCnt}\n";
# 		print QSUBFILE
# 		  "#PBS -l nodes=1:ppn=${ppn},mem=${mem},walltime=${walltime}:00:00\n";
# 
# 		print QSUBFILE "cd $gQsubDir\n\n";
# 
# 		#print QSUBFILE "echo \"$gQsubCnt:$cmdLine\"\n\n";
# 		print QSUBFILE "${cmdLine}\n\n";
# 		close(QSUBFILE);
# 		system("$cmdLine");
# 		$jobId = $gQsubCnt;
# 
# 		#sleep(2);
# 	}
# 
# 	if ($gDebug) {
# 		print "qsub1[$jobDsc]----MYJOBID[$jobId]\n";
# 	}
# 	return $jobId;
# }
