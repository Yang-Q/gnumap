#include <fstream>
#include <iostream>
#include <sstream>
#include <string>
#include <cstring>
#include "const_include.h"

unsigned int gBISULFITE = 0;
unsigned int gBISULFITE2 = 0;
unsigned int gATOG = 0;
unsigned int gATOG2 = 0;
unsigned int gRNA = 0;
unsigned int gSNP = 0;
float gSNP_PVAL = 0.001;
bool gSNP_MONOP = false;
bool gPRINT_VCF = false;
unsigned char gINT2BASE[16];
unsigned char g_bs_CONVERSION[256];
unsigned char g_gen_CONVERSION[256];
std::string gNtConvFlag = "nn";

int gVERBOSE=1;
#define DEF_MER_SIZE	10
int gMER_SIZE=0;
unsigned int gMAX_MATCHES = 20;
unsigned int gMAX_PRE_MATCHES_HALF = PRE_MATCH_BACKUP_FOLDS_STRAND * gMAX_MATCHES; //when a total number of hits in hash is greater than this cutoff, we believe that this is one of repeated region or non-informative redundant mapping so we skip mapping on this read to speed up.

//unsigned int gMAX_PRE_MATCHES = PRE_MATCH_BACKUP_FOLDS * gMAX_MATCHES; //this variable becomes effective when -t is used. If there are too many top hits in prematch process, we think the number of top matches in refined matching process will be greater than the number of max list, which we will not consider this mapping. As of 09/24/2012, I disabled this feature!

bool gUNIQUE = false;
const unsigned int gREAD_BUFFER = 1024 * 5;
unsigned int gBUFFER_SIZE = 1024*512;
unsigned int gMEM_BUFFER_SIZE = 67108864; //1024*1024*64;
//unsigned int gMEM_BUFFER_SIZE = 536870912; //1024*1024*512;
//unsigned int gMEM_BUFFER_SIZE = 1048576; //1024*1024;
//unsigned int gMEM_BUFFER_SIZE = 1024*1024*64; //1024*1024;
//unsigned int gBUFFER_SIZE = 1024;
unsigned int gHASH_SIZE = 1048576;	//gHASH_SIZE is set as 4^10
unsigned int gMAX_HASH_SIZE = 1000;
//unsigned int gMAX_HASH_SIZE = 0;
float gALIGN_SCORE = 0.9;
bool perc = true;
float gCUTOFF_SCORE = 0.0;

float gADJUST = 0.25;
float gMATCH = 3;
float gTRANSITION=-2;
float gTRANSVERSION=-3;
float gTRANSVERSION_bs=-4;
float gMATCH_offset=1;
float gGAP = -4;

float gGAP_SW = 2*gGAP;
float gTRANSITION_SW = gTRANSITION;
float gTRANSVERSION_SW = gTRANSVERSION;

// These are the defaults
float gPHMM_match=0.98;
float gPHMM_syn=0.01;
float gPHMM_nsyn=0.005;
//float gPHMM_match=0.97;
//float gPHMM_syn=0.01;
//float gPHMM_nsyn=0.01;

int gMAX_GAP = 3;
//int gMAX_GAP_SW = 2*gMAX_GAP;
//unsigned int gSwMinLen = 18;
float gNW_Eval_cutoff = 1.0;
float gSW_Eval_cutoff = 1e-30;
float gEvalLambda = 0.4173;

char* g_adaptor = NULL;
int gPRINT_FULL=0;
float gALIGN_SCORES[256][4];
float gPHMM_ALIGN_SCORES[256][4];
int gSEQ_LENGTH = 101;
unsigned int gMinReadL = 22;

bool gFAST = false;	// Do we perform a fast alignment or a more thorough one?
bool gMATCH_POS_STRAND = true;
bool gMATCH_NEG_STRAND = true;

//the following four lines will be not used and will be deleted soon!
// bool gMATCH_POS_STRAND_RE = true;
// bool gMATCH_NEG_STRAND_RE = false;
// bool gMATCH_POS_STRAND_BS = true;
// bool gMATCH_NEG_STRAND_BS = false;

bool gSAM2GMP = false; // added by CJ on 04/23/2012, an option for printing gmp file in an output
bool gSOFT_CLIP = false; // added by CJ to support smith-waterman local alignment for remote homology detection
bool gSKIP_TOO_MANY_REP = true;

//trimming adaptor is not fully implemented so please do preprocess including trimming low quality or adaptor
float gMIN_ADAPTOR_DIFF = 0.75;	//used to define the minimum difference from the given g_adaptor sequence to the consensus
unsigned int gMIN_CHOPPED_BASES = 4;	//if there's less than 4 bases that match the g_adaptor, we're not going to chop them.
int gGEN_SKIP = 0;	// The number of bases we skip as we hash the genome
unsigned int gGEN_SIZE = 8;	// The number of bases each GEN segment contains (per each unsigned int)
unsigned int gJUMP_SIZE = 1;	// The number of bases we jump in a read (per hash)
//unsigned int gMIN_JUMP_MATCHES = 2; // CJ won't use this param any more!

int gSGREX = false;

// This is for making MPI work
bool gMPI_LARGEMEM = false;

//bool gILLUMINA = false;
int gPHRED_OFFSET = SANGER_PHRED_OFFSET;

int gTOP_K_HASH = 1;

//I'm not sure what BIN_SIZE means...It's not used anymore--only appearance is in Genome.cpp:365
//int BIN_SIZE = 10;
int gBIN_SIZE = 1;
 
// gSAVE and gREAD are global variables used to indicate whether the hashed genome should be
// read and saved or created from fasta files.
bool gSAVE=false;
bool gREAD=false;
char* gSAVE_FN;
char* gREAD_FN;

int gPAIR_ID = 0;
std::string gPAIR_ID_STR = "0";

unsigned long gGenSizeBp;
unsigned long gGenSizeBp0 = 0;

int bad_files =0;

void printReadPWM(Read* read);

int iproc=0,nproc=1;

void InitProg() {
	memset(g_gen_CONVERSION,4,256);
    g_gen_CONVERSION[(int)'a'] = g_gen_CONVERSION[(int)'A'] = 0;
    g_gen_CONVERSION[(int)'c'] = g_gen_CONVERSION[(int)'C'] = 1;
    g_gen_CONVERSION[(int)'g'] = g_gen_CONVERSION[(int)'G'] = 2;
    g_gen_CONVERSION[(int)'t'] = g_gen_CONVERSION[(int)'T'] = 3;
    g_gen_CONVERSION[(int)'n'] = g_gen_CONVERSION[(int)'N'] = 4;
    g_gen_CONVERSION[(int)'\n'] = g_gen_CONVERSION[(int)'\r'] = 5;
    // Strange white-space characters
    g_gen_CONVERSION[10] = g_gen_CONVERSION[11] = g_gen_CONVERSION[12] = g_gen_CONVERSION[13] = 5;
    g_gen_CONVERSION[(int)'\0'] = 6;
    g_gen_CONVERSION[(int)'>'] = 7;

    memset(g_bs_CONVERSION,4,256);
    g_bs_CONVERSION[(unsigned int)'a'] = g_bs_CONVERSION[(unsigned int)'A'] = 0;
    g_bs_CONVERSION[(unsigned int)'c'] = g_bs_CONVERSION[(unsigned int)'C'] = 1;
    g_bs_CONVERSION[(unsigned int)'g'] = g_bs_CONVERSION[(unsigned int)'G'] = 2;
    g_bs_CONVERSION[(unsigned int)'t'] = g_bs_CONVERSION[(unsigned int)'T'] = 3;
    g_bs_CONVERSION[(unsigned int)'n'] = g_bs_CONVERSION[(unsigned int)'N'] = 4;
    g_bs_CONVERSION[(unsigned int)'\n'] = g_bs_CONVERSION[(unsigned int)'\r'] = 5;
    // Strange white-space characters
    g_bs_CONVERSION[10] = g_bs_CONVERSION[11] = g_bs_CONVERSION[12] = g_bs_CONVERSION[13] = 5;
    g_bs_CONVERSION[(unsigned int)'\0'] = 6;
    g_bs_CONVERSION[(unsigned int)'>'] = 7;

    memset(gINT2BASE,'?',16);
    gINT2BASE[0] = 'a';
    gINT2BASE[1] = 'c';
    gINT2BASE[2] = 'g';
    gINT2BASE[3] = 't';
    gINT2BASE[4] = 'n';	
}
