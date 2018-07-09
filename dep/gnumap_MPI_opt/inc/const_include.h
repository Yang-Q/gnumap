#ifndef CONST_INCLUDE_H
#define CONST_INCLUDE_H

#define gVERSION "3.1.0 BETA"

#include <cmath>
#include <string>
#include <sys/time.h>
#include <ctype.h>

typedef unsigned long long BIT64_t;

//If it's wanted to have a base-pair resolution, un-comment these lines and comment the next two.
// NOTE:  CANNOT CHANGE THIS TYPE WITHOUT RUINING MPI
typedef unsigned int GEN_TYPE;
#define GEN_PACK	8

extern unsigned int gGEN_SIZE;

//#define MAX_MISMATCH 100

// This flag will create a lot of verbose output.
//#define DEBUG

#define MAX_LN_SZ	1024
#define MAX_NAME_SZ	1024
#define MAX_CIGAR_SZ	1024
#define MAX_MER_SIZE	32

#define NUM_CHAR_VALS	5
#ifdef _INDEL
#define NUM_SNP_VALS	NUM_CHAR_VALS+2
#else
#define NUM_SNP_VALS	NUM_CHAR_VALS
#endif
#define A_POS 0
#define C_POS 1
#define G_POS 2
#define T_POS 3
#define N_POS 4
#define DEL_POS 5 /* read gap, genome insertion */
#define INS_POS 6 /* genome gap, read insertion */

#define READS_PER_PROC	1024*2

#define MUTEX_LOCK(t) pthread_mutex_lock(t)
#define MUTEX_UNLOCK(t) pthread_mutex_unlock(t)

#if defined( INTDISC )
#define FLTFROM255(t, c) ((float)(c)/255.0f*(t))
#define CHARFROMFLT(t, f) (unsigned char)(round((f)/(t)*255))

inline void adjust255(float prevV, float newV, unsigned char *curvals) {
	int i;
	if(prevV != 0 && prevV!=newV && newV>1)
		i=0;
	for(i=0; i<NUM_SNP_VALS; i++) {
		curvals[i] = CHARFROMFLT(newV, FLTFROM255(prevV, curvals[i]));
	}
}
#endif

#define ILLUMINA_PHRED_OFFSET 64
#define SANGER_PHRED_OFFSET 33
#define NO_FASTQ 0
#define EVAL_K 0.711

/******** Needed for Probcons *********/
//int VERSION = 1;
const int NumInsertStates = 2;
/**************************************/

extern unsigned char g_bs_CONVERSION[256];
extern unsigned char g_gen_CONVERSION[256];

extern unsigned char gINT2BASE[16];
extern int gVERBOSE;
extern int gMER_SIZE;
extern unsigned int gBUFFER_SIZE;
extern const unsigned int gREAD_BUFFER;
extern unsigned int gHASH_SIZE;	
extern unsigned int gMAX_HASH_SIZE;
extern bool gUNIQUE;
extern unsigned int gMAX_MATCHES;
extern unsigned int gMAX_PRE_MATCHES_HALF;

extern bool gMATCH_POS_STRAND;
extern bool gMATCH_NEG_STRAND;
// extern bool gMATCH_POS_STRAND_RE;
// extern bool gMATCH_NEG_STRAND_RE;
// extern bool gMATCH_POS_STRAND_BS;
// extern bool gMATCH_NEG_STRAND_BS;

extern float gALIGN_SCORE;
extern float gALIGN_SCORES[256][4];
extern float gPHMM_ALIGN_SCORES[256][4];

extern int gMAX_GAP;
//extern int gMAX_GAP_SW;
extern int gSEQ_LENGTH;
extern float gGAP;
extern float gADJUST;
extern float gMATCH;
extern float gTRANSITION;
extern float gTRANSVERSION;
extern float gTRANSVERSION_bs;
extern float gMATCH_offset;
extern float gGAP_SW; // added by CJ
extern float gTRANSITION_SW;
extern float gTRANSVERSION_SW;
extern unsigned int gSwMinLen;
extern unsigned long gGenSizeBp;
extern unsigned long gGenSizeBp0;

extern float gNW_Eval_cutoff;
extern float gSW_Eval_cutoff;
extern float gEvalLambda;

extern int gGEN_SKIP;
extern unsigned int gBISULFITE;
extern unsigned int gBISULFITE2;
extern std::string gNtConvFlag;
extern unsigned int gATOG;
extern unsigned int gATOG2;
extern unsigned int gRNA;
extern unsigned int gSNP;
extern float gSNP_PVAL;
extern bool gSNP_MONOP;
extern bool gPRINT_VCF;
extern unsigned int gMinReadL;

extern int gBIN_SIZE;

extern char* g_adaptor;
extern float gMIN_ADAPTOR_DIFF;
extern unsigned int gMIN_CHOPPED_BASES;

//extern bool gILLUMINA;
extern int gPHRED_OFFSET;
extern int gTOP_K_HASH;
extern bool gSAVE;
extern bool gREAD;
extern char* gSAVE_FN;
extern char* gREAD_FN;

extern int iproc, nproc;
extern bool gMPI_LARGEMEM;
extern bool gSAM2GMP; //added by CJ
extern bool gSOFT_CLIP; //added by CJ
extern bool gSKIP_TOO_MANY_REP; //added by CJ
extern int gPAIR_ID;
extern std::string gPAIR_ID_STR;

enum FileType {
	CHROM=0,
	SEQUENCE=1
};

enum SeqType {
	INT,
	PRB,
	FASTA,
	FASTQ
};

struct Read {
	char* name;
	Read() { name = 0; length=0; }
	
	Read(float** p, int l) {
		pwm = p; length = l;
	}

	unsigned int length;
	float** pwm;
	std::string seq;
	std::string fq;
};

struct GenomeLocation {
	GenomeLocation() {
		amount = 0.0;
		packed_char = 0;
	}

	GenomeLocation(GEN_TYPE pc, float amt)
		: amount(amt), packed_char(pc) {}

	float amount;
	GEN_TYPE packed_char;
};

#define READ_FAILED		-1
#define READ_TOO_SHORT	-2
#define READ_TOO_POOR	-3
#define READ_TOO_MANY	999999
#define SAME_DIFF		0.00001
#define SUBOPT_TOP_R 0.01

#define POS_STRAND	0
#define NEG_STRAND	1

#define POS_UNIQUE 0
#define POS_MULTIPLE 1
#define NEG_UNIQUE 2
#define NEG_MULTIPLE 3

#define TOO_MANY_MAP -2
#define MAP_NEED2OPT 1

#define SEN_GEN_START 0
#define SEN_GEN_LAST 1

#define NW_ALIGN 0
#define SW_ALIGN 1
#define SPLICE_JT 2

#define PRE_TOP_R 0.25

#define EARLY_STOP_SC_R 0.44
#define EARLY_STOP_SC_R_AFTER_HIT 0.49
#define PRE_MATCH_BACKUP_FOLDS_STRAND 500
#define SW_ALIGN_R_CUTOFF 0.50

#define BASE3_PVAL_OFFSET 6.0
#define BASE4_PVAL_OFFSET 3.0

struct TopReadOutput {
	char READ_NAME[MAX_NAME_SZ];
	char CHR_NAME[MAX_NAME_SZ];
	unsigned long CHR_POS;
	int strand;
	int MAPQ;
	float A_SCORE;
	float POST_PROB;
	int SIM_MATCHES; //it represents an instance id mapped by a given query
	float DENOM;
	char CIGAR[MAX_CIGAR_SZ];
	unsigned int readIndex;
	std::string consensus;
	std::string qual;
	unsigned int genStrLen;
};

/**********************************************
 * These next are for Command-Line Parsing
 **********************************************/
#define INVALID_ARG -1
#define NO_MATCHING_ARG -2
#define UNKNOWN_ARG -3
#define PARSE_ERROR -10

struct CMDLine {
	int errnum;
	int errpos;
};
/**********************************************/

// Needed mostly for testing, but also other stuff
/* Return the current time in seconds, using a float precision number. */
inline double When() {
	struct timeval tp;
	gettimeofday(&tp, NULL);
	return ((double) tp.tv_sec + (double) tp.tv_usec * 1e-6);
}

#endif //CONST_INCLUDE_H
