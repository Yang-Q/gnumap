#include <iostream>
#include <math.h>

#include "gsl/gsl_cdf.h"
#include "GenomeMem.h"

#include "const_include.h"
#include "const_define.h"

using namespace std;

int arg_max(float a, float b, float c) {
    if(a >= b) {
        if(a >= c) // a >= c, a >= b
            return 1;
        else  // c > a >= b
            return 3;
    }
    else { // b > a
        if(b >= c) // b > a, b >= c
            return 2;
        else // c > b > a
            return 3;
    }
}

float max(float a, float b, float c) {
	int max_a = arg_max(a, b, c);
	if(max_a == 1)
		return a;
	if(max_a == 2)
		return b;
	return c;
}

bool is_indel(float ins_count, float del_count, float norm_count) {
	return false;
}

float indelLRT(float ins_count, float del_count, float norm_count, int &type) {
	type = arg_max(ins_count, del_count, norm_count);
	float max_val = 0;
	if(type == 1)
		max_val = ins_count;
	else if(type == 2)
		max_val = del_count;
	else
		// It's definitely not an indel if norm is the highest
		return 1;
		//max_val = norm_count;

	float sum = ins_count+del_count+norm_count;

	double ratio = (1.0/pow(3.0f, sum)) /
					( pow((sum-max_val)/sum/2, sum-max_val) );

	double pval = 1-gsl_cdf_chisq_P(-2*log(ratio), 1);

	return pval;
}

// Make these global variables
string chromosome_name;
Genome gen;

#define UNIFORM_PRIOR	0.2
#define BUFF_SIZE	1024
float MAX_PVAL = 0.1;

int prev_ins_counts_n[BUFF_SIZE];
int   ins_counts_n[BUFF_SIZE];
float ins_counts  [BUFF_SIZE];
float del_counts  [BUFF_SIZE];
float norm_counts [BUFF_SIZE];


// Process the array up to, but not including, nex_pos
void process_array_piece(char* chr_name, float max_pval, unsigned int cur_pos, unsigned int next_pos) {
#ifdef DEBUG
	for(unsigned int i=0; i<BUFF_SIZE; i++) {
		fprintf(stderr, "%d\t", cur_pos+i);
	}
	fprintf(stderr, "\n");
	for(unsigned int i=0; i<BUFF_SIZE; i++) {
		fprintf(stderr, "%02.3f\t", norm_counts[ (cur_pos+i) % BUFF_SIZE ]);
	}
	fprintf(stderr, "\n");
	for(unsigned int i=0; i<BUFF_SIZE; i++) {
		fprintf(stderr, "%02.3f\t", del_counts[ (cur_pos+i) % BUFF_SIZE ]);
	}
	fprintf(stderr, "\n");
	for(unsigned int i=0; i<BUFF_SIZE; i++) {
		fprintf(stderr, "%02.3f\t", ins_counts[ (cur_pos+i) % BUFF_SIZE ]);
	}
	fprintf(stderr, "\n");
	for(unsigned int i=0; i<BUFF_SIZE; i++) {
		fprintf(stderr, "%05d\t", ins_counts_n[ (cur_pos+i) % BUFF_SIZE ]);
	}
	fprintf(stderr, "\n");
	for(unsigned int i=0; i<BUFF_SIZE; i++) {
		fprintf(stderr, "%05d\t", prev_ins_counts_n[ (cur_pos+i) % BUFF_SIZE ]);
	}
	fprintf(stderr, "\n\n");
#endif

	unsigned int arr_pos;
	int indel_type;
	float pval;

	for(unsigned int i=cur_pos; i<next_pos; i++) {
		arr_pos = i % BUFF_SIZE;
		ins_counts[arr_pos] += UNIFORM_PRIOR;
		del_counts[arr_pos] += UNIFORM_PRIOR;
		norm_counts[arr_pos] += UNIFORM_PRIOR;
		pval = indelLRT(ins_counts[arr_pos], del_counts[arr_pos], norm_counts[arr_pos], 
					indel_type);
		if(pval <= max_pval) {
			// Print out a record
			if(indel_type == 1) {
				/*
				fprintf(stdout, "%s\t%u\t.\t%c\t%c?\t.\t.\tINS-IntoReads(%d);pval=%f\n",
						chr_name, i, gen.GetChar(i-1), gen.GetChar(i-1), 
						ins_counts_n[arr_pos], pval);
				*/
				fprintf(stdout, "%s\t%u\t.\t%c\t%c", chr_name, i, gen.GetChar(i-1), gen.GetChar(i-1));
				for(int j=0; j<ins_counts_n[arr_pos]; j++) {
					fprintf(stdout, "?");
				}
				fprintf(stdout, "\t.\t.\tINS-IntoReads(%d);pval=%f\n", ins_counts_n[arr_pos], pval);
			}
			else {
				fprintf(stdout, "%s\t%u\t.\t%c%c\t%c\t.\t.\tDEL-FromReads;pval=%f\n",
						chr_name, i, gen.GetChar(i-1), gen.GetChar(i), 
						gen.GetChar(i), pval);
			}
		}

		prev_ins_counts_n[arr_pos] = 0;
		ins_counts_n[arr_pos] = 0;
		ins_counts  [arr_pos] = 0;
		del_counts  [arr_pos] = 0;
		norm_counts [arr_pos] = 0;
	}
}

#define MAX(a,b) a > b ? a : b

void process_single_sam(unsigned int cur_pos, unsigned int record_pos, char* record_cigar, float record_posterior) {
	unsigned int arr_pos = record_pos % BUFF_SIZE;

	char* temp = record_cigar;
	int count;
	char type;
	bool first = true;

	while(strlen(temp) > 1) {
		sscanf(temp, "%d%c", &count, &type);
		if(count >= 1000)
			temp += 4;
		else if(count >= 100)
			temp += 4;
		else if(count >= 10)
			temp += 3;
		else
			temp += 2;

		if(type == 'M') {
			for(int i=0; i<count; i++) {
				norm_counts[ (arr_pos+i) % BUFF_SIZE ] += record_posterior;
			}
			arr_pos += count;
		}
		else if(type == 'I') {
			if(!first) {
				ins_counts[ arr_pos ] += count*record_posterior;
				prev_ins_counts_n[ arr_pos ] = ins_counts_n[ arr_pos ];
				ins_counts_n[ arr_pos ] =
							MAX(ins_counts_n[ arr_pos ], count);	// Set it as the max of this or
																	// any other record
			}
		}
		else if(type == 'D') {
			if(!first) {
				for(int i=0; i<count; i++) {
					del_counts[ (arr_pos+i) % BUFF_SIZE ] += record_posterior;
				}
				arr_pos += count;
			}
		}
		first = false;
	}

	if(type == 'I') {	// need to roll back the changes if the last one was an insertion/deletion
		ins_counts[ arr_pos ] -= record_posterior;
		ins_counts_n[ arr_pos ] = prev_ins_counts_n[ arr_pos ];
	}
	else if(type == 'D') {
		for(int i=0; i<count; i++) {
			del_counts[ (arr_pos-i) % BUFF_SIZE ] -= record_posterior;
		}
	}
}

int main(int argc, char* argv[]) {
	if(argc < 6) {
		fprintf(stderr, "usage: sam_indel <INFILE.sam> <SEQ_LENGTH> <PVAL> <GENOME> <CHR_NAME>\n");
		return -1;
	}
	
	// Global vars necessary for GenomeMem
	InitProg();

	char* gINPUT = argv[1];
	unsigned int gSEQ_LEN = atoi(argv[2]);
	float gPVAL = atof(argv[3]);
	char* gGENFILE = argv[4];
	char* gCHR_NAME = argv[5];

	gen.use(gGENFILE);
	gen.StoreGenome(false);	// don't need to make extra arrays--just want genome stored

	unsigned int num_processed = 0;
	unsigned int cur_pos = 0;
	string line;
	bool first = true;
	ifstream infile(gINPUT);
	if(infile.is_open()) {
		char record_cigar[1024];
		char chr_name[16];
		unsigned int record_pos;
		float record_posterior;
		while(infile.good()) {
			getline(infile, line);
			if(line.size() == 0)
				continue;

			int ret = sscanf(line.c_str(), "%*s %*d %s %u %*d %s = 0 0 %*s %*s XA:f:%*f XP:f:%f", 
					chr_name, &record_pos, record_cigar, &record_posterior);
			if(ret != 4) {
				fprintf(stderr, "ERROR (couldn't find two elements) on line %s\n", line.c_str());
				// Don't break here... :/
				continue;
			}
			if(strcmp(chr_name, gCHR_NAME) != 0)
				continue;
			
			if(first) {
				cur_pos = record_pos;
				first = false;
			}
			if(record_pos+gSEQ_LEN >= cur_pos+BUFF_SIZE) {
				// process the array up to record_pos
				process_array_piece(gCHR_NAME, gPVAL, cur_pos, record_pos);
				cur_pos = record_pos;
			}
			process_single_sam(cur_pos, record_pos, record_cigar, record_posterior);
			num_processed++;
		}
		process_array_piece(gCHR_NAME, gPVAL, cur_pos, record_pos+gSEQ_LEN+1);
	}
	else {
		// Didn't open
		fprintf(stderr, "Couldn't open gINPUT file: %s\n", gINPUT);
		perror("");
		return -1;
	}

	fprintf(stderr, "Finished.  %u sequences processed\n", num_processed);
}
