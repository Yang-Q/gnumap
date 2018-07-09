#ifndef BIN_SEQ_H
#define BIN_SEQ_H

#include "UnitTest.h"
#include "const_include.h"
#include "Reader.h"
#include <cassert>
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <sstream>
#include <algorithm>

using namespace std;

#define DEF_ARR_SIZE	101*108
#define NEG_INF			-100000

enum CIGAR_TYPE {
	MATCH, INS, DEL, NOT_ALIGNED
};

#define CIGAR_TO_STR(k) (k == MATCH) ? 'M' : (k == INS) ? 'I' : (k == DEL) ? 'D' : 'S'

class bin_seq {
	private:
		// Pair HMM Values
		float PHMM_q;
		float PHMM_t;
		float PHMM_d;
		float PHMM_e;
		float PHMM_Tmm;
		float PHMM_Tgm;
		float PHMM_Tmg;
		float PHMM_Tgg;


	public:
		bin_seq() {
			PHMM_q = 0.25;
			PHMM_t = 0.05;
			PHMM_d = 0.0025; // usual GAP score
			//PHMM_d = 0.005; // testing GAP score
			PHMM_e = 0.5;
			PHMM_Tmm = 1-2*PHMM_d-PHMM_t;
			PHMM_Tgm = 1-PHMM_d-PHMM_t;
			PHMM_Tmg = PHMM_d;
			PHMM_Tgg = PHMM_e;

			MXY_size = DEF_ARR_SIZE;
			fMXY = new double[MXY_size];
			bMXY = new double[MXY_size];
			pMXY = new double[MXY_size];
			def_arr = new float[DEF_ARR_SIZE];
			
			//def_arr_lite = new float[DEF_ARR_SIZE];
			for(unsigned int i=0; i<MXY_size; i++) {
				fMXY[i] = 0;
				bMXY[i] = 0;
				pMXY[i] = 0;
				def_arr[i] = 0;
				//def_arr_lite[i] = 0;
			}
		
			def_arr_size = DEF_ARR_SIZE;
			def_moves_arr = 0;
			def_moves_size = 0;
			
			//def_arr_size_lite = DEF_ARR_SIZE;
			//def_moves_arr_lite = 0;
			//def_moves_size_lite = 0;
			
		};

		~bin_seq() {
			delete[] def_arr;
			//delete[] def_arr_lite;
			delete[] fMXY;
			delete[] bMXY;
			delete[] pMXY;
			def_arr_size = 0;
			
			if(def_moves_arr)
				delete[] def_moves_arr;
			def_moves_size = 0;
			
			//def_arr_size_lite = 0;
			//if(def_moves_arr_lite)
			//	delete[] def_moves_arr_lite;
			//def_moves_size_lite = 0;
			
		}
		
		static unsigned int str2int();
		static unsigned int str2int(string &a);
		static string hash2str(const unsigned int h, int length);
		
		/*
		 * get_hash will get the hashing value for either a given string or
		 * a given location on genome.  The return value of get_hash is a pair
		 * of bool and BIT64_t, representing whether the hash was successful, and,
		 * if it was, the value for the hash.
		 * get_hash will return false only if there is an unknown character found in
		 * the middle of the sequence.  For example, actgNgcta will not be a unique 
		 * 9-mer, so it is not added to the hash.  In addition, when trying to map
		 * back to a genome, if the read has an unknown character in the given K-mer,
		 * it should not map back exactly to any portion of the genome.
		 */
		static pair<bool,unsigned long> get_hash(string str);
		static pair<bool,string> get_hash_str(string str);
		static pair<bool,unsigned long> get_hash(const char* str);
		static pair<int,unsigned long> get_hash(unsigned int &pos, const unsigned char* genome);
		static pair<int,unsigned long> get_hash(unsigned int &pos, const unsigned char* genome, const unsigned int max_pos);
		static pair<int,string> get_hash_str(unsigned int &pos, const unsigned char* genome, const unsigned int max_pos);
		
		inline unsigned int base2int(char c);
		inline char int2base(unsigned int i);
	
		/*
		 * Make sure we have enough space on the heap allocated for these sequences
		 */
		void check_pointer_length(unsigned int, unsigned int);

		/*
		 * Alignment without the optimizations
		 */
		float get_align_score(const Read&, const string&, int align_mode, int iHalf, float earlyTerminationScore);
		
		/*
		 * Optimize the alignment by guaranteeing there is a matching portion
		 */
		float get_align_score(const Read&, const string&, unsigned int, unsigned int);
		float get_align_score_begin(const Read &, const string &, unsigned int, int, float);
		float smith_waterman_score_begin(const Read &, const string &, const int);
		float get_align_score_mid(const Read &, const string &, unsigned int, unsigned int);
		float get_align_score_end(const Read &, const string &, unsigned int);
		float get_align_score(Read &read, string &gen);

		/*
		 * pairHMM will perform the pair HMM on the two strings, the one being the read and
		 * the other the genome.  The result will come in a double array of floats representing
		 * the amount at each position.  This can then be mapped to the genome.
		 */
		float** pairHMM(const Read&, const string&, const string&);
		inline float p(char, char);
		inline float pam_p(int, char);
		inline float p_seq(float x[4], char);
		
		//string get_align_score_w_traceback(const Read &read, const string &consense, const string &gen);
		/**
		 * The paired version of the get_align_score_w_traceback will align the two sequences and
		 * return strings that represent the alignment.  The first string will be the result of the
		 * alignment in the first sequence.  The second sequence will be the SAM-style CIGAR 
		 * alignment values, where only the following characters are used:
		 * 
		 * M => Alignment Match (can be match or mismatch)
		 * I => Insertion into Reference
		 * D => Deletion from the Reference
		 *
		 * Unlike the SAM output format where the resulting sequence does not have any gaps, the
		 * first sequence will be the gapped sequence so the .sgrex format can still be used.
		 *
		 * @param read - The pwm representing the read to be aligned
		 * @param consense - The consensus sequence of the read
		 * @param gen - The genomic sequence to which this sequence will be mapped
		 * 
		 * @return pair<string,string> where pair.first is the alignment sequence (with dashes for
		 *		gaps) and pair.second is the CIGAR of the alignment
		 */
		pair<string,string> get_align_score_w_traceback(const Read &read, const string &consense, const string &gen);
		pair<string,string> get_align_score_w_traceback2(const Read &read, const string &consense, const string &gen, double &align_score); // added by CJ
		pair<string,string> smith_waterman_w_track(const Read &read, const string &consense, const string &gen, double &align_score, int &gen_start);
			
		int get_align_score(string seq1, string seq2);	//archaic
		void detect_nw_ends(const Read &read, const string &gen, int &gen_start, int &gen_last, int sensing_end_mode); // added by CJ
		
		static bool Test_arr_eq(float[][5], int, int, float**, int, int);
		static bool Test(ostream &os, unsigned int&);
		
		
	private:	
		float get_val(char read, char genome);
		float get_val(const float* a, const char b);

		int max(int a, int b, int c);
		float max_flt(float a, float b, float c);
		
		float max_flt(float a, float b, float c, float z); // this is for SW
		unsigned int max_flt(float array[], unsigned int array_size);
		float max_flt(char&, float, float, float);
		float max_flt(char&, float, float, float, float); //this is for SW
		void print_align(string str1, string align, string str2);
		
		inline char max_char(const float* chr) const; //temp to be removed!
		
		float* def_arr;
		unsigned int def_arr_size;
		char* def_moves_arr;
		unsigned int def_moves_size;
		
		//the following template is used for detect_nw_ends
		//float* def_arr_lite;
		//unsigned int def_arr_size_lite;
		//char* def_moves_arr_lite;
		//unsigned int def_moves_size_lite;
		

		double* fMXY;
		double* bMXY;
		double* pMXY;
		unsigned int MXY_size;
};

#endif
