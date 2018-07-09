#ifndef SEQUENCE_OPERATIONS_H
#define SEQUENCE_OPERATIONS_H

#include <string>
#include <algorithm>
#include <cstring>
#include "const_include.h"

using namespace std;

inline string reverse_comp(string &s) {
	string temp = "";
	int len = s.length();
	
	for(int i=1; i<=len; i++) {
		switch(s[len-i]) {
			case 'a':
				temp += 't';
				break;
			case 'A':
				temp += 'T';
				break;
			case 't':
				temp += 'a';
				break;
			case 'T':
				temp += 'A';
				break;
			case 'c':
				temp += 'g';
				break;
			case 'C':
				temp += 'G';
				break;
			case 'g':
				temp += 'c';
				break;
			case 'G':
				temp += 'C';
				break;
			case '-':
				temp += '-';
				break;
			default:
				temp += 'n';
				break;
		}
	}
	
	return temp;
}

inline bool str2Lower(const string &s, const string &r) {
	
	bool needProb = false;
	
	//debug-todo
	//cout << "read:" << s << endl;
	//cout << "ref:" << r << endl;
	
	for (unsigned int j=0; j<s.length(); j++) {
		if (g_bs_CONVERSION[(unsigned int)s[j]] != g_bs_CONVERSION[(unsigned int)r[j]]){
			needProb = true;
			break;
		}
	}
	return needProb;
}

inline string reverse_qual(string &s, int len) {
	string temp = "";
	//int len = s.length();
	
	for(int i=1; i<=len; i++) {
		temp += s[len-i];
	}
	
	return temp;
}

inline string reverse_CIGAR(const char* s) {
	string temp_CIGAR;
	string number = "";
	for(unsigned int i=0; i<strlen(s); i++) {
		if(48 <= (int)s[i] && (int)s[i] <= 58) {
			number += s[i];
		}
		else {	// found a character.  Reset number and prepend to CIGAR
			temp_CIGAR = number+s[i] + temp_CIGAR;
			number = "";
		}
	}

	return temp_CIGAR;
}

inline void reverse_comp(float** pwm, unsigned int length) {
	float temp_pwm[length][4];


	//for(r_vit = v.rbegin(); r_vit != v.rend(); ++r_vit) {
	for(unsigned int i = 0; i < length; ++i) {
		//need to switch places with the characters: a-t, g-c
		//also need to rev-com it.
		temp_pwm[length-1-i][A_POS] = pwm[i][T_POS];
		temp_pwm[length-1-i][C_POS] = pwm[i][G_POS];
		temp_pwm[length-1-i][G_POS] = pwm[i][C_POS];
		temp_pwm[length-1-i][T_POS] = pwm[i][A_POS];
	}
		
	for(unsigned int i=0; i<length; i++) {
		for(unsigned int j=0; j<4; j++) {
			pwm[i][j] = temp_pwm[i][j];
		}
	}
	
	//cout << "temp's size: " << temp.size() << endl;
	//cout << "originals size: " << v.size() << endl; 
}

inline float** reverse_comp_cpy(float** pwm, unsigned int length) {
	float** temp_pwm = new float*[length];

	for(unsigned int i=0; i<length; i++) {
		temp_pwm[length-1-i] = new float[4];
		temp_pwm[length-1-i][A_POS] = pwm[i][T_POS];
		temp_pwm[length-1-i][C_POS] = pwm[i][G_POS];
		temp_pwm[length-1-i][G_POS] = pwm[i][C_POS];
		temp_pwm[length-1-i][T_POS] = pwm[i][A_POS];
	}

	return temp_pwm;
}

inline float** create_pwm(int readL) {
	int i;
	float** temp_pwm = new float*[readL+1];
	for (i=0; i<=readL; i++) {
		temp_pwm[i] =  new float[4];
		temp_pwm[i][A_POS] = 0.0;
		temp_pwm[i][C_POS] = 0.0;
		temp_pwm[i][G_POS] = 0.0;
		temp_pwm[i][T_POS] = 0.0;
	}
	return temp_pwm;
}

inline float** read_create_copy(float **pwm, int first, int last) {
	int i, i0;
	float** temp_pwm = new float*[last-first+1];
	for (i=first; i<=last; i++) {
		i0 = i - first;
		temp_pwm[i0] =  new float[4];
		temp_pwm[i0][A_POS] = pwm[i][A_POS];
		temp_pwm[i0][C_POS] = pwm[i][C_POS];
		temp_pwm[i0][G_POS] = pwm[i][G_POS];
		temp_pwm[i0][T_POS] = pwm[i][T_POS];
	}
	return temp_pwm;
}

inline void read_copy(float **pwm_dest, float **pwm, int first, int last) {
	//todo: check the dimension of pwm_dest
	int i, i0;
	for (i=first; i<=last; i++) {
		i0 = i - first;
		pwm_dest[i0][A_POS] = pwm[i][A_POS];
		pwm_dest[i0][C_POS] = pwm[i][C_POS];
		pwm_dest[i0][G_POS] = pwm[i][G_POS];
		pwm_dest[i0][T_POS] = pwm[i][T_POS];
	}
}

/*
infile float** read_ext_ends_cpy(float **pwm, unsigned int length) {
	float** temp_pwm = new float*[length+2*gMAX_GAP2];
	
	for(unsigned int i=0; i<gMAX_GAP2; i++) {
		temp_pwm[i] = new float[4];
		temp_pwm[length+i] = new float[4];
	}
	
	for(unsigned int i=0; i<length; i++) {
		temp_pwm[i+gMAX_GAP2] = new float[4];
		temp_pwm[i+gMAX_GAP2][A_POS] = pwm[i][T_POS];
		temp_pwm[i+gMAX_GAP2][C_POS] = pwm[i][G_POS];
		temp_pwm[i+gMAX_GAP2][G_POS] = pwm[i][C_POS];
		temp_pwm[i+gMAX_GAP2][T_POS] = pwm[i][A_POS];
	}
	return temp_pwm;
}
*/

inline float** reverse_comp_cpy_phmm(float** pwm, unsigned int length) {
	float** temp_pwm = new float*[length];

	for(unsigned int i=0; i<length; i++) {
		temp_pwm[length-1-i] = new float[NUM_SNP_VALS];
		temp_pwm[length-1-i][A_POS] = pwm[i][T_POS];
		temp_pwm[length-1-i][C_POS] = pwm[i][G_POS];
		temp_pwm[length-1-i][G_POS] = pwm[i][C_POS];
		temp_pwm[length-1-i][T_POS] = pwm[i][A_POS];
		temp_pwm[length-1-i][N_POS] = pwm[i][N_POS];
#ifdef _INDEL
		temp_pwm[length-1-i][INS_POS] = pwm[i][INS_POS];
		temp_pwm[length-1-i][DEL_POS] = pwm[i][DEL_POS];
#endif
	}

	return temp_pwm;
}

inline unsigned int max_flt(float array[], unsigned int array_size) {
	float* pos = max_element(array, array+array_size);
	unsigned int max_pos = pos - &array[0];
	return max_pos;
}

// Need the MAX_PRB for fasta files where there's a 100% call at each base
// Otherwise, it won't print out any characters
#define MAX_PRB	0.9999

inline string str2qual(Read &r) {
	if(r.fq.size() > 0)
		return r.fq;

	string qual = "";

	unsigned int i;
	for(i=0; i<r.length; i++) {
		int max_pos = max_flt(r.pwm[i],4);
		if(r.pwm[i][max_pos] > MAX_PRB)
			//if(gILLUMINA) 	qual += (char)((-10 * log((1-MAX_PRB)/MAX_PRB)/log(10.))+33);
			//else 			
			qual += (char)((-10 * log(1-MAX_PRB)/log(10.))+33);
			
		else {
			//fprintf(stderr,"max_pos of %d and resulting character of %c\n",
			//		max_pos, (char)((-10 * log(1-r.pwm[i][max_pos])/log(10))+33));
			//if(gILLUMINA) 	qual += (char)((-10 * log((1-r.pwm[i][max_pos])/r.pwm[i][max_pos])/log(10.))+33);
			//else			
			qual += (char)((-10 * log(1-r.pwm[i][max_pos])/log(10))+33);
		}
	}

	return qual;
}

inline void delete_read(Read* temp_read) {
	//done with the read.  Delete it.
	for(unsigned int j=0; j<temp_read->length; j++) {
		delete[] temp_read->pwm[j];
	}
	delete[] temp_read->pwm;
	if(temp_read->name)
		delete[] temp_read->name;
	
	delete temp_read;	
}

#endif
