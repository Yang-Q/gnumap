#include "bin_seq.h"
	
#define ABS(v) (v)>0 ? (v) : -(v)

inline float bin_seq::p(char x, char y) {
	//if(x == y) return 10f*4f/5f; else return 10*1/5;
	if(x == y) return 8; else return 2;
}

/*
 * Use the built-in PAM we have created for the alignment here.
 * We can multiply by 10 to keep the values in check--we just need to make sure
 * the scores are proportional.
 */
inline float bin_seq::pam_p(int x, char y) {
	//return 20*gPHMM_ALIGN_SCORES[(int)x][(int)y];
	return gPHMM_ALIGN_SCORES[(int)y][x];
}
	
inline float bin_seq::p_seq(float x[4], char y) {
	float sum = 0;
/*
	sum += x[A_POS] * p('a',tolower(y));
	sum += x[C_POS] * p('c',tolower(y));
	sum += x[G_POS] * p('g',tolower(y));
	sum += x[T_POS] * p('t',tolower(y));
//*/
///*
	sum += x[A_POS] * pam_p(A_POS,tolower(y));
	sum += x[C_POS] * pam_p(C_POS,tolower(y));
	sum += x[G_POS] * pam_p(G_POS,tolower(y));
	sum += x[T_POS] * pam_p(T_POS,tolower(y));
//*/	
	// Just a scaling parameter to prevent over- or underflow
	return 3*sum;
}
							
							
float** bin_seq::pairHMM(const Read &r, const string &consensus, const string &genome) {
	//fprintf(stderr,"\naligning read %s: \n+%s\n=%s\n",r.name,consensus.c_str(),genome.c_str());
	
	int i,j; //i is for read index and j is gen_string_index. i is column index and j is for row index
	int n = r.length;
	int m = genome.size();
	int col_size = (m+1)*3;
	int p_col_size = m*3;

	unsigned int max_size = (n+1)*(m+1)*3;
	if(max_size > MXY_size) {
		// need to realloc them
		MXY_size = max_size;

		delete[] fMXY;
		delete[] bMXY;
		delete[] pMXY;
		
		fMXY = new double[MXY_size];
		bMXY = new double[MXY_size];
		pMXY = new double[MXY_size];
	}

	// Forward
	// Create Matrices
	//double fMXY[n+1][(m+1)*3];
	// backward matrices
	//double bMXY[n+1][(m+1)*3];
	/*
	double ** fMXY = new double*[n+1];
	double ** bMXY = new double*[n+1];
	for(i=0; i<n+1; i++) {
		fMXY[i] = new double[(m+1)*3];
		bMXY[i] = new double[(m+1)*3];
	}*/

	// set them to zero
	/*
	for(i=0; i<(int)max_size; i++) {
		fMXY[i] = 0;
		bMXY[i] = 0;
		pMXY[i] = 0;
	}
	*/
	memset(fMXY,0,max_size*sizeof(double));
	memset(bMXY,0,max_size*sizeof(double));
	memset(pMXY,0,max_size*sizeof(double));
	
	// final matrices
	float** pGenScore = new float*[m];
	//double pMXY[n][3*m];
	
	/*
	for(i=0; i<n+1; i++) {
		for(j=0; j<m+1; j++) {
			fMXY[i][3*j] = 0;
			fMXY[i][3*j+1] = 0;
			fMXY[i][3*j+2] = 0;

			bMXY[i][3*j] = 0;
			bMXY[i][3*j+1] = 0;
			bMXY[i][3*j+2] = 0;

		}
	}
	*/
	for(j=0; j<m; j++) {
		//pGenScore[i] = new float[5];
		pGenScore[j] = new float[NUM_SNP_VALS];
		pGenScore[j][A_POS] = 0;
		pGenScore[j][C_POS] = 0;
		pGenScore[j][G_POS] = 0;
		pGenScore[j][T_POS] = 0;
		pGenScore[j][N_POS] = 0;
#ifdef _INDEL
		pGenScore[j][INS_POS] = 0;
		pGenScore[j][DEL_POS] = 0;
#endif // _INDEL
		
	}

	// Forward
	fMXY[0] = 1;
	for(i=1; i<n+1; i++) {
		for(j=1; j<m+1; j++) {
		
			//fM[i][j] = p_seq(r.pwm[i-1],genome[j-1]) * (PHMM_Tmm*fM[i-1][j-1]
			//	+ PHMM_Tgm*fX[i-1][j-1]
			//	+ PHMM_Tgm*fY[i-1][j-1]);
			//fX[i][j] = PHMM_q*(PHMM_Tmg*fM[i-1][j]+PHMM_Tgg*fX[i-1][j]);
			//fY[i][j] = PHMM_q*(PHMM_Tmg*fM[i][j-1]+PHMM_Tgg*fY[i][j-1]);
			//not quite sure if this recursion is correct!
			fMXY[i*col_size + 3*j] = p_seq(r.pwm[i-1],genome[j-1]) * (PHMM_Tmm*fMXY[(i-1)*col_size + 3*(j-1)]
				+ PHMM_Tgm*fMXY[(i-1)*col_size + 3*(j-1) +1]
				+ PHMM_Tgm*fMXY[(i-1)*col_size + 3*(j-1) +2]);
			fMXY[i*col_size + 3*j +1] = PHMM_q*(PHMM_Tmg*fMXY[(i-1)*col_size + 3*j]+PHMM_Tgg*fMXY[(i-1)*col_size + 3*j +1]); //insertion in read
			fMXY[i*col_size + 3*j +2] = PHMM_q*(PHMM_Tmg*fMXY[i*col_size + 3*(j-1)]+PHMM_Tgg*fMXY[i*col_size + 3*(j-1) +2]); //deletion in read
		}
	}

	//double fE=PHMM_t*(fM[n][m]+fX[n][m]+fY[n][m]);
	double fE = PHMM_t*(fMXY[n*col_size + 3*m]+fMXY[n*col_size + 3*m +1]+fMXY[n*col_size + 3*m +2]);

	// Backward
	//bM[n-1][m-1] = bX[n-1][m-1] = bY[n-1][m-1] = PHMM_t;
	bMXY[(n-1)*col_size + 3*(m-1)] = bMXY[(n-1)*col_size + 3*(m-1) +1] = bMXY[(n-1)*col_size + 3*(m-1) +2] = PHMM_t;

	for(i=n-1; i>=0; i--) {
		for(j=m-1; j>=0; j--) {
			if(j == (m-1) && i == (n-1))
				continue;
			if(j == m-1) {
				//bM[i][j] = PHMM_q*PHMM_Tmg*bX[i+1][j];
				//bX[i][j] = PHMM_q*PHMM_Tgg*bX[i+1][j];
				//bY[i][j] = 0;
				bMXY[i*col_size + 3*j] = PHMM_q*PHMM_Tmg*bMXY[(i+1)*col_size + 3*j +1];
				bMXY[i*col_size + 3*j +1] = PHMM_q*PHMM_Tgg*bMXY[(i+1)*col_size + 3*j +1];
				bMXY[i*col_size + 3*j +2] = 0;
				continue;
			}
			if(i == n-1) {
				//bM[i][j] = PHMM_q*PHMM_Tmg*bY[i][j+1];
				//bY[i][j] = PHMM_q*PHMM_Tgg*bY[i][j+1];
				//bX[i][j] = 0;
				bMXY[i*col_size + 3*j] = PHMM_q*PHMM_Tmg*bMXY[i*col_size + 3*(j+1) +2];
				bMXY[i*col_size + 3*j +2] = PHMM_q*PHMM_Tgg*bMXY[i*col_size + 3*(j+1) +2];
				bMXY[i*col_size + 3*j +1] = 0;
				continue;
			}
			
			//bM[i][j]= p_seq(r.pwm[i+1],genome[j+1])*PHMM_Tmm*bM[i+1][j+1]
			//		+ PHMM_q*PHMM_Tmg*bX[i+1][j] + PHMM_q*PHMM_Tmg*bY[i][j+1];
			//bX[i][j]= p_seq(r.pwm[i+1],genome[j+1])*PHMM_Tgm*bM[i+1][j+1]
			//		+ PHMM_q*PHMM_Tgg*bX[i+1][j];
			//bY[i][j]= p_seq(r.pwm[i+1],genome[j+1])*PHMM_Tgm*bM[i+1][j+1]
			//		+ PHMM_q*PHMM_Tgg*bY[i][j+1];

			bMXY[i*col_size + 3*j]= p_seq(r.pwm[i+1],genome[j+1])*PHMM_Tmm*bMXY[(i+1)*col_size + 3*(j+1)]
					+ PHMM_q*PHMM_Tmg*bMXY[(i+1)*col_size + 3*j +1] + PHMM_q*PHMM_Tmg*bMXY[i*col_size + 3*(j+1) +2];
			bMXY[i*col_size + 3*j +1]= p_seq(r.pwm[i+1],genome[j+1])*PHMM_Tgm*bMXY[(i+1)*col_size + 3*(j+1)]
					+ PHMM_q*PHMM_Tgg*bMXY[(i+1)*col_size + 3*j +1];
			bMXY[i*col_size + 3*j +2]= p_seq(r.pwm[i+1],genome[j+1])*PHMM_Tgm*bMXY[(i+1)*col_size + 3*(j+1)]
					+ PHMM_q*PHMM_Tgg*bMXY[i*col_size + 3*(j+1) +2];

		}
	}

	for(i=1; i<n+1; i++) {
		for(j=1; j<m+1; j++) {
			//pM[i-1][j-1] = fM[i][j]*bM[i-1][j-1]/fE;
			//pX[i-1][j-1] = fX[i][j]*bX[i-1][j-1]/fE;
			//pY[i-1][j-1] = fY[i][j]*bY[i-1][j-1]/fE;
			
			pMXY[(i-1)*p_col_size + 3*(j-1)] = fMXY[i*col_size + 3*j]*bMXY[(i-1)*col_size + 3*(j-1)]/fE;
			pMXY[(i-1)*p_col_size + 3*(j-1) +1] = fMXY[i*col_size + 3*j +1]*bMXY[(i-1)*col_size + 3*(j-1) +1]/fE;
			pMXY[(i-1)*p_col_size + 3*(j-1) +2] = fMXY[i*col_size + 3*j +2]*bMXY[(i-1)*col_size + 3*(j-1) +2]/fE;
		}
	}


	for(j=0; j<m; j++) {
		for(i=0; i<n; i++) {
			//pGenScore[i][g_gen_CONVERSION[(unsigned)consensus[j]]] += pY[j][i] + pM[j][i];

#ifdef _INDEL
			///*
			pGenScore[j][g_gen_CONVERSION[(unsigned)consensus[i]]] 
					+= pMXY[i*p_col_size + j*3];
			pGenScore[j][INS_POS] += pMXY[i*p_col_size + j*3 + 2];
			pGenScore[j][DEL_POS] += pMXY[i*p_col_size + j*3 + 1];
			//*/
#else
			
//			pGenScore[j][g_gen_CONVERSION[(unsigned)consensus[i]]] 
//					+= pMXY[i*p_col_size + j*3 +2]+pMXY[i*p_col_size + j*3];
				pGenScore[j][g_gen_CONVERSION[(unsigned)consensus[i]]] 
					+= pMXY[i*p_col_size + j*3];
			
			/*
			pGenScore[i][INS_POS] += pMXY[j*p_col_size + i*3 + 2];
			pGenScore[i][DEL_POS] += pMXY[j*p_col_size + i*3 + 1];
			//*/			
#endif // _INDEL

		}
	}


	/*
	fprintf(stderr,"pMXY:\n");
	for(i=0; i<n+1; i++) {
		for(j=0; j<m+1; j++) {
			//printf("%.2f\t",pGenScore[i][j]);
			//fprintf(stderr,"%.2e,%.2e,%.2e\t",fMXY[i*col_size + 3*j],fMXY[i*col_size + 3*j +1],fMXY[i*col_size + 3*j +2]);
			fprintf(stderr,"%.2f,%.2f,%.2f\t",pMXY[i*p_col_size + 3*j],pMXY[i*p_col_size + 3*j +1],pMXY[i*p_col_size + 3*j +2]);
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	fprintf(stderr,"pGenScore:\n");
	for(i=0; i<n; i++) {
		for(j=0; j<NUM_CHAR_VALS; j++) {
			fprintf(stderr,"%.2f\t",pGenScore[i][j]);
		}
		fprintf(stderr,"\n");
	}
	fprintf(stderr,"\n");
	//*/
	
	return pGenScore;
}


string bin_seq::hash2str(const unsigned int h, int length) {
	char to_return[length];
	unsigned int temp = h;

	for(int i=1; i<=length; i++) {
		//to_return[length-i] = int2base((temp % 4));
		to_return[length-i] = gINT2BASE[(temp % 4)];
		temp >>= 2;
	}
	
	to_return[length] = '\0';
	
	return to_return;
}

pair<bool,unsigned long> bin_seq::get_hash(const char* str) {
	string str_str = str;
	
	return get_hash(str_str);
}

pair<bool,unsigned long> bin_seq::get_hash(string str) {
	unsigned long temp = 0;

	//ensure that the sequence doesn't overflow the possible size of an unsinged int,
	//as well as not overflowing the gMER_SIZE.
	for(unsigned int i=0; i<(sizeof(unsigned long) * 4) && i<str.length(); i++) {
		
		//int spot = base2int(str[i]); 
		int spot = g_bs_CONVERSION[(int)str[i]]; 
		if(spot > 3)	//if it's not a valid, hashable character
			return pair<bool,unsigned long>(false,0);
			
		temp <<= 2;
		temp |= spot;
	}
		
	//temp %= HASH_SIZE;
	return pair<bool,unsigned long>(true,temp);
}

pair<bool,string> bin_seq::get_hash_str(string str) {
	string toret = "";

	for(unsigned int i=0; i<str.length(); i++) {
		int spot = g_bs_CONVERSION[(int)str[i]];
		if(spot > 3)
			return pair<bool,string>(false,"");
		toret += str[i];
	}

	return pair<bool,string>(true,toret);
}


pair<int,unsigned long> bin_seq::get_hash(unsigned int &pos, const unsigned char* genome, const unsigned int max_pos) {
	unsigned long temp = 0;
	unsigned int spot = 0;

	while(genome[pos] == '\n')
		pos++;

	int found = 0;
	// Ensure that the sequence doesn't overflow the possible size of a unsigned int,
	// as well as not overflowing the MER_SIZE.  (We won't allow MER_SIZE to
	// overflow the size of an unsigned long--64 bits.  This means we can allow up to 32bp MER_SIZE)
	// We want to make sure we find exactly MER_SIZE bases (found), but the int i helps us determine
	// where on the genome we look for them--while ignoring newline characters.
	//for(int i=0,found=0; found<gMER_SIZE && found<64 && pos+i<BUFFER_SIZE; found++, i++) {
	int i;
	for(i=0; found<gMER_SIZE && pos+i<max_pos; i++) {
		//unsigned int spot = base2int(genome[pos+i]);
		spot = g_bs_CONVERSION[(unsigned int)genome[pos+i]];
		
		if(spot < 4) {
			temp <<= 2;
			temp |= spot;
		
			found++;
		}
		else if(spot == 4) {
			//cout << "Unknown character: <" << genome[pos+i] << ">" << endl;
			return pair<int,unsigned long>(0,0);
		}
		else if(spot == 5) {	//if there's an newline reached, this isn't an error: just ignore.
			if(i == 0) {	//if there's a newline reached on the first of a sequence...
				pos++;
				i=-1;
			}
		}
		else if(spot == 6) {
			//pos = BUFFER_SIZE;	// 5 comes from the flag for the end of the stream.
			return pair<int,unsigned long>(-1,0);
		}
		else if(spot == 7) {
			pos = pos+i;
			//return a flag for a new genome, with a number representing the number of bases skipped.
			return pair<int,unsigned long>(-2,found); 
		}
	}
	
	// This is because it couldn't read enough bases
	if(pos+i >= max_pos) {	
		//cout << "Couldn't read enough bases, pos:"<<pos<<" i:"<<i<<" max_pos:"<<max_pos << endl;
		return pair<bool,unsigned long>(0,0);
	}
	

	if(spot == 7) {
		return pair<int,unsigned long>(-2,found);
	}

	return pair<bool,unsigned long>(1,temp);
}

pair<int,string> bin_seq::get_hash_str(unsigned int &pos, const unsigned char* genome, const unsigned int max_pos) {
	string toret = "";
	unsigned int spot = 0;

	while(genome[pos] == '\n')
		pos++;

	int found = 0;
	// Ensure that the sequence doesn't overflow the possible size of a unsigned int,
	// as well as not overflowing the MER_SIZE.  (We won't allow MER_SIZE to
	// overflow the size of an unsigned long--64 bits.  This means we can allow up to 32bp MER_SIZE)
	// We want to make sure we find exactly MER_SIZE bases (found), but the int i helps us determine
	// where on the genome we look for them--while ignoring newline characters.
	//for(int i=0,found=0; found<gMER_SIZE && found<64 && pos+i<BUFFER_SIZE; found++, i++) {
	int i;
	for(i=0; found<gMER_SIZE && pos+i<max_pos; i++) {
		//unsigned int spot = base2int(genome[pos+i]);
		spot = g_bs_CONVERSION[(unsigned int)genome[pos+i]];
		
		if(spot < 4) {
			toret += genome[pos+i];
		
			found++;
		}
		else if(spot == 4) {
			//cout << "Unknown character: <" << genome[pos+i] << ">" << endl;
			return pair<int,string>(0,"");
		}
		else if(spot == 5) {	//if there's an newline reached, this isn't an error: just ignore.
			if(i == 0) {	//if there's a newline reached on the first of a sequence...
				pos++;
				i=-1;
			}
		}
		else if(spot == 6) {
			//pos = BUFFER_SIZE;	// 5 comes from the flag for the end of the stream.
			return pair<int,string>(-1,"");
		}
		else if(spot == 7) {
			pos = pos+i;
			//return a flag for a new genome, with a number representing the number of bases skipped.
			return pair<int,string>(-2,toret); 
		}
	}
	
	// This is because it couldn't read enough bases
	if(pos+i >= max_pos) {	
		//cout << "Couldn't read enough bases, pos:"<<pos<<" i:"<<i<<" max_pos:"<<max_pos << endl;
		return pair<bool,string>(0,"");
	}
	

	if(spot == 7) {
		return pair<int,string>(-2,toret);
	}

	return pair<bool,string>(1,toret);
}


/*
 * The next few functions are for producing a Needleman-Wunsch alignment
 * for two given sequences.
 */
void bin_seq::print_align(string str1, string align, string str2) {
	int s1_size = str1.size();
	for(int i=1; i<=s1_size; i++) 
		cout << str1[s1_size-i]; 
	cout << endl;
	
	int a_size = align.size(); 
	for(int i=1; i<=a_size; i++)
		cout << align[a_size-i];
	cout << endl;
	
	int s2_size = str2.size();
	for(int i=1; i<=s2_size; i++)
		cout << str2[s2_size-i];
	cout << endl;
}

/*
 * This is meant to align the first piece of the read.  Aligning from the (exact) matching
 * portion which was produced from the hash piece to the beginning of the read.
 *
 *  *
 *	  *		<- we want to match this first piece.
 *		*
 *        *
 *          X	<- end.  We know it must be an exact match
 *			  X
 *
 * 
 
float bin_seq::get_align_score_begin_rect(const Read &read, const string &gen, const unsigned int end, const unsigned int gen_end) {
	if(end == 0)	//trivial case
		return 0;

	float m_mm1, m_mm2, m_mm, gap1, gap2, max_score;
	int i,j;
	unsigned ui,uj;
	
	int SEQ_LENGTH = read.length;

	unsigned int seq1Length = SEQ_LENGTH;
	unsigned int genSegLen = gen_end;

	// Check to see if we need to extend the default array
	unsigned int squareLen = (seq1Length+1) * (genSegLen+1);
	if(squareLen > def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		//memset(def_arr,NEG_INF,squareLen*sizeof(float));
		for(ui=0; ui<squareLen; ui++) {
			//def_arr[i] = NEG_INF;
			def_arr[ui] = -9;
		}
		def_arr_size = squareLen;
	}
	

	if(squareLen > def_moves_size) {
		if(def_moves_arr)
			delete[] def_moves_arr;
		def_moves_arr = new char[squareLen];
		
		for(unsigned int i=0; i<squareLen; i++) {
			def_moves_arr[i] = '-';
		}
		def_moves_size = squareLen;
	}

	
	float* nm = def_arr;
	//char* moves = def_moves_arr;
	
	// The size of a single row
	unsigned int row_size = seq1Length+1;
	
	//Initialize the entire array to zeros.
	for(i=0; i<=SEQ_LENGTH; i++)
		//for(int j=0; j<=gen.size(); j++) {
		for(int j=i-gMAX_GAP-1; j<=i+gMAX_GAP+1; j++) {
			if(j < 0 || j > (int)gen.size() )
				continue;
			//nm[i*row_size + j] = NEG_INF;
			nm[i*row_size + j] = -9;
			//moves[i*row_size + j] = 'D';
		}
	
	//Initialize the first row to have the gap penalty.
	for(i=0; i<=gMAX_GAP+2; i++) {
		nm[i*row_size] = gGAP * i;
		//moves[i*row_size] = 'U';
	}
	//Initialize the first column to have the gap penalty.
	for(j=0; j<=gMAX_GAP+2; j++) {
		nm[j] = gGAP * j;
		//moves[j] = 'L';
	}
	//moves[0] = 'D';
	
//#ifdef DEBUG
			if(gVERBOSE > 2)
				for(unsigned int q=0; q<=end; q++) {
					for(unsigned int p=0; p<=gen_end; p++) {
						fprintf(stderr,"%.2f\t",nm[q*row_size + p]);
					}
					fprintf(stderr,"\n");
				}
				fprintf(stderr,"\n");
//#endif
	
	//Set up the NW matrix	
	i=0;
	while(i<SEQ_LENGTH) { //finish DP right after reaching to the bottom of row
		i++;
		//for(unsigned int j=1; j<=gen.size(); j++) {
		for(int j=i-gMAX_GAP; j<=i+gMAX_GAP; j++) {
			if(j <= 0)
				continue;
			if (j>(int)gen.size())
				break;
			
			m_mm1 = nm[(i-1)*row_size + (j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			gap1 = nm[(i-1)*row_size + j] + gGAP;
			gap2 = nm[i*row_size + (j-1)] + gGAP;
			
			max_score = max_flt(m_mm,gap1,gap2);
			nm[i*row_size + j] = max_score;

		}	
	}

//#ifdef DEBUG
			if(gVERBOSE > 2)
				for(unsigned int q=0; q<=end; q++) {
					for(unsigned int p=0; p<=gen_end; p++) {
						fprintf(stderr,"%.2f\t",nm[q*row_size + p]);
					}
					fprintf(stderr,"\n");
				}
				fprintf(stderr,"\n");
//#endif

	float align_score = nm[0]-j*gGAP;
	// Don't delete.  Keep in memory
	//delete[] nm;
	return align_score;

}
*/


inline char bin_seq::max_char(const float* chr) const {
	if((chr[0] == chr[1]) && (chr[0] == chr[2]) && (chr[0] == chr[3]))
		return 'n';
	if(chr[0] >= chr[1]) {
		if(chr[0] >= chr[2]) {
			if(chr[0] >= chr[3])
				return 'a';
			else
				return 't';
			}
		else {
			if(chr[2] >= chr[3])
				return 'g';
			else
				return 't';
		}
	}
	else {
		if(chr[1] >= chr[2]) {
			if(chr[1] >= chr[3])
				return 'c';
			else
				return 't';
		}
		else {
			if(chr[2] >= chr[3])
				return 'g';
			else return 't';
		}
	}

	return 'n';
}

void bin_seq::detect_nw_ends(const Read &read,const string &gen, int &gen_start, int &gen_last, int sensing_end_mode) {
	
	string read_consensus = "";

	for(unsigned int i=0; i<read.length; i++) {
		read_consensus += max_char(read.pwm[i]);
	}
/*  	
	cerr << "\naligning " << read_consensus << endl
	     <<   "with     " << gen << endl;

	cerr << "Scores: " << get_val(read.pwm[0],'a') << ' '
		 			   << get_val(read.pwm[0],'c') << ' '
		 			   << get_val(read.pwm[0],'g') << ' '		 			   
		 			   << get_val(read.pwm[0],'t') << endl; 
*/

// Note that this version is different from the one in get_align_score_begin_rect which is backward DP. But this is forward DP. This is nothing but for printing CIGAR.

	float m_mm1, m_mm2, m_mm, gap1, gap2, max_score, sw_max_score=0.0;
	
	int seqLen = read.length;
	int genSegLen = gen.size();
	
	int i,j,i_max,j_max,i_start = -1,j_start = -1;

	// Check to see if we need to extend the default array
	int squareLen = (seqLen+1) * (genSegLen+1);
	
	if(squareLen > (int)def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		//memset(def_arr,NEG_INF,squareLen*sizeof(float));
		for(i=0; i<squareLen; i++) {
			def_arr[i] = 0;
		}
		def_arr_size = squareLen;
	}
	
	if(squareLen > (int)def_moves_size) {
		if(def_moves_arr)
			delete[] def_moves_arr;
		def_moves_arr = new char[squareLen];
		
		for(i=0; i<squareLen; i++) {
			def_moves_arr[i] = '-';
		}
		def_moves_size = squareLen;
	}
	
	float* nm = def_arr;
	char* moves = def_moves_arr;
	
	// The size of a single row
	int row_size = genSegLen+1;
	
	//Initialize the entire array to zeros.
	for(i=0; i<=seqLen; i++)
		for(j=0; j<=genSegLen; j++) {
			if(j < 0 || j > genSegLen )
				continue;
			nm[i*row_size + j] = 0;
			moves[i*row_size + j] = 'D';
	}
	
	//Initialize the first row to have the gap penalty.
	for(i=0; i<=seqLen; i++)
		moves[i*row_size] = 'S'; //insertion w.r.t query_read
	
	//Initialize the second half column to have the gap penalty.
	for(j=0; j<=genSegLen; j++)
		moves[j] = 'S'; //deletion w.r.t query_read
	
	moves[0] = 'S'; //means start cell in s-w alignment

#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=seqLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
	
	// Print the "moves" matrix
	for(i=0; i<=seqLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%c\t",moves[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif

	//Set up the NW matrix	
	
	for (i=1; i<seqLen; i++) { //finish DP right after reaching to the bottom of row
		
		//for(unsigned int j=1; j<=gen.size(); j++) {
		for( j=1; j<=genSegLen; j++) {

			m_mm1 = nm[(i-1)*row_size + (j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			gap1 = nm[(i-1)*row_size + j] + gGAP_SW;
			gap2 = nm[i*row_size + (j-1)] + gGAP_SW;
			
			max_score = max_flt(moves[i*row_size + j],m_mm,gap1,gap2,0.0);
			nm[i*row_size + j] = max_score;

		}
	}

//this is the last iteration to obtain a cell whose score is max
	i = seqLen;
	i_max = i;
	sw_max_score = 0;
	for(j=1; j<=genSegLen; j++) {

		m_mm1 = nm[(i-1)*row_size + (j-1)];
		m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
		m_mm = m_mm1 + m_mm2;
		gap1 = nm[(i-1)*row_size + j] + gGAP_SW;
		gap2 = nm[i*row_size + (j-1)] + gGAP_SW;
			
		max_score = max_flt(moves[i*row_size + j],m_mm,gap1,gap2,0.0);
		nm[i*row_size + j] = max_score;
		
		if (max_score >= sw_max_score) {
			sw_max_score = max_score;
			j_max=j;
		}
	}
	
#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=seqLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
	
	// Print the "moves" matrix
	for(i=0; i<=seqLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%c\t",moves[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif

	//Trace back though the NW matrix to find the best alignment.

	i = i_max;
	j = j_max;

	//cerr << "imax:" << i << "\njmax:" << j << "\n";

	while( (i != 0) && (j != 0) ) {
		switch(moves[i*row_size + j]) {
			case 'D': // Made a diagonal transition
				i -= 1;
				j -= 1;
				break;
				
			case 'U': // Made a transition from the upper sequence
				i -= 1;
				break;
				
			case 'L': // Made a transition from the left
				j -= 1;
				break;
				
			case 'S': // stop here
				//i'm about to be out!
				i_start=i;
				j_start=j;
				i=0;
				j=0;
				break;
			
			default: 
				cerr << "ERROR at " << i << "," << j << "." << endl;
				cerr << "\tValue is: " << moves[i*row_size+j] << "(" << (int)moves[i*row_size+j] << ")" << endl;
		}
	}
	
	if ((i_start == -1) || (j_start == -1)) {
		i_start = i;
		j_start = j;
	}
	
	if (sensing_end_mode == SEN_GEN_START) {
		if (i_start>0) {
			j_start = j_start-i_start;
			if (j_start<0)
				j_start = 0;
		}
		gen_start = j_start;
		gen_last = j_max-1;
	}
	else {
		if (i_max<seqLen) {
			j_max = j_max+seqLen-i_max;
			if (j_max>genSegLen)
				j_max = genSegLen;
		}
		gen_start = j_start;
		gen_last = j_max-1;
	}

#ifdef DEBUG
	cerr << "gen_start:" << gen_start << "\n";
	cerr << "gen_last:" << gen_last << "\n";
#endif
}

pair<string,string> bin_seq::smith_waterman_w_track(const Read &read, const string &consense, const string &gen, double &align_score, int &gen_start) {
// Note that this version is different from the one in get_align_score_begin_rect which is backward DP. But this is forward DP. This is nothing but for printing CIGAR.

	float m_mm1, m_mm2, m_mm, gap1, gap2, max_score, sw_max_score=0.0;

	int readLen = read.length;
	int genSegLen = gen.size();

	int i,j,i_max,j_max;
	unsigned int ui;
	// Check to see if we need to extend the default array
	int squareLen = (readLen+1) * (genSegLen+1);

	int max_gap_sw = (int)floor(2*gMAX_GAP*(readLen+1)/gSEQ_LENGTH+0.5);
	if (max_gap_sw < 2*gMAX_GAP)
		max_gap_sw = 2*gMAX_GAP;
	
	if(squareLen > (int)def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		//memset(def_arr,NEG_INF,squareLen*sizeof(float));
		for(i=0; i<squareLen; i++) {
			def_arr[i] = 0;
		}
		def_arr_size = squareLen;
	}
	
	if(squareLen > (int)def_moves_size) {
		if(def_moves_arr)
			delete[] def_moves_arr;
		def_moves_arr = new char[squareLen];
		
		for(i=0; i<squareLen; i++) {
			def_moves_arr[i] = '-';
		}
		def_moves_size = squareLen;
	}
	
	float* nm = def_arr;
	char* moves = def_moves_arr;
	
	// The size of a single row
	int row_size = genSegLen+1;
	
	//Initialize the entire array to zeros.
	for(i=0; i<=readLen; i++) {
		for(j=i-max_gap_sw-1; j<=i+max_gap_sw+1 && j<=genSegLen; j++) {
			if(j < 0 )
				continue;
			nm[i*row_size + j] = 0;
			moves[i*row_size + j] = 'D';
		}
	}
	
	//Initialize the first column w/ S
	for(i=0; i<=max_gap_sw+2; i++)
		moves[i*row_size] = 'S'; //insertion w.r.t query_read
	
	//Initialize the first row w/ S
	for(j=0; j<=max_gap_sw+2; j++)
		moves[j] = 'S'; //deletion w.r.t query_read
	
	moves[0] = 'S'; //means start cell in s-w alignment

#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=readLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
	
	// Print the "moves" matrix
	for(i=0; i<=readLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%c\t",moves[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif
	
	//Set up the NW matrix	
	
	for (i=1; i<=readLen; i++) { //finish DP right after reaching to the bottom of row
		
		//for(unsigned int j=1; j<=gen.size(); j++) {
		for( j=i-max_gap_sw; j<=i+max_gap_sw; j++) {

			if (j<=0) 
				continue;
			if (j>genSegLen)
				break;
			
			m_mm1 = nm[(i-1)*row_size + (j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			gap1 = nm[(i-1)*row_size + j] + gGAP_SW;
			gap2 = nm[i*row_size + (j-1)] + gGAP_SW;
			
			max_score = max_flt(moves[i*row_size + j],m_mm,gap1,gap2,0.0);
			nm[i*row_size + j] = max_score;
			
			if (max_score>=sw_max_score) {
				sw_max_score = max_score;
				i_max=i;
				j_max=j;
			}
		}
	}
	
#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=readLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
	
	// Print the "moves" matrix
	for(i=0; i<=readLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%c\t",moves[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif
	
	//Trace back though the NW matrix to find the best alignment.
	string aligned = "";
	char CIGAR[1024];
	CIGAR[0]=0;

	int gaps = 0;
	
	i = i_max;
	
	//print 'S' on CIGAR
	CIGAR_TYPE c_type = NOT_ALIGNED;
	char temp[1024];
	
	if (i<readLen) {
		
		sprintf(CIGAR,"%d%c",(readLen-i),CIGAR_TO_STR(c_type));
		strcpy(temp,CIGAR);
		
		while (i<readLen) {
			aligned += "s";
			i++;
		}
	
	}
	
	i = i_max;
	j = j_max;
	int c_counter = 0;
	c_type = MATCH;
	
	bool stopTraceBack = false;
	
	while( (!stopTraceBack) && (i != 0) && (j != 0) ) {
		switch(moves[i*row_size + j]) {
			case 'D': // Made a diagonal transition
				aligned += consense[i-1];
				
				// If it's a mismatch or a match
				if(c_type == MATCH)
					c_counter++;
				else {
					// Only do this if we've already started 
					if(c_counter) {
						// prepend the CIGAR format
						strcpy(temp,CIGAR);
						//sprintf(CIGAR,"%c%d%s",'M',j,temp);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = MATCH;
					c_counter = 1;
				}

				i -= 1;
				j -= 1;

				break;
			case 'U': // Made a transition from the upper sequence
				//aligned += tolower(consense[i]);
				aligned += 'i';

				// prepend the CIGAR format
				if(c_type == INS)
					c_counter++;
				else {
					if(c_counter) {
						strcpy(temp,CIGAR);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = INS;
					c_counter = 1;
				}
				
				i -= 1;
				gaps++;

				break;
			case 'L': // Made a transition from the left
				//aligned_gen += gen[j];
				aligned += '-';

				if(c_type == DEL)
					c_counter++;
				else {
					if(c_counter) {
						// prepend the CIGAR format
						strcpy(temp,CIGAR);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = DEL;
					c_counter = 1;
				}
				j -= 1;
				gaps++;
				
				break;
			case 'S': // stop here
				//i'm about to be out!
				stopTraceBack = 1;
				break;
				
			default: 
				cerr << "ERROR at " << i << "," << j << "." << endl;
				cerr << "\tValue is: " << moves[i*row_size+j] << "(" << (int)moves[i*row_size+j] << ")" << endl;
				return pair<string,string>("","");
		}
	} // end of traceback
	
	if (c_counter > 0){
		strcpy(temp,CIGAR);
		sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
	}

	if (i>0) {
		gen_start = i;
		c_type = NOT_ALIGNED;
		strcpy(temp,CIGAR);
		sprintf(CIGAR,"%d%c%s",i,CIGAR_TO_STR(c_type),temp);
		
		while(i>0) {
			aligned += "s"; //pretend insertion so it will skip when bookkeeping
			i--;
		}
	}
	
	// Reverse the sequences
	string temp_aligned = aligned; 
	for(ui=0; ui<aligned.size(); ui++) {
		aligned[ui] = temp_aligned[(aligned.size()-1)-ui];
	}

	return pair<string,string>(aligned,CIGAR);
}


/*
 * this is a Smith-Waterman (Local) Alignment to align two remote homologies
 * only obtain an alignment score
 * always forward DP
 */


float bin_seq::smith_waterman_score_begin(const Read &read, const string &gen, const int readLen) {

// Note that this version is different from the one in get_align_score_begin_rect which is backward DP. But this is forward DP. This is nothing but for printing CIGAR.
// usually, readLen == gen.size(); square matrix

	float m_mm1, m_mm2, m_mm, gap1, gap2, max_score, sw_max_score=0.0;

	int genSegLen = gen.size();
	int max_gap_sw = (int)floor(2*gMAX_GAP*(readLen+1)/gSEQ_LENGTH+0.5);
	if (max_gap_sw < 2*gMAX_GAP)
		max_gap_sw = 2*gMAX_GAP;
	
	int i,j;
	//int i_max,j_max;

	// Check to see if we need to extend the default array
	int squareLen = (readLen+1) * (genSegLen+1);
	
	if(squareLen > (int)def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		//memset(def_arr,NEG_INF,squareLen*sizeof(float));
		for(i=0; i<squareLen; i++) {
			def_arr[i] = 0;
		}
		def_arr_size = squareLen;
	}
	
	float* nm = def_arr;
	
	// The size of a single row
	int row_size = genSegLen+1;
	
	//Initialize the entire array to zeros.
	for(i=0; i<=readLen; i++)
		for(j=i-max_gap_sw-1; j<=i+max_gap_sw+1 && j<=genSegLen; j++) {
			if (j < 0)
				continue;
			nm[i*row_size + j] = 0;
	}
	
#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=readLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif
	
	//Set up the SW matrix	
	for (i=1; i<=readLen; i++) { //finish DP right after reaching to the bottom of row
		for( j=i-max_gap_sw; j<=i+max_gap_sw; j++) {

			if (j<=0) 
				continue;
			if (j>genSegLen)
				break;
			
			m_mm1 = nm[(i-1)*row_size + (j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			gap1 = nm[(i-1)*row_size + j] + gGAP_SW;
			gap2 = nm[i*row_size + (j-1)] + gGAP_SW;
			
			max_score = max_flt(m_mm,gap1,gap2,0.0);
			nm[i*row_size + j] = max_score;
			
			if (max_score>=sw_max_score) {
				sw_max_score = max_score;
				//i_max=i;
				//j_max=j;
			}
		}
	}
	
#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=readLen; i++) {
		for(j=0; j<=genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif
	
	return sw_max_score;
}

//will return the alignment string
//Note: gen string is usually len(read)+2*gMAX_GAP-1 and to resolve the alignment problem at the both ends.
//positions and align_score will be updated. Finally, cigar will be generated too.
//this NW is a forward DP which a narrow rectange. j(row) x i(col) = gen x query_read

pair<string,string> bin_seq::get_align_score_w_traceback2(const Read &read, const string &consense, const string &gen, double &align_score) {
	/* 
	cerr << "\naligning " << gen << endl
	     <<   "with     " << consense << endl;
	cerr << "Scores: " << get_val(read.pwm[0],'a') << ' '
		 			   << get_val(read.pwm[0],'c') << ' '
		 			   << get_val(read.pwm[0],'g') << ' '		 			   
		 			   << get_val(read.pwm[0],'t') << endl; 
	//*/

// Note that this version is different from the one in get_align_score_begin_rect which is backward DP. But this is forward DP. This is nothing but for printing CIGAR.

	float m_mm1, m_mm2, m_mm, gap1, gap2, max_score;
	
	int seqLen = read.length;
	int genSegLen = gen.size();
	int i,j,lenDiff;
	int row_size = genSegLen+1;
	int col_size = seqLen+1;
	unsigned int ui;
	bool columnLead = false;

	int squareLen = row_size*col_size;

	int max_gap_nw = (int)floor(gMAX_GAP*(seqLen+1)/gSEQ_LENGTH+0.5);
	if (max_gap_nw < gMAX_GAP)
		max_gap_nw = gMAX_GAP;
	
	// Check to see if we need to extend the default array	// The size of a single row
	if(squareLen > (int)def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		//memset(def_arr,NEG_INF,squareLen*sizeof(float));
		for(i=0; i<squareLen; i++) {
			def_arr[i] = NEG_INF;
			//def_arr[ui] = -9;
		}
		def_arr_size = squareLen;
	}
	if(squareLen > (int)def_moves_size) {
		if(def_moves_arr)
			delete[] def_moves_arr;
		def_moves_arr = new char[squareLen];
		
		for(i=0; i<squareLen; i++) {
			def_moves_arr[i] = '-';
		}
		def_moves_size = squareLen;
	}
	
	float* nm = def_arr;
	char* moves = def_moves_arr;
	
	// determine if this is narrow or wide DP matrix
	if (genSegLen>=seqLen) {//narrow?
		lenDiff = genSegLen - seqLen;
		columnLead = true;
	}
	else {//wide? then rowLead
		lenDiff = seqLen - genSegLen;
	}
	
	if (columnLead) {
		//Initialize the entire array to zeros.
		for(i=0; i<=seqLen; i++)
			for(j=i-gMAX_GAP-1; j<=i+gMAX_GAP+1+lenDiff; j++) {
				if(j < 0)
					continue;
				if (j>genSegLen)
					break;
					
				nm[i*row_size + j] = NEG_INF;
				//nm[i*row_size + j] = -9;
				moves[i*row_size + j] = 'D';
			}
		
		//Initialize the first row to have the gap penalty.
		for(i=0; i<=gMAX_GAP+2; i++) {
			nm[i*row_size] = gGAP * i;
			moves[i*row_size] = 'U'; //insertion w.r.t query_read
		}
		
		//Initialize the second half column to have the gap penalty.
		for(j=0; j<=gMAX_GAP+2+lenDiff; j++) {
			nm[j] = gGAP * j;
			moves[j] = 'L'; //deletion w.r.t query_read
		}
		moves[0] = 'D';
	}
	else { // row based DP initialization
		//Initialize the entire array to zeros.
		for(j=0; j<=genSegLen; j++)
			for(i=j-gMAX_GAP-1; i<=j+gMAX_GAP+1+lenDiff; i++) {
				if(i < 0)
					continue;
				if (i>seqLen)
					break;
					
				nm[i*row_size + j] = NEG_INF;
				//nm[i*row_size + j] = -9;
				moves[i*row_size + j] = 'D';
			}
		
		//Initialize the first row to have the gap penalty.
		for(i=0; i<=gMAX_GAP+2; i++) {
			nm[i*row_size] = gGAP * i;
			moves[i*row_size] = 'U'; //insertion w.r.t query_read
		}
		
		//Initialize the second half column to have the gap penalty.
		for(j=0; j<=gMAX_GAP+2+lenDiff; j++) {
			nm[j] = gGAP * j;
			moves[j] = 'L'; //deletion w.r.t query_read
		}
		moves[0] = 'D';
	}
	
#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=SEQ_LENGTH; i++) {
		for(j=0; j<=(int)genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
	
	// Print the "moves" matrix
	for(i=0; i<=SEQ_LENGTH; i++) {
		for(j=0; j<=(int)genSegLen; j++) {
			fprintf(stderr,"%c\t",moves[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif

	if (columnLead) {
		//Set up the NW matrix
		for (i=1;i<=seqLen;i++) { //finish DP right after reaching to the bottom of row
			for( j=i-gMAX_GAP; j<=i+gMAX_GAP+lenDiff; j++) {
				if(j <= 0)
					continue;
				if (j>genSegLen)
					break;
				
				m_mm1 = nm[(i-1)*row_size + (j-1)];
				m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
				m_mm = m_mm1 + m_mm2;
				gap1 = nm[(i-1)*row_size + j] + gGAP;
				gap2 = nm[i*row_size + (j-1)] + gGAP;
				
				max_score = max_flt(moves[i*row_size + j],m_mm,gap1,gap2);
				nm[i*row_size + j] = max_score;
			}
		}
	}
	else {
		for (j=1;j<=genSegLen;j++) { //finish DP right after reaching to the bottom of row
			for( i=j-gMAX_GAP; i<=j+gMAX_GAP+lenDiff; i++) {
				if(i <= 0)
					continue;
				if (i>seqLen)
					break;
				
				m_mm1 = nm[(i-1)*row_size + (j-1)];
				m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
				m_mm = m_mm1 + m_mm2;
				gap1 = nm[(i-1)*row_size + j] + gGAP;
				gap2 = nm[i*row_size + (j-1)] + gGAP;
				
				max_score = max_flt(moves[i*row_size + j],m_mm,gap1,gap2);
				nm[i*row_size + j] = max_score;
			}
		}
	}
	
#ifdef DEBUG
	// Print the NW matrix
 	for(i=0; i<=seqLen; i++) {
		for(j=0; j<=(int)genSegLen; j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
	
	// Print the "moves" matrix
	for(i=0; i<=seqLen; i++) {
		for(j=0; j<=(int)genSegLen; j++) {
			fprintf(stderr,"%c\t",moves[i*row_size + j]);
		}
		cerr << endl;
	}
	cerr << endl;
#endif
	
	//now, modify max_score
	align_score = max_score;
	
	//Trace back though the NW matrix to find the best alignment.
	string aligned = "";
	//string aligned_gen = "";
	char CIGAR[1024];
	CIGAR[0] = 0;
	
	int gaps=0;

	i = seqLen;
	j = genSegLen;

	CIGAR_TYPE c_type = MATCH;
	int c_counter = 0;

	char temp[1024];
	while( (i != 0) && (j != 0) ) {
		switch(moves[i*row_size + j]) {
			case 'D': // Made a diagonal transition
				aligned += consense[i-1];

				// If it's a mismatch or a match
				if(c_type == MATCH)
					c_counter++;
				else {
					// Only do this if we've already started 
					if(c_counter) {
						// prepend the CIGAR format
						strcpy(temp,CIGAR);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = MATCH;
					c_counter = 1;
				}

				i -= 1;
				j -= 1;

				break;
			case 'U': // Made a transition from the upper sequence
				//aligned += tolower(consense[i]);
				aligned += 'i';
				
				// prepend the CIGAR format
				if(c_type == INS)
					c_counter++;
				else {
					if(c_counter) {
						strcpy(temp,CIGAR);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = INS;
					c_counter = 1;
				}
				
				i -= 1;
				gaps++;

				break;
			case 'L': // Made a transition from the left
				//aligned_gen += gen[j];
				aligned += '-';

				if(c_type == DEL)
					c_counter++;
				else {
					if(c_counter) {
						// prepend the CIGAR format
						strcpy(temp,CIGAR);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = DEL;
					c_counter = 1;
				}
				j -= 1;
				gaps++;
				break;
				
			default: 
				cerr << "ERROR at " << i << "," << j << "." << endl;
				cerr << "\tValue is: " << moves[i*row_size+j] << "(" << (int)moves[i*row_size+j] << ")" << endl;
				return pair<string,string>("","");
		}
	}
	// wrap up the current c_count
	
	// make sure both i and j are zero.  If not, add insertions and deletions
	// until they are both zero

	// At this point, either i or j must be zero.
	// decrementing i inserts an insertion
	// this will seldom happen!
	while(i > 0) { //hit the top wall
		//fprintf(stderr,"For sequence %s, i is %d\n",aligned.c_str(),i);
		//aligned_gen += "-";
		aligned += consense[i];
		if(c_type == INS)
			c_counter++;
		else { 
			strcpy(temp,CIGAR);
			sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);

			c_type = INS;
			c_counter = 1;
		}
		i--;
		gaps++;
	}
	
	// Decrementing j corresponds to inserting a deletion
	while(j > 0) {
		//fprintf(stderr,"For sequence %s, j is %d\n",aligned.c_str(),j);
		aligned += "-";
		//aligned_gen += gen[j];
		if(c_type == DEL)
			c_counter++;
		else {
			strcpy(temp,CIGAR);
			sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);

			c_type = DEL;
			c_counter = 1;
		}
		j--;
		gaps++;
	}

	// Only print out if there are unwritten items
	if(c_counter > 0) {
		char temp[1024];
		// Clean up the CIGAR from the front of the sequence
		strcpy(temp,CIGAR);
		sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
	}

	// Reverse the sequences
	string temp_aligned = aligned; 
	for(ui=0; ui<aligned.size();ui++) {
		aligned[ui] = temp_aligned[(aligned.size()-1)-ui];
	}
	
	//fprintf(stderr,"For sequences: \n%s\n%s\nCIGAR is %s\n",aligned.c_str(),aligned_gen.c_str(),CIGAR);
	//fprintf(stderr,"Alignment score is %f\n\n\n",nm[row_size*(SEQ_LENGTH+1)]);
	
	return pair<string,string>(aligned,CIGAR);
}

//will return the alignment string
pair<string,string> bin_seq::get_align_score_w_traceback(const Read &read, const string &consense, const string &gen) {
	/* 
	cerr << "\naligning " << gen << endl
	     <<   "with     " << consense << endl;
	cerr << "Scores: " << get_val(read.pwm[0],'a') << ' '
		 			   << get_val(read.pwm[0],'c') << ' '
		 			   << get_val(read.pwm[0],'g') << ' '		 			   
		 			   << get_val(read.pwm[0],'t') << endl; 
	//*/

	float m_mm1, m_mm2, m_mm, gap1, gap2, max_score;
	
	int SEQ_LENGTH = read.length;

	unsigned int seq1Length = SEQ_LENGTH;
	unsigned int genSegLen = gen.size();
	
	int max_gap_nw = (int)floor(gMAX_GAP*(seq1Length+1)/gSEQ_LENGTH+0.5);
	if (max_gap_nw < gMAX_GAP)
		max_gap_nw = gMAX_GAP;

	// Check to see if we need to extend the default array
	unsigned int squareLen = (seq1Length+1) * (genSegLen+1);
	if(squareLen > def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		//memset(def_arr,NEG_INF,squareLen*sizeof(float));
		for(unsigned int i=0; i<squareLen; i++) {
			def_arr[i] = NEG_INF;
		}
		def_arr_size = squareLen;
	}
	if(squareLen > def_moves_size) {
		if(def_moves_arr)
			delete[] def_moves_arr;
		def_moves_arr = new char[squareLen];
		
		for(unsigned int i=0; i<squareLen; i++) {
			def_moves_arr[i] = '-';
		}
		def_moves_size = squareLen;
	}
	
	float* nm = def_arr;
	char* moves = def_moves_arr;
	
	// The size of a single row
	unsigned int row_size = seq1Length+1;

	
	//Initialize the entire array to zeros.
	for(int i=0; i<=SEQ_LENGTH; i++)
		//for(int j=0; j<=gen.size(); j++) {
		for(int j=i-max_gap_nw-1; j<=i+max_gap_nw+1; j++) {
			if(j < 0 || j > (int)gen.size() )
				continue;
			nm[i*row_size + j] = NEG_INF;
			moves[i*row_size + j] = 'D';
		}
	
	//Initialize the first row to have the gap penalty.
	for(int i=0; i<=max_gap_nw+2; i++) {
		nm[i*row_size] = gGAP * i;
		moves[i*row_size] = 'U';
	}
	//Initialize the first column to have the gap penalty.
	for(int j=0; j<=max_gap_nw+2; j++) {
		nm[j] = gGAP * j;
		moves[j] = 'L';
	}
	moves[0] = 'D';

	//Set up the NW matrix	
	for(int i=1; i<=SEQ_LENGTH; i++) {
		//for(unsigned int j=1; j<=gen.size(); j++) {
		for(int j=i-max_gap_nw; j<=i+max_gap_nw; j++) {
			if(j <= 0)
				continue;
			if (j>(int)gen.size())
				break;
			
			m_mm1 = nm[(i-1)*row_size + (j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			gap1 = nm[(i-1)*row_size + j] + gGAP;
			gap2 = nm[i*row_size + (j-1)] + gGAP;
			
			max_score = max_flt(moves[i*row_size + j],m_mm,gap1,gap2);
			
			nm[i*row_size + j] = max_score;

		}	
	}

#ifdef DEBUG
	// Print the NW matrix
 	for(int i=0; i<=SEQ_LENGTH; i++) {
		for(unsigned int j=0; j<=gen.size(); j++) {
			fprintf(stderr,"%.2f\t",nm[i*row_size + j]);
			
		}
		cerr << endl;
	}
	cerr << endl;
	
	// Print the "moves" matrix
	for(int i=0; i<=SEQ_LENGTH; i++) {
		for(unsigned int j=0; j<=gen.size(); j++) {
			fprintf(stderr,"%c\t",moves[i*row_size + j]);
			
		}
		cerr << endl;
	}
	cerr << endl;
#endif

	//Trace back though the NW matrix to find the best alignment.
	string aligned = "";
	//string aligned_gen = "";
	char CIGAR[1024];
	CIGAR[0] = 0;
	
	int gaps=0;

	int i = SEQ_LENGTH;
	int j = gen.size();

	CIGAR_TYPE c_type = MATCH;
	int c_counter = 0;

	char temp[1024];		
	while( (i != 0) && (j != 0) ) {
		switch(moves[i*row_size + j]) {
			case 'D': // Made a diagonal transition
				aligned += consense[i-1];
				//aligned_gen += gen[j-1];

				// If it's a mismatch or a match
				//if(consense[i-1] != gen[j-1]) {
					if(c_type == MATCH)
						c_counter++;
					else {
						// Only do this if we've already started 
						if(c_counter) {
							// prepend the CIGAR format
							strcpy(temp,CIGAR);
							//sprintf(CIGAR,"%c%d%s",'M',j,temp);
							sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
						}
						c_type = MATCH;
						c_counter = 1;
					}
				//}

				i -= 1;
				j -= 1;

				break;
			case 'U': // Made a transition from the upper sequence

				//aligned_gen += "-";
				aligned += consense[i];

				// prepend the CIGAR format
				if(c_type == INS)
					c_counter++;
				else {
					if(c_counter) {
						strcpy(temp,CIGAR);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = INS;
					c_counter = 1;
				}
				
				i -= 1;
				gaps++;

				break;
			case 'L': // Made a transition from the left
				//aligned_gen += gen[j];
				aligned += '-';

				if(c_type == DEL)
					c_counter++;
				else {
					if(c_counter) {
						// prepend the CIGAR format
						strcpy(temp,CIGAR);
						sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
					}
					c_type = DEL;
					c_counter = 1;
				}

				j -= 1;
				gaps++;
				
				break;
			default: 
				cerr << "ERROR at " << i << "," << j << "." << endl;
				cerr << "\tValue is: " << moves[i*row_size+j] << "(" << (int)moves[i*row_size+j] << ")" << endl;
				return pair<string,string>("","");
		}
	}
	
	// make sure both i and j are zero.  If not, add insertions and deletions
	// until they are both zero

	// At this point, either i or j must be zero.
	// decrementing i inserts an insertion
	while(i > 0) {
		//fprintf(stderr,"For sequence %s, i is %d\n",aligned.c_str(),i);
		//aligned_gen += "-";
		aligned += consense[i];
		if(c_type == INS)
			c_counter++;
		else {
			strcpy(temp,CIGAR);
			sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
			c_type = INS;
			c_counter = 1;
		}

		i--;
		gaps++;
	}

	// Decrementing j corresponds to inserting a deletion
	while(j > 0) {
		//fprintf(stderr,"For sequence %s, j is %d\n",aligned.c_str(),j);
		aligned += "-";
		//aligned_gen += gen[j];
		if(c_type == DEL)
			c_counter++;
		else {
			strcpy(temp,CIGAR);
			sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
			c_type = DEL;
			c_counter = 1;
		}

		j--;
		gaps++;
	}
	
	// Only print out if there are unwritten items
	if(c_counter > 0) {
		char temp[1024];
		// Clean up the CIGAR from the front of the sequence
		strcpy(temp,CIGAR);
		sprintf(CIGAR,"%d%c%s",c_counter,CIGAR_TO_STR(c_type),temp);
	}

	// Reverse the sequences
	string temp_aligned = aligned; 
	for(unsigned int i=0; i<aligned.size(); i++) {
		aligned[i] = temp_aligned[(aligned.size()-1)-i];
	}

	/*
	// Don't need to reverse the genome because we don't use it
	string temp_gen = aligned_gen;
	for(unsigned int i=0; i<aligned_gen.size(); i++) {
		aligned_gen[i] = temp_gen[aligned_gen.size()-1-i];
	}
	 //*/
	
	//fprintf(stderr,"For sequences: \n%s\n%s\nCIGAR is %s\n",aligned.c_str(),aligned_gen.c_str(),CIGAR);
	//fprintf(stderr,"Alignment score is %f\n\n\n",nm[row_size*(SEQ_LENGTH+1)]);
	
	return pair<string,string>(aligned,CIGAR);
}

void bin_seq::check_pointer_length(unsigned int read_len, unsigned int gen_len) {

	// Check to see if we need a larger pointer length
	unsigned int squareLen = (read_len+1) * (gen_len+2*gMAX_GAP+1);
	if(squareLen > def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		//memset(def_arr,NEG_INF,squareLen*sizeof(float));
		for(unsigned int i=0; i<squareLen; i++) {
			def_arr[i] = NEG_INF;
		}

		def_arr_size = squareLen;
	}
}

/*
 * Puts all the begin, mid and end functions together and returns the value incorporated here(to get max align score).
 */
float bin_seq::get_align_score(const Read &read, const string &gen,
		const unsigned int begin, const unsigned int end) {
		
#ifdef DEBUG
	fprintf(stderr,"Aligning read named %s to %s\n",read->name,gen.c_str());
#endif

	assert(begin<=end);
	
	float value = 0;
	
	// Make sure we have enough space on the heap
	check_pointer_length(read.length, gen.size());
	
	value += get_align_score_begin(read,gen,begin,0,0.0);
	value += get_align_score_mid(read,gen,begin,end);
	value += get_align_score_end(read,gen,end);
	
	return value;
		
}

float bin_seq::get_align_score(const Read &read, const string &gen, int align_mode, int readLenMid, float earlyTerminationScore) {
	
	// Make sure we have enough space on the heap
	check_pointer_length(read.length, gen.size());
	//todo: do we have to use a different gap_cost depending on align_mode?
	if (align_mode==NW_ALIGN)
		return get_align_score_begin(read, gen, read.length, readLenMid, earlyTerminationScore);
	else
		return smith_waterman_score_begin(read, gen, read.length);

	return -1;
}


/*
 * This is meant to align the first piece of the read.  Aligning from the (exact) matching
 * portion which was produced from the hash piece to the beginning of the read.
 *
 *  *
 *	  *		<- we want to match this first piece.
 *		*
 *        *
 *          X	<- end.  We know it must be an exact match
 *			  X
 *
 */
float bin_seq::get_align_score_begin(const Read &read, const string &gen, const unsigned int end, int iHalf, float earlyTerminationScore) {
	if(end == 0)	//trivial case
		return 0;

	int dpLen=end; //dpLen=read.len=gen.len
	int i,j;
	int dp_size=dpLen+1;
	
	int max_gap_nw = (int)floor(gMAX_GAP*dp_size/gSEQ_LENGTH+0.5);
	if (max_gap_nw < gMAX_GAP)
		max_gap_nw = gMAX_GAP;
	
	int squareLen=dp_size*dp_size;
	float m_mm1,m_mm2,m_mm,gap1,gap2,max_score;
	float align_score=0.0;
	float mid_max_score=0.0;
	
	// Check to see if we need to extend the default array	// The size of a single row
	if(squareLen > (int)def_arr_size) {
		delete[] def_arr;
		def_arr = new float[squareLen];
		
		for(i=0; i<squareLen; i++) {
			def_arr[i] = NEG_INF;
		}
		def_arr_size = squareLen;
	}
	
	float* nm = def_arr;
	for(i=0; i<=dpLen; i++) {
		for(j=i-max_gap_nw-1; j<=i+max_gap_nw+1; j++) {
			if(j < 0)
				continue;
			if (j>dpLen)
				break;
				
			nm[i*dp_size + j] = NEG_INF;
		}
	}
	
	//Initialize the first row to have the gap penalty.
	for(i=0; i<=max_gap_nw+2; i++)
		nm[i*dp_size] = gGAP * i;
	
	//Initialize the second half column to have the gap penalty.
	for(j=0; j<=max_gap_nw+2; j++)
		nm[j] = gGAP * j;


	//I. first half part
	for (i=1;i<iHalf;i++) { //finish DP right after reaching to the bottom of row
		for( j=i-max_gap_nw; j<=i+max_gap_nw; j++) {
			if(j <= 0)
				continue;
			if (j>dpLen)
				break;
			
			m_mm1 = nm[(i-1)*dp_size + (j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			gap1 = nm[(i-1)*dp_size + j] + gGAP;
			gap2 = nm[i*dp_size + (j-1)] + gGAP;
			
			max_score = max_flt(m_mm,gap1,gap2);
			nm[i*dp_size + j] = max_score;
		}
	}
	
	//II. half in DP
	for( j=iHalf-max_gap_nw; j<=iHalf+max_gap_nw; j++) {
		if(j <= 0)
			continue;
		if (j>dpLen)
			break;
		
		m_mm1 = nm[(i-1)*dp_size + (j-1)];
		m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
		m_mm = m_mm1 + m_mm2;
		gap1 = nm[(i-1)*dp_size + j] + gGAP;
		gap2 = nm[i*dp_size + (j-1)] + gGAP;
		
		max_score = max_flt(m_mm,gap1,gap2);
		nm[i*dp_size + j] = max_score;
		
		if (max_score>=mid_max_score) {
			mid_max_score = max_score;
		}
	}

	if (mid_max_score < earlyTerminationScore)
		return align_score;
	
	//III. second half part
	for (i=iHalf+1;i<=dpLen;i++) { //finish DP right after reaching to the bottom of row
		for( j=i-max_gap_nw; j<=i+max_gap_nw; j++) {
			if(j <= 0)
				continue;
			if (j>dpLen)
				break;
			
			m_mm1 = nm[(i-1)*dp_size + (j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			gap1 = nm[(i-1)*dp_size + j] + gGAP;
			gap2 = nm[i*dp_size + (j-1)] + gGAP;
			
			max_score = max_flt(m_mm,gap1,gap2);
			nm[i*dp_size + j] = max_score;
		}
	}
	
	//float align_score = nm[0][0];
	align_score = nm[dpLen*dp_size+dpLen];
	// Don't delete.  Keep in memory
	//delete[] nm;
	return align_score;
}

/*
 * We know there's an exact match, so we just compute the alignment score with itself
 * from start to end.
 * start is the first position of the exact match.  end is the last position of the 
 * exact match (not the first position of a possible match).
 * 
 * Note:  This is the Hamming Distance (with fastq modifications)
 */
float bin_seq::get_align_score_mid(const Read &read, const string &consense, 
			const unsigned int start, const unsigned int end) {
	float score = 0;
	for(unsigned int i=start; i<=end; i++) {

		score += (read.pwm[i][0] * gALIGN_SCORES[(unsigned int)consense[i]][0])
		   + (read.pwm[i][1] * gALIGN_SCORES[(unsigned int)consense[i]][1])
		   + (read.pwm[i][2] * gALIGN_SCORES[(unsigned int)consense[i]][2])
		   + (read.pwm[i][3] * gALIGN_SCORES[(unsigned int)consense[i]][3]);

#ifdef DEBUG
		if(VERBOSE > 2) {
		fprintf(stderr, "%c/%c %.4f %.4f %.4f %.4f = %.4f\n",
			read.seq[i], consense[i],
			read.pwm[i][0] * gALIGN_SCORES[(unsigned int)consense[i]][0],
			read.pwm[i][1] * gALIGN_SCORES[(unsigned int)consense[i]][1],
			read.pwm[i][2] * gALIGN_SCORES[(unsigned int)consense[i]][2],
			read.pwm[i][3] * gALIGN_SCORES[(unsigned int)consense[i]][3],
			score);
		fprintf(stderr, "\t%.4f %.4f %.4f %.4f\n",
			read.pwm[i][0],
			read.pwm[i][1],
			read.pwm[i][2],
			read.pwm[i][3]);
		fprintf(stderr, "\talign_score(match=%f, adj=%f) : %.4f %.4f %.4f %.4f\n",
			gMATCH,gADJUST,
			gALIGN_SCORES[(unsigned int)consense[i]][0],
			gALIGN_SCORES[(unsigned int)consense[i]][1],
			gALIGN_SCORES[(unsigned int)consense[i]][2],
			gALIGN_SCORES[(unsigned int)consense[i]][3]);
		}
#endif
	}
	
	return score;
}

/*
 * This is meant to align the last piece of the read.  Aligning from the (exact) matching
 * portion which was produced from the hash piece to the end of the read.
 *
 *  X
 *	  X		<- start.  We know it must be an exact match
 *		*
 *        *
 *          *	<- we want to match this last piece.
 *			  *
 *
 */
float bin_seq::get_align_score_end(const Read &read, const string &gen, const unsigned int start) {
	if(start == read.length-1)	//trivial case
		return 0;
		
	unsigned int length = read.length - start;	//start+length == read.length
	
	int max_gap_nw = (int)floor(gMAX_GAP*length/gSEQ_LENGTH+0.5);
	if (max_gap_nw < gMAX_GAP)
		max_gap_nw = gMAX_GAP;
	
	// Use the default array
	float* nm = def_arr;
	
	//Initialize the middle of the array to zeros.
	for(unsigned int i=0; i<length; i++) {
		for(int j=i-max_gap_nw-1; j<=(signed int)i+max_gap_nw+1 && j < (signed int)length; j++) {	//j needs to be signed
			if(j < 0)
				continue;
			//nm[i][j] = NEG_INF;
			nm[i*length+j] = NEG_INF;
		}
	}
	
	//Initialize the first row to have the gap penalty.
	for(unsigned int i=0; i<=(unsigned int)max_gap_nw+1 && i<length; i++) {
		//nm[i][0] = gGAP * i;
		nm[i*length+0] = gGAP * i;
	}
	//Initialize the first column to have the gap penalty.
	for(unsigned int j=0; j<=(unsigned int)max_gap_nw+1 && j<length; j++) {
		//nm[0][j] = gGAP * j;
		nm[0*length+j] = gGAP * j;
	}
	
	
	float m_mm1=0, m_mm2=0, m_mm=0, gap1=0, gap2=0, max_score=0;
	
	//Set up the NW matrix	
	for(int i=1; i<(int)length; i++) {
		for(int j=i-max_gap_nw; j<=i+max_gap_nw && j<(int)length; j++) {
			if(j <= 0)
				continue;
			if (j>=(int)length)
				break;
			
			//float m_mm1 = nm[i-1][j-1];
			m_mm1 = nm[(i-1)*length+(j-1)];
			m_mm2 = get_val(read.pwm[i-1],gen[j-1]);
			m_mm = m_mm1 + m_mm2;
			//float gap1 = nm[i][j-1] + gGAP;
			gap1 = nm[i*length+(j-1)] + gGAP;
			//float gap2 = nm[i-1][j] + gGAP;
			gap2 = nm[(i-1)*length+j] + gGAP;

			max_score = max_flt(m_mm,gap1,gap2);
			
			//nm[i][j] = max_score;
			nm[i*length+j] = max_score;

		}	
	}
	
#ifdef DEBUG
/*
	cerr << "\nEND\nstart: " << start << endl;
	for(unsigned int i=0; i<length; i++) {
		for(unsigned int j=0; j<length; j++) {
			fprintf(stderr,"%.2f\t",nm[i*length+j]);
			
		}
		cerr << endl;
	}
	cerr << endl;
	fprintf(stderr,"Returning %f\n",nm[(length-1)*length+(length-1)]);
//*/
#endif	

	//float align_score = nm[length-1][length-1];
	float align_score = nm[(length-1)*length+(length-1)];
	// Keep the default array in memory
	//delete[] nm;
	return align_score;
}


inline float bin_seq::get_val(const float* a, const char b) {

	float score = 0;
	float* scoreptr = gALIGN_SCORES[(unsigned int)b];
	
	if (gBISULFITE) {
		if ((b == 'c') && (a[1]<a[3])) {
			scoreptr = gALIGN_SCORES[(unsigned int)'d'];
		}
	}
	else if (gBISULFITE2) {
		if ((b == 'g') && (a[2]<a[0])) {
			scoreptr = gALIGN_SCORES[(unsigned int)'h'];
		}
	}
	else if (gATOG) {
		if ((b=='a') && (a[0]<a[2])){
			scoreptr = gALIGN_SCORES[(unsigned int)'b'];
		}
	}
	else if (gATOG2) {
		if ((b == 't') && (a[3]<a[1])) {
			scoreptr = gALIGN_SCORES[(unsigned int)'u'];
		}
	}
	score = (a[0] * scoreptr[0]
		   + a[1] * scoreptr[1]
		   + a[2] * scoreptr[2]
		   + a[3] * scoreptr[3]);
	return score;
}

float bin_seq::max_flt(char &path, float diag, float upgap, float leftgap, float sw_zero) {
	//cout << diag << ' '<< upgap << ' ' << c << endl;
	if(diag >= upgap) {
		if(diag >= leftgap) {// now upgap out and diag vs. leftgap
			if (diag >= sw_zero){
				path = 'D';
				return diag;
			}
			else{
				path = 'S';
				return sw_zero;
			}
		}
		else {
			if (leftgap >= sw_zero){
				path = 'L';
				return leftgap;
			}
			else{
				path = 'S';
				return sw_zero;
			}
		}
	}
	else {	//we know upgap>diag.
		if(upgap >= leftgap) {
			if (upgap >= sw_zero){
				path = 'U';
				return upgap;
			}
			else{
				path = 'S';
				return sw_zero;
			}
		}
		else {
			if (leftgap >= sw_zero){
				path = 'L';
				return leftgap;
			}
			else{
				path = 'S';
				return sw_zero;
			}
		}
	}
}

float bin_seq::max_flt(char &path, float diag, float upgap, float leftgap) {
	//cout << diag << ' '<< upgap << ' ' << c << endl;
	if(diag >= upgap) {
		if(diag >= leftgap) {
			path = 'D';
			return diag;
		}
		else {
			path = 'L';
			return leftgap;	//implies leftgap>diag, which means leftgap>upgap.
		}
	}
	else {		//we know know upgap>diag.
		if(upgap >= leftgap) {
			path = 'U';
			return upgap;
		}
		else {
			path = 'L';
			return leftgap;
		}
	}
}

inline float bin_seq::max_flt(float a, float b, float c) {
	//cout << a << ' '<< b << ' ' << c << endl;
	if(a >= b) {
		if(a >= c)
			return a;
		else return c;	//implies c>a, which means c>b.
	}
	else {		//we know know b>a.
		if(b >= c)
			return b;
		else
			return c;
	}
}

inline float bin_seq::max_flt(float a, float b, float c, float z) {
	//cout << a << ' '<< b << ' ' << c << endl;
	if(a >= b) {
		if(a >= c) {
			if(a >= z)
				return a;
			else
				return z;
		}	
		else {
			if (c >= z)	//implies c>a, which means c>b.
				return c;
			else
				return z;
		}
	}
	else {	//we know know b>a.
		if(b >= c) {
			if(b >= z)
				return b;
			else
				return z;
		}
		else {
			if (c >= z)
				return c;
			else
				return z;
		}
	}
}


#define CLOSE_ENOUGH 0.01

bool bin_seq::Test_arr_eq(float a[][5], int a_n, int a_m, float** b, int b_n, int b_m) {
	if(a_n != b_n || a_m != b_m)
		return false;	// already know they're not equal
		
	for(int i=0; i<a_n; i++) {
		for(int j=0; j<a_m; j++) {
			if(abs(a[i][j] - b[i][j]) > CLOSE_ENOUGH) {
				fprintf(stderr,"Failed at i=%d,j=%d with a=%f b=%f\n",i,j,a[i][j],b[i][j]);
				return false;
			}
		}
	}
	
	return true;
}

bool bin_seq::Test(ostream &os, unsigned int &warnings) {
	bin_seq bs;
	bool success = true;
	
	Read r1;
	r1.length = 35;
	string consense = "acgtcgatcgtggctaatcgttcgtagatcgatta";
	string gen_str =  "acgtcgatcgtggctaatcgttgtagatccgatta";
	r1.pwm = new float*[35];
	for(unsigned int i=0; i<35; i++) {
		r1.pwm[i] = new float[4];
	}
	unsigned int row = 0;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 1.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 1.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 1.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 1.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 1.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 1.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 1.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 1.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 0.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 1.0; row++;
	r1.pwm[row][0] = 1.0; r1.pwm[row][1] = 0.0; r1.pwm[row][2] = 0.0; r1.pwm[row][3] = 0.0; row++;
	
	pair<string,string> result = bs.get_align_score_w_traceback(r1, consense, gen_str);
	printf("CONSENSE: %s\nGENOME:   %s\nRESULT:   %s\n          %s\n",
		consense.c_str(),gen_str.c_str(),result.first.c_str(),result.second.c_str());

	TEST(result.first == "acgtcgatcgtggctaatcgttggtagat-cgatta");
	//TEST(result.second == "acgtcgatcgtggctaatcgtt-gtagatccgatta");
	TEST(result.second == "22M1I6M1D6M");
	fprintf(stderr,"Expected is 22M1I6M1D6M and actual is %s\n",result.second.c_str());

	//consense = "acgtcgatcgtggctaatcgttcgtagatcgatta";
	// gen_str = "acgtcgatcgtggctaatcgttcggagatcgatta";
	string gen_str2 =  "acgtcgtttatcgtggctaatcgttccgattaccc";
	// Check timing on this alignment for traceback
	/*
	for(unsigned int i=0; i<1000; i++) {
		result = bs.get_align_score_w_traceback(r1, consense, gen_str2);
	}
	*/
	result = bs.get_align_score_w_traceback(r1, consense, gen_str2);
	fprintf(stderr,"CONSENSE: %s\nGENOME:   %s\nRESULT:   %s\n CIGAR:   %s\n",
		consense.c_str(),gen_str2.c_str(),result.first.c_str(),result.second.c_str());
	TEST(result.first == "acgtcg---atcgtggctaatcgttcttggatcaatta");
	TEST(result.second == "6M3D17M1I1M1I4M1I4M");

	double start = When();
	float align_score;
	// Check timing on this alignment for normal
	for(unsigned int i=0; i<1000000; i++) {
		align_score = bs.get_align_score(r1,gen_str,10u,20u);
	}
	double end = When();
	fprintf(stderr,"Total time is %f seconds for 10M alignments, score is %f\n",end-start,align_score/gADJUST);
	TEST(align_score == 26.0f*gADJUST);

	/************PHMM****************/
	float phmm_answer[][5] = {
		{1.00,0.00,0.00,0.00,0.00},
		{0.00,1.00,0.00,0.00,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,1.00,0.00,0.00,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{1.00,0.00,0.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,1.00,0.00,0.00,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{0.00,1.00,0.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{1.00,0.00,0.00,0.00,0.00},
		{1.00,0.00,0.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,1.00,0.00,0.00,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,1.00,0.00,0.00,0.00},
		{0.00,0.00,0.99,0.01,0.00},
		{0.99,0.00,0.00,0.01,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{1.00,0.00,0.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,1.00,0.00,0.00,0.00},
		{0.00,0.00,1.00,0.00,0.00},
		{1.00,0.00,0.00,0.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.00,0.00,0.00,1.00,0.00},
		{0.50,0.00,0.00,0.50,0.00},
		{1.00,0.00,0.00,0.00,0.00}};

	///*
	fprintf(stderr,"PHMM Score sequence:\n");
	fprintf(stderr,"_\tA\tC\tG\tT\n");
	fprintf(stderr,"A\t");
	for(unsigned int i=0; i<4; i++) {
		fprintf(stderr,"%f\t",gPHMM_ALIGN_SCORES[(int)'a'][i]);
	}
	fprintf(stderr,"\nC\t");
	for(unsigned int i=0; i<4; i++) {
		fprintf(stderr,"%f\t",gPHMM_ALIGN_SCORES[(int)'c'][i]);
	}
	fprintf(stderr,"\nG\t");
	for(unsigned int i=0; i<4; i++) {
		fprintf(stderr,"%f\t",gPHMM_ALIGN_SCORES[(int)'g'][i]);
	}
	fprintf(stderr,"\nT\t");
	for(unsigned int i=0; i<4; i++) {
		fprintf(stderr,"%f\t",gPHMM_ALIGN_SCORES[(int)'t'][i]);
	}
	fprintf(stderr,"\nN\t");
	for(unsigned int i=0; i<4; i++) {
		fprintf(stderr,"%f\t",gPHMM_ALIGN_SCORES[(int)'n'][i]);
	}
	fprintf(stderr,"\n\n");
	
	float** phmm = bs.pairHMM(r1, consense, gen_str);
	for(unsigned int i=0; i<r1.length; i++) {
		printf("[%u, %c]\t",i,gen_str[i]);
		for(int j=0; j<5; j++) {
			printf("%3.2f\t",phmm[i][j]);
		}
		printf("\n");
		delete[] phmm[i];
	}

	printf("\n\n");

	gen_str = "acgtcgatcgtggctaatcgttcgagatcgattaa";
	phmm = bs.pairHMM(r1, consense, gen_str);
	for(unsigned int i=0; i<r1.length; i++) {
		printf("[%u, %c]\t",i,gen_str[i]);
		for(int j=0; j<5; j++) {
			printf("%3.2f\t",phmm[i][j]);
		}
		printf("\n");
	}

	// see if they're equal
	TEST_W(Test_arr_eq(phmm_answer,r1.length,5,phmm,r1.length,5));

	//*/

	// Clean up after ourselves
	for(unsigned int i=0; i<35; i++) {
		delete[] r1.pwm[i];
	}
	delete[] r1.pwm;
	for(unsigned int i=0; i<r1.length; i++) {
		delete[] phmm[i];
	}
	delete[] phmm;
	
	
	unsigned int GEN_SIZE = 100000;
	unsigned char* genome = new unsigned char[GEN_SIZE];

	ifstream in;	
	in.open("test/short.fa");
	in.seekg(7);
	in.read((char*)genome, GEN_SIZE);
	
	string dna = "aaaagaaag";
	string dna2 = "aaaaaaaaa";
	string dnat = "ttttttttt";

	unsigned int bs2_t_good = 6270;
	unsigned int bs2_t_bad = 33;
	unsigned int good = 7057;
	bin_seq bs2;
	
	//string loads = "ccccccccc";
	//cout << "highest " <<  bs.str2int(loads) << endl;
	
	gMER_SIZE = 11;
	const char* test_genome = "aaaaaaaaaataaaatNtttt";
	unsigned int pos = 0;
	TEST(bs.get_hash(pos,(unsigned char*)test_genome,21u).second == 3);
	TEST(bs.hash2str(3,11) == "aaaaaaaaaat");
	os << bs.hash2str(3,11) << endl;
	
	gMER_SIZE = 9;
	//fprintf(stderr,"Result is %d and should be 0 with genome=%.30s and pos=%u\n",bs.get_hash(bs2_t_bad,genome).first,genome+bs2_t_bad,bs2_t_bad);
	TEST(!bs.get_hash(bs2_t_bad,genome,bs2_t_bad+50u).first);
	TEST(bs.get_hash(bs2_t_good,genome,bs2_t_good+50u).first);

	//fprintf(stderr,"Result is %lx and should be %lx with genome=%.30s and pos=%u\n",bs.get_hash(good,genome).second,65712,genome+good,good);
	TEST(bs.get_hash(good,genome,good+50u).second == 65712);
	//TEST(bs2.get_hash(dna) == bs.get_hash(bs2_t_good,genome));
	cout << "get hash: " << bs2.get_hash(bs2_t_good,genome,50u).second << endl;
	cout << "dna: " << bs2.get_hash(dna).second << endl;
	TEST(bs.hash2str(514,9) == dna);
	TEST(bs.hash2str(0,9) == dna2);
	TEST(bs.hash2str(262143,9) == dnat);
	TEST(bs.get_hash(dnat).second == 262143);
	TEST(bs.get_hash("caaaagtaa").second == 65712);
	TEST(bs.get_hash("atcgccaga").second == 55624);
	
	delete[] genome;
	in.close();
	return success;
}
