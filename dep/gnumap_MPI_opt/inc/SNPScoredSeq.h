#ifndef SNP_SCORED_SEQ_H
#define SNP_SCORED_SEQ_H

#include "ScoredSeq.h"

class SNPScoredSeq : public ScoredSeq {
	
	public:
		SNPScoredSeq() :ScoredSeq() {}
	
		SNPScoredSeq(string seq, double as, unsigned long pos, bool us, int am, bool pe1, bool pe2) 
			: ScoredSeq(seq,as,pos,us,am,pe1,pe2) {}
		
		SNPScoredSeq(const ScoredSeq &other) : ScoredSeq(other) {}
		
		SNPScoredSeq& operator =(const ScoredSeq &other) {
			Init(other);
			return *this;
		}
		
		virtual ~SNPScoredSeq() {}
	
		/*
		 * score will add the log-odds score (adjusted by denom) to each
		 * position on the Genome gen.  Additionally, this can be passed
		 * a flag to specify the number of positions after the start position
		 * to add to (for adding to different lengths of sequences).
		 */
		virtual void score(double, Genome &, unsigned int, Read &, map<unsigned long, unsigned int> &, double &, pthread_mutex_t&);
		
		
	private:
		inline string GetConsensus(const Read &read) const {
			string seq = "";
		
			for(unsigned int i=0; i<read.length; i++) {
				seq += max_char(read.pwm[i]);
			}
		
			return seq;
		}
		
		/*
		 * Because this is used to produce the consensus sequence, we need ambiguity characters...
		 */
		inline char max_char(const float* chr) const {
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


};

#endif
