#ifndef NORMAL_SCORED_SEQ_H
#define NORMAL_SCORED_SEQ_H

#include "ScoredSeq.h"

class NormalScoredSeq : public ScoredSeq {
	
	public:
		NormalScoredSeq() :ScoredSeq() {}
	
		NormalScoredSeq(string seq, double as, unsigned long pos, bool us, int am, bool pe1, bool pe2) 
			: ScoredSeq(seq,as,pos,us,am,pe1,pe2) {}
		
		NormalScoredSeq(const ScoredSeq &other) : ScoredSeq(other) {}
		
		NormalScoredSeq& operator =(const ScoredSeq &other) {
			Init(other);
			return *this;
		}
		
		virtual ~NormalScoredSeq() {}
		
		/*
		 * score will add the log-odds score (adjusted by denom) to each
		 * position on the Genome gen.  Additionally, this can be passed
		 * a flag to specify the number of positions after the start position
		 * to add to (for adding to different lengths of sequences).
		 */
		virtual void score(double, Genome &, unsigned int, Read &, map<unsigned long, unsigned int> &, double &, pthread_mutex_t&);
		
		//virtual void score2(double, double &, Genome &, unsigned int, Read &, pthread_mutex_t&);
		
};

#endif

