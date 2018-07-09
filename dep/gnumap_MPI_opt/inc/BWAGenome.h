/*
 *  BWAGenome.h
 *
 *  Created by Jamison Dance on 12/4/09.
 *
 *  Subclass of GenomeMem to interface with the bwa code. It allows
 *  you to align a sequence against a bwa index and get a list of locations
 *  in the genome back.
 */

#ifndef BWAGENOME_H
#define BWAGENOME_H

#include "GenomeMem.h"
#include "UnitTest.h"
#include "Reader.h"
#include "bwtaln.h"

//Just a wrapper to call methods on the bwa code
class BWAGenome : public GemomeMem {
public:
	BWAGenome();
	
	BWAGenome(const char * bwaIndexLocation);
	
	~BWAGenome();
	string GetString(const unsigned long begin, const unsigned int size);
	HashLocation* GetMatches(unsigned int hash);
private:
	bwt_t *bwt[2];
};
#endif