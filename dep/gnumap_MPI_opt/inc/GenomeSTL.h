#ifndef GENOME_STL_H
#define GENOME_STL_H

#ifdef	MPI_RUN
#include <mpi.h>
#endif

#include <map>
#include <string>
#include <vector>
#include <fstream>

#include "UnitTest.h"
#include "gvector.h"
#include "bin_seq.h"
#include "const_include.h"
#include "Reader.h"
#include "Genome.h"

typedef unsigned long LOC_ARR_TYPE;


/* 
 * The GenomeSTL class is used as an extension of the Genome with an STL hash
 */

class GenomeSTL : public Genome {
	
	public:
		/*
		 * Constructor.  Creates new Genome with filename.
		 * fn - the filename of the file to be parsed into the genome.
		 */
		GenomeSTL(const char* fn);
		GenomeSTL();
		
		~GenomeSTL();
		
		/*
		 * GetMatches will return the vector that represents the given hash.
		 *
		 * @return the vector of the given hash.
		 */
		HashLocation* GetMatches(unsigned int hash);
		HashLocation* GetMatches(string &);
		
		/*
		 * hash_and_store() will cause the genome to be hashed--from the file fn given above.
		 * it will then store the genome in a smaller binary array.
		 */
		void hash_and_store();	

		//! Mostly used in testing...
		void print_to_file(map<MAP_TYPE, HashLocation> & gh, char* ofn);

	private:
		
		/*!
		 * fix_hash will remove all the high-density matches in the hash.  Anything over the
		 * threshold (MAX_HASH_SIZE) will be removed.
		 */
		void fix_hash();
		void fix_hash(map<MAP_TYPE, vector<unsigned long> > *,
							 gvector<GEN_TYPE> *);

		// Reads the genome file from memory
		void readHash(FILE* readF, bool ldebug);
		// Saves the genome file to memory
		void saveHash(FILE* saveF, bool ldebug);
		
		/**
		 * convertToVector will convert the hashed-and-stored genome to the non-STL version
		 *
		 * @post gh and g will be deleted
		 */
		void convertToVector(map<MAP_TYPE, vector<unsigned long> > *gh,
		                     //gvector<GenomeLocation> *g);
		                     gvector<GEN_TYPE> *pg);

		//gen_piece is used to store a portion of the genome as we read it in.
		//unsigned char* gen_piece;
		map<MAP_TYPE, HashLocation> gen_hash;
		LOC_ARR_TYPE* hash_data;
		unsigned long num_hash_elements;
};

#endif //GENOME_STL_H
