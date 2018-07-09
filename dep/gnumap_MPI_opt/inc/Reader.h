#ifndef READER_H
#define READER_H

#include "UnitTest.h"
#include "const_include.h"
#include <string>
#include <cstring>
#include <iostream>
#include <fstream>
#include <cstdio>

using namespace std;

/* 
 * The ChromReader class is used to store the genome that is loaded in a decent size.
 */
class Reader {
	public:
		Reader();
		~Reader();
		
		/*
		 * Constructor.  Creates new ChromReader with filename.
		 * fn - the filename of the file to be parsed into the genome.
		 */
		Reader(const char* fn);
		
		void use(string str) {
			use(str.c_str());
		}
		void use(const char* fn);
		
		Reader& operator =(const Reader &other);
		bool read(unsigned char* genome);
		
		/*
		 * ShiftOffset is used by the Genome when it's absolutely necessary to shift the offset.  It's not
		 * usually a good move, because it screws up the reading, but in the case of multiple chromosome
		 * names in files, it works out great.
		 * We're including a signed int here so we can shift it up or down as needed.
		 */
		void ShiftOffset(int to_shift);

		/* 
		 * ReadName is intended to read the next line in the input stream
		 * It's used for the Genome when it comes accross a genomie name and needs to get the rest of the name
		 */
		char* ReadLine();
		string GetName() const;
		long GetFileSize() const;
		
		static bool Test(ostream &os, unsigned int &);

	private:
		void Init(const char* fn);

		char* filename;
		BIT64_t buffer_read;
		BIT64_t current_offset;
		bool finished_reading;
		string chr_name;
		long file_size;
		
		ifstream infile;
		FILE *in;

};

#endif
