#include "Reader.h" 

/*!
 * Constructor.  Creates new Genome with filename.
 * @param fn - the filename of the file to be parsed into the genome.
 */
Reader::Reader(const char* fn) :
	filename(0) {
	Init(fn);
}

Reader::Reader() :
	filename(0), file_size(0) {
	//Init("");
}

Reader::~Reader() {
	delete[] filename;
}

void Reader::Init(const char* fn) {
	if(filename)
		delete[] filename;
	filename = new char[strlen(fn)+1];
	strcpy(filename,fn);

	buffer_read = 0;
	finished_reading = false;
	
	ifstream infile;
	infile.open(fn);
	if(!infile) {
		perror("");
		char* errstr = new char[150];
		strcpy(errstr,"Bad File: ");
		strncat(errstr,fn,150-strlen(errstr));
		throw errstr;
	}
	
	getline(infile,chr_name);
	if(chr_name[0] == '>') {
		chr_name = chr_name.substr(1,chr_name.size());
		if(gVERBOSE)
			cout << "Reading: " << chr_name << endl;

		current_offset = chr_name.size() + 1;
	}
	else {
		current_offset = 0;
		chr_name = filename;
	}
	
	FILE * fle = fopen(fn,"r");
	fseek(fle,0,SEEK_END);
	file_size = ftell(fle);
	fclose(fle);
	
	infile.close();
}

long Reader::GetFileSize() const {
	return file_size;
}

void Reader::use(const char* fn) {
	Init(fn);
}

Reader& Reader::operator =(const Reader &other) {
	filename = other.filename;
	buffer_read = other.buffer_read;
	current_offset = other.current_offset; 
	finished_reading = other.finished_reading;
	chr_name = other.chr_name;
	
	return *this;
}


bool Reader::read(unsigned char genome[]) {
  
	if (finished_reading) 
    	return false;

	ifstream in(filename);

	if (!in)
		return false;
  
	//genome = ""; 
	in.seekg(current_offset);

	in.read((char*)genome, gBUFFER_SIZE);
	buffer_read = in.gcount();

  	/*if(current_offset < 5500000 && current_offset > 4800000) {
  		cout << "current offset: " << current_offset << endl;
  		in.seekg(current_offset - 10);
		char buff[110];
  		in.read(buff, 110);
  		buff[110] = '\0';
		//strncpy(buff,genome,9);
		//buff[10] = '\0';
		//cout << "Buffer: " << buff << endl;
		//for(int i=0; i<110; i++)
		//	cout << (int)buff[i] <<' ';
		//cout << endl;
	}*/
	
  	// If it didn't pull out gBUFFER_SIZE characters, we can assume it has reached an eof
  	// place an \0 at the end so it will be recognized by the hash code.
  	//if (buffer_read < gBUFFER_SIZE)
		//genome[buffer_read+1] = '\0';
		genome[buffer_read] = '\0';
	
	if(in.eof()) {
		finished_reading = true;
		return false;
	}
  	
  	in.close();
  	
	//current_offset needs to ensure that we have the MER_SIZE bases at the end of each read.
	//The next read will begin and catch the MER_SIZE mers.
	current_offset += buffer_read-gMER_SIZE;
	return true;
}

void Reader::ShiftOffset(int to_shift) {
	current_offset += to_shift;
}

/**
 * ReadLine will read the next line, starting at the current offset
 */
char* Reader::ReadLine() {
	int LINE_SZ = 1024;

	char* buffer = new char[LINE_SZ+2];

	ifstream in(filename);
	in.seekg(current_offset+gMER_SIZE);

	buffer[0] = '\0';
	in.getline(buffer,LINE_SZ);
	if(in.fail()) {
		buffer[LINE_SZ] = '\0';
		printf("WARNING: Reading Failed!\n");
		return buffer;
	}
	
	in.close();
	streamsize bases_read = in.gcount();
	if(gVERBOSE > 1)
		printf("Read %d bases just barely.  Going to increment by %d: %s\n",
			(int)bases_read,(int)(strlen(buffer) + gMER_SIZE),buffer);

	//we want to increment the current offset the number of bases we've just read.  But the
	//next time around, we don't want to return those same bases, so we need to make up for
	//missing bases the last time.
	current_offset += gMER_SIZE + bases_read;	
	return buffer;
}

string Reader::GetName() const {
	return chr_name;
}

bool Reader::Test(ostream &os, unsigned int &warnings) {
	bool success = true;
	
	os << "**EMPTY**" << endl;
	Reader cr;
	
	return success;
}


