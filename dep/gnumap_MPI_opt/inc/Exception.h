#ifndef EXCEPTION_H
#define EXCEPTION_H

#include <string>
#include <cstring>

using namespace std;

struct Exception {
	public:
		Exception(const string msg) {
			errstr = new char[msg.size() + 1];
			strcpy(errstr,msg.c_str());
		}

		Exception(const char* msg) {
			errstr = new char[strlen(msg)];
			strcpy(errstr,msg);
		}

		Exception(const char* msg, const char* FN, const unsigned int LN) {
			errstr = new char[strlen(msg) + strlen(FN) + 25];
			sprintf(errstr,"!!%s:%u %s",FN,LN,msg);
		}

		~Exception() {
			if(errstr)
				delete[] errstr;
		}

		char* GetMessage() {
			return errstr;
		}

		void SetMessage(const char* msg) {
			errstr = new char[strlen(msg)];
			strcpy(errstr,msg);
		}

	private:
		char* errstr;
};

#endif

