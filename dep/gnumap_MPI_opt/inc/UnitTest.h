#ifndef UNIT_TEST_H
#define UNIT_TEST_H

#define TEST(cond) \
	if (!(cond)) { \
		success = false; \
		os << "!!  Test FAILED at [" << __FILE__ << ", " << __LINE__ << "]  !!" << endl; \
	} \
	else { \
		os << "==  Test Passed at [" << __FILE__ << ", " << __LINE__ << "]  ==" << endl; \
	}

#define TEST_W(cond) \
	if (!(cond)) { \
		warnings++;	\
		os << "!!  Test WARNING at [" << __FILE__ << ", " << __LINE__ << "]  !!" << endl; \
	} \
	else { \
		os << "==  Test Passed at [" << __FILE__ << ", " << __LINE__ << "]  ==" << endl; \
	}

#endif

