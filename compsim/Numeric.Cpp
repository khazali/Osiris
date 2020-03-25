#ifndef NUMERIC_H
#define NUMERIC_H

//#include "Globals.h"
#include "Functions.h"
typedef union {
	long long i64;
	double d64;
} dbl_64;

//FType machine_eps (void);
//FType machine_epsv(FType);


FType MachineEpsilon(void) {
	FType machEps = 1.0;

	do {
		machEps /= 2.0;
	} while ((FType)(1.0 + (machEps / 2.0)) != 1.0);

	return machEps;
}


FType machine_eps(void) {
	dbl_64 s;
	s.d64 = 1.0;	//value
	s.i64++;
	return s.d64 - 1.0;
}

FType machine_epsv(FType myValue) {
	dbl_64 s;
	s.d64 = myValue;	//value
	s.i64++;
	return (s.i64 < 0 ? myValue - s.d64 : s.d64 - myValue);
}

#endif