// ACCIV.cpp : Defines the entry point for the console application.
//

#include "ACCIVPass.h"

int main(int argc, char * argv[])
{

	ACCIVPass pass;
	UString folder(".");
	if(argc > 1)
		folder = argv[1];

	bool result = pass.doPass(folder);
	if(!result)
		fprintf(stderr, "Pass failed!\n");

	return result? 0 : 1;
}

