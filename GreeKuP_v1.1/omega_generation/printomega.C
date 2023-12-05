#include <iostream>
#include <sstream>
using namespace std;
#include <stdlib.h>
#include "keys.h"
#include "omegainout.h"

int main (int argc, char* argv[])
{
	string programname = "Printomega";

	int omegapresent = 0;		//defines whether -omega key is present
	string omeganame;
	int inputresult = 0;

	int outputpresent = 0;		//defines whether -output key is present
	string outputname;
	
	keyarray keyarr;
	int keynum;
	int keycounter;

	omegaclass omega;

	keyarr.find_parse_keys(argc, argv);
	keynum = keyarr.keynum;
	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str()=="omega")
		{
			omegapresent = 1;
			omeganame = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="output")
		{
			outputpresent = 1;
			outputname = keyarr.arr[keycounter].param.str();
		}
		else
		{
			cerr << programname << ": Unknown key -" << keyarr.arr[keycounter].name.str() << " found, program terminated" << endl;
			exit(1);
		}
	}

	if (omegapresent == 0)
	{
		omeganame = "OMEGA";
	}

	if (outputpresent == 0)
	{
		outputname = omeganame + "FORM";
	}

	inputresult = omega.input_binary(omeganame.c_str());
	if (inputresult)
	{
		cerr << programname << ": ERROR: Omega file " << omeganame << " was not opened, program terminated" << endl;
		exit(1);
	}

	omega.output_form(outputname.c_str(), "Omega values printed:");


	omega.memory_remove();
	keyarr.memory_remove();

	return 0;
}
