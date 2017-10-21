//version GreeKuP_v1.0
#include <iostream>
using namespace std;
#include <complex>
#include "nablainout.h"
#include "keys.h"
#include <stdlib.h>

int main (int argc, char* argv[])
{
	string programname = "Printnabla";

	nablaclass nabla;
	
	keyarray keyarr;
	int keynum;
	int keycounter;
	
	int nablafound = 0;		//defines whether -nabla key is found
	string nablaname;		//the name of the NABLA file

	int outputfound = 0;		//defines whether -output key is found
	string outputname;		//the name of the output file

	int inputresult = 0;		//used to check the result of the opening of the NABLA file

	keyarr.find_parse_keys(argc, argv);
	keynum = keyarr.keynum;
	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str() == "nabla")
		{
			nablafound = 1;
			nablaname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str() == "output")
		{
			outputfound = 1;
			outputname = keyarr.arr[keycounter].param.str();
		}

	}


	if (nablafound == 0)
	{
		nablaname = "NABLA";
	}
	
	if (outputfound == 0)
	{
		outputname = nablaname + "FORM";
	}

	inputresult = nabla.input_binary(nablaname.c_str());
	if (inputresult)
	{
		cerr << programname << ": ERROR: Nabla file " << nablaname << " was not opened, program terminated" << endl;
		exit(1);
	}
	
	nabla.output_form(outputname.c_str(), "Nabla matrix elements printed: ");

	nabla.memory_remove();

	keyarr.memory_remove();

	return 0;
}
