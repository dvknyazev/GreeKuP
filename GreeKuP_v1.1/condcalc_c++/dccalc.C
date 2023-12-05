#include <iostream>
using namespace std;
#include <sstream>
#include <fstream>
#include <stdlib.h>
#include "myinout.h"
#include "keys.h"

int main (int argc, char* argv [])
{
	string programname = "Dccalc";

	conductivity inputcond, outputcond;

	keyarray keyarr;
	int keynum = 0;
	int keycounter;

	int inputpresent = 0;	//shows whether -input=inputname key is present
	string inputname;	//the name of the input conductivity file
	int inputresult = 0;	//0 - input was successful, 1 - othercase

	int dcpresent = 0;		//shows whether -dc key is present
	string dcname;		//the name of the output conductivity file with DC conductivity value

	int dcformpresent = 0;	//shows whether -dcform key is present
	string dcformname;	//the name of the the output formatted conductivity file with DC conductivity value

	int point1;
	int point2;
	
	int deltaecounter;

	keyarr.find_parse_keys(argc, argv);
	keynum = keyarr.keynum;
	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str()=="input")
		{
			inputpresent = 1;
			inputname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="dc")
		{
			dcpresent = 1;
			dcname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="dcform")
		{
			dcformpresent = 1;
			dcformname = keyarr.arr[keycounter].param.str();
		}
		else
		{
			cerr << programname << ": ERROR: Unknown key -" << keyarr.arr[keycounter].name.str() << " found, program terminated" << endl;
			exit(1);
		}
	}

	if (inputpresent == 0)
	{
		cerr << programname << ": ERROR: -input key should be present, program terminated" << endl;
		exit(1);
	}

	if (dcpresent == 0)
	{
		dcname = inputname + "_DC";
	}

	if (dcformpresent == 0)
	{
		dcformname = dcname + "FORM";
	}

	inputresult=inputcond.input_binary(inputname.c_str());
	if (inputresult)
	{
		cerr << programname << ": ERROR: input file " << inputname << " was not opened, program terminated" << endl;
		exit(1);
	}
	
	outputcond.memory_create(inputcond.deltaenum, 3);

	point1 = 0;
	point2 = 1;

	for (deltaecounter = 0; deltaecounter < inputcond.deltaenum; deltaecounter++)
	{
		outputcond.deltae[deltaecounter] = inputcond.deltae[deltaecounter];
		outputcond.omega[0] = 0;
		outputcond.sigma[deltaecounter][0] = inputcond.sigma[deltaecounter][point1] + ((inputcond.sigma[deltaecounter][point2] - inputcond.sigma[deltaecounter][point1]) / (inputcond.omega[point2] - inputcond.omega[point1])) * (0 - inputcond.omega[point1]);
		outputcond.omega[1] = inputcond.omega[point1];
		outputcond.sigma[deltaecounter][1] = inputcond.sigma[deltaecounter][point1];
		outputcond.omega[2] = inputcond.omega[point2];
		outputcond.sigma[deltaecounter][2] = inputcond.sigma[deltaecounter][point2];
	}


	outputcond.output_binary(dcname.c_str());
	outputcond.output_form(dcformname.c_str(), "DC value calculated:", "conductivity, 1 / (Ohm * m)");
	
	outputcond.memory_remove();
	inputcond.memory_remove();

	keyarr.memory_remove();

	return 0;
}
