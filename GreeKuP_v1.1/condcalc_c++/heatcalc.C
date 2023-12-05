#include <iostream>
using namespace std;
#include <fstream>
#include <sstream>
#include <stdlib.h>
#include "myinout.h"
#include "keys.h"

int main (int argc, char* argv[])
{
	char* programname = "Heatcalc";

	keyarray keyarr;
	int keynum;
	int keycounter;
	
	int resigmafound = 0;
	string resigmaname = "";		//name of input ReCONDUCTIVITY (binary) file

	int L12found = 0;
	string L12name = "";			//name of input L12 file
	int L21found = 0;
	string L21name = "";			//name of input L21 file
	int L22found = 0;
	string L22name = "";			//name of input L22 file

	int thermalfound = 0;
	string thermalname = "";
	int thermalformfound = 0;
	string thermalformname = "";

	int ratiofound = 0;
	string rationame = "";
	int ratioformfound = 0;
	string ratioformname = "";

	conductivity resigma;
	conductivity L12;
	conductivity L21;
	conductivity L22;
	
	conductivity thermal;
	conductivity ratio;

	keyarr.find_parse_keys(argc, argv);
	keynum = keyarr.keynum;
	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str()=="resigma")
		{
			resigmafound = 1;
			resigmaname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="L12")
		{
			L12found = 1;
			L12name = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="L21")
		{
			L21found = 1;
			L21name = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="L22")
		{
			L22found = 1;
			L22name = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="thermal")
		{
			thermalfound = 1;
			thermalname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="thermalform")
		{
			thermalformfound = 1;
			thermalformname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="ratio")
		{
			ratiofound = 1;
			rationame = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="ratioform")
		{
			ratioformfound = 1;
			ratioformname = keyarr.arr[keycounter].param.str();
		}
		else
		{
			cerr << programname << ": ERROR: Unknown key -" << keyarr.arr[keycounter].name.str() << " found; program terminated" << endl;
			exit(1);
		}
	}	

	if (resigmafound == 0)
	{
		cerr << programname << ": ERROR: Necessary key -resigma is not present; input ReCONDUCTIVITY file is not specified; program termionated" << endl; 
		exit(1);
	}

	if (L12found == 0)
	{
		cerr << programname << ": ERROR: Necessary key -L12 is not present; input L12 file is not specified; program terminated" << endl;
		exit(1);
	}
	if (L21found == 0)
	{
		cerr << programname << ": ERROR: Necessary key -L21 is not present; input L21 file is not specified; program terminated" << endl;
		exit(1);
	}
	if (L22found == 0)
	{
		cerr << programname << ": ERROR: Necessary key -L22 is not present; input L22 file is not specified; program terminated" << endl;
		exit(1);
	}

	if (thermalfound == 0)
	{
		thermalname = "THERMAL";
	}
	if (thermalformfound == 0)
	{
		thermalformname = thermalname + "FORM";
	}

	if (ratiofound == 0)
	{
		rationame = "RATIO";
	}
	if (ratioformfound == 0)
	{
		ratioformname = rationame + "FORM";
	}

	if (resigma.input_binary(resigmaname.c_str()) != 0)
	{
		cerr << programname << ": ERROR: Input file " << resigmaname << " was not opened; program terminated" << endl;
		exit(1);
	}
	if (L12.input_binary(L12name.c_str()) != 0)
	{
		cerr << programname << ": ERROR: Input file " << L12name << " was not opened; program terminated" << endl;
		exit(1);
	}
	if (L21.input_binary(L21name.c_str()) != 0)
	{
		cerr << programname << ": ERROR: Input file " << L21name << " was not opened; program terminated" << endl;
		exit(1);
	}
	if (L22.input_binary(L22name.c_str()) != 0)
	{
		cerr << programname << ": ERROR: Input file " << L22name << " was not opened; program terminated" << endl;
		exit(1);
	}

	thermal.copy_frame(resigma);
	thermal = L22 - L12 * L21 / resigma;
	
	thermal.output_binary(thermalname.c_str());
	thermal.output_form(thermalformname.c_str(), "Thermal conductivity calculated:", "Thermal conductivity, W / (m * K)");

	thermal.memory_remove();

	ratio.copy_frame(resigma);

	ratio = L12 * L21 / resigma;
	ratio = ratio / L22;

	ratio.output_binary(rationame.c_str());
	ratio.output_form(ratioformname.c_str(), "Ratio (L12 * L21) / (L11 * L22) calculated: ", "         Ratio");

	ratio.memory_remove();

	resigma.memory_remove();
	L12.memory_remove();
	L21.memory_remove();
	L22.memory_remove();

	keyarr.memory_remove();

	return 0;
}
