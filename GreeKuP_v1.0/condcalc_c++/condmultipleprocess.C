#include <iostream>
using namespace std;
#include <fstream>
#include <sstream>
#include <math.h>
#include <stdlib.h>
#include "myinout.h"
#include "keys.h"

int main (int argc, char* argv[])
{
	char* programname = "Condmultipleprocess";

	int filesreadreport = 0; //defines whether report about reading CONDUCTIVITYFILE files should be printed
	int avercondreport = 1; //defines whether report about conductivity averaging should be printed

	int filecounter; //input CONDUCTIVITYFILE counter

	int averfound = 0;		//defines whether input files should be averaged, default is zero (-aver=avername key)
	string avername = "";		//used for the name of output file with average conductivity (binary)
	string avernameform = "";	//used for the name of output file with average conductivity (formatted)

	int namebasefound = 0;		//defines whether -namebase key is found and names may be formed on the base of namebase;
	string namebase = "";		//the base of file name; other names (by default) will be formed with suffixes _AVER, _DISP

	int dispfound = 0;		//defines whether dispersion (of single file) should be calculated (-disp=dispname key)
	string dispname = "";		//used for the name of output file with dispersion (of single file) (binary)
	string dispnameform = "";	//used for the name of output file with dispersion (of single file) (formatted)

	int dispaverfound = 0;		//defines whether dispersion (of average value) should be calculated
					// (-dispaver=dispavername key)
	string dispavername = "";	//used for the name of output file with dispersion (of average value) (binary)
	string dispavernameform = "";	//used for the name of output file with dispersion (of average value) (formatted)

	int reldispaverfound = 0;	//defines whether relative dispersion (of average value) should be calculated
					// (reldispaver=dispavername key)
	string reldispavername = "";	//used for the name of output file with relative dispersion (of average value) (binary)
	string reldispavernameform = "";	//used for the name of output file with relative dispersion (of average value) (formatted)
	
	int wsumfound = 0;		//defines whether weighted sum of input files should be calculated (-wsum=wsumname key)
	string wsumname;		//used for the name of output file with weighted sum (binary)
	string wsumnameform;		//used for the name of output file with weighted sum (formatted)
	int weightfound = 0;		//defines whether input file with weights should be used (-weight=weightname)
					//default weights will equal 1
	string weightname;		//used for the name of input text file with weights
	ifstream weightfile;		//file stream for weight file
	int weightnum;			//number of weights in weightfile (must coincide with the number of input files)

	conductivity outputcond;
	conductivity avercond;
	conductivity dispcond;
	conductivity dispavercond;
	conductivity tempcond1;
	conductivity reldispavercond;

	keyarray keyarr;
	int keynum;		//number of keys
	int keycounter;

	double temp1;

	keyarr.find_parse_keys(argc, argv);
	keynum = keyarr.keynum;
	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str()=="aver")
		{
			averfound = 1;
			avername = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="disp")
		{
			dispfound = 1;
			dispname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="dispaver")
		{
			dispaverfound = 1;
			dispavername = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="reldispaver")
		{
			reldispaverfound = 1;
			reldispavername = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="wsum")
		{
			wsumfound = 1;
			wsumname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="weight")
		{
			weightfound = 1;
			weightname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="namebase")
		{
			namebasefound = 1;
			namebase = keyarr.arr[keycounter].param.str();
		}
		else
		{
			cerr << programname << ": Unknown key " << *(argv + keycounter + 1) << "  Program terminated" << endl;
			exit(1);
		}
	}

	if (argc < (2 + keynum))
	{
		cerr << programname << ": At least one input CONDUCTIVITYFILE should be present, program terminated" << endl;
		exit(1);
	}
	
	if (namebasefound == 0)
	{
		namebase = "ReCONDUCTIVTY";
	}

	if ((averfound) && (avername == ""))
	{
		avername = namebase + "_AVER";
	}

	if ((dispfound) && (dispname == ""))
	{
		dispname = namebase + "_DISP";
	}

	if ((dispaverfound) && (dispavername == ""))
	{
		dispavername = namebase + "_DISPAVER";
	}

	if ((reldispaverfound) && (reldispavername == ""))
	{
		reldispavername = namebase + "_RELDISPAVER";
	}

						//initialising of an array of input conductivity objects
	conductivity* inputcond = new conductivity[argc - 1 - keynum];
						//reading from the files tgo input conductivty objects
	for (filecounter = 0; filecounter < (argc - keynum - 1); filecounter++)
	{
		if (inputcond[filecounter].input_binary(*(argv + filecounter + 1 + keynum)) != 0)
		{
			cerr << programname <<": ERROR: Input file " << argv[filecounter + 1 + keynum] << " was not opened, program terminated" << endl;
			exit (1);
		}
						//checking the dimensions of the files
		if (filecounter != 0)
		{
			if (inputcond[filecounter].check_dimensions(inputcond[0]) == 0)
			{
				cerr << programname << ": Dimensions of " << *(argv + 1 + keynum) << "do not coincide with dimensions of " << *(argv + filecounter + 1 + keynum) << "   program terminated" << endl;
				exit(1);
			}
		}
	}

	if ((averfound) || (dispfound) || (dispaverfound))
	{
		avercond.copy_frame(inputcond[0]);
		for (filecounter = 0; filecounter < (argc - 1 - keynum); filecounter++)
		{
			avercond = avercond + inputcond[filecounter];
		}
		avercond = div_double(avercond, (argc - 1 - keynum));
		if (averfound)
		{
			avercond.output_binary(avername.c_str());
			avernameform = avername + "FORM";
			avercond.output_form(avernameform.c_str(), "Averaged:", "conductivity");
		}
		if ((dispfound) || (dispaverfound) || (reldispaverfound))
		{
			if ((argc - 1 - keynum) == 1)
			{
				cerr << programname << ": WARNING: Dispersion can not be calculated with only one input conductivity file" << endl;
			}
			else
			{
				dispcond.copy_frame(inputcond[0]);
				tempcond1.copy_frame(inputcond[0]);
				for (filecounter = 0; filecounter < (argc - 1 - keynum); filecounter++)
				{
					tempcond1 = inputcond[filecounter] - avercond;
					tempcond1 = tempcond1 * tempcond1;
					dispcond = dispcond + tempcond1;
				}
				tempcond1 = div_double(dispcond, (argc - 1 - keynum - 1));
				dispcond.square_root(tempcond1);
				if (dispfound)
				{
					dispcond.output_binary(dispname.c_str());
					dispnameform = dispname + "FORM";
					dispcond.output_form(dispnameform.c_str(), "Dispersion of single element calculated:", "dispersion");
				}

				if ((dispaverfound) || (reldispaverfound))
				{
					dispavercond.copy_frame(dispcond);
					dispavercond = div_double(dispcond, sqrt(argc - 1 - keynum));
					if (dispaverfound)
					{
						dispavercond.output_binary(dispavername.c_str());
						dispavernameform = dispavername + "FORM";
						dispavercond.output_form(dispavernameform.c_str(), "Dispersion of average calculated:", "dispersion of average");
					}

					if (reldispaverfound)
					{
						reldispavercond.copy_frame(dispavercond);
						reldispavercond = dispavercond / avercond;
						reldispavercond.output_binary(reldispavername.c_str());
						reldispavernameform = reldispavername + "FORM";
						reldispavercond.output_form(reldispavernameform.c_str(), "Relative dispersion of average calculated:", "Relative dispersion of average");
						reldispavercond.memory_remove();
					}
					dispavercond.memory_remove();
				}

				dispcond.memory_remove();
				tempcond1.memory_remove();
			}
		}		
		avercond.memory_remove();
	}

	if (wsumfound)
	{
		double* weight = new double[argc - 1 - keynum];
		if (weightfound)
		{
			weightfile.open(weightname.c_str(), ios_base::in);
			if (! weightfile)
			{
				cerr << programname << ": file with weights" << weightname << " was not opened, program terminated" << endl;
				exit(1);
			}
			weightnum = 0;
			while (weightfile >> temp1)
				weightnum++;
			if (weightnum != (argc - 1 - keynum))
			{
				cerr << programname << ": Number of weights in WEIGHT file does not coincide with number of input CONDUCTIVITY files; program terminated" << endl;
				exit(1);
			}
			
			weightfile.clear();
			weightfile.seekg(0);
			for (filecounter = 0; filecounter < (argc - 1 - keynum); filecounter++)
			{
				weightfile >> weight[filecounter];
			}
			weightfile.close();

		}
		else
		{
			for (filecounter = 0; filecounter < (argc - 1 - keynum); filecounter++)
				weight[filecounter] = 1;
		}

		outputcond.copy_frame(inputcond[0]);
		for (filecounter = 0; filecounter < (argc - 1 - keynum); filecounter++)
		{
			outputcond = outputcond + mul_double(inputcond[filecounter], weight[filecounter]);
		}

		outputcond.output_binary(wsumname.c_str());
		wsumnameform = wsumname + "FORM";
		outputcond.output_form(wsumnameform.c_str(), "Weighted sum of conductivities:", "conductivity (1 / (Ohm * m))");

		outputcond.memory_remove();			
		delete [] weight;
	}

	for (filecounter = 0; filecounter < argc - 1 - keynum; filecounter++)
	{
		inputcond[filecounter].memory_remove();
	}
	delete [] inputcond;

	keyarr.memory_remove();

	return 0;
}

