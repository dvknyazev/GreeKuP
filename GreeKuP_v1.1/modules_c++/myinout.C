#ifndef MYINOUT_CPP
#define MYINOUT_CPP
#include <iostream>
#include <fstream>
#include <math.h>
using namespace std;
#include <stdlib.h>
#include <cstring>
#include "myinout.h"

conductivity::conductivity()
{
	this->memoryexists = 0;
	this->deltaenum = 0;
	this->omeganum = 0;
}

int conductivity::memory_create(int newdeltaenum, int newomeganum)
{
	int deltaecounter;

	if ((this->memoryexists != 0) && ((this->deltaenum == newdeltaenum) && (this->omeganum == newomeganum)))
	{
		return 0;
	}
	else if ((this->memoryexists != 0) && ((this->deltaenum != newdeltaenum) || (this->omeganum != newomeganum)))
	{
		this->memory_remove();
	}
	
	this->deltaenum = newdeltaenum;
	this->deltae = new double [newdeltaenum];

	this->omeganum = newomeganum;
	this->omega = new double [newomeganum];

	this->sigma = new double* [newdeltaenum];
	for (deltaecounter = 0; deltaecounter < newdeltaenum; deltaecounter++)
	{
		this->sigma[deltaecounter] = new double[newomeganum];
	}

	this->memoryexists = 1;

	return 0;
}

int conductivity::input_binary(const char* inputfilename)
{
	ifstream inputfile;

	int condidum;		//record length in bytes

	double rcondidum;	//same as condidum, but double
	double rdeltaenum;	//same as this->deltaenum, but double
	double romeganum;	//same as this->omeganum, but double

	int deltaecounter;
	int omegacounter;

	inputfile.open(inputfilename, ios_base::in|ios_base::binary);
	if (! inputfile)
	{
		return 1;
	}

	inputfile.seekg(0);
	inputfile.read((char*)&rcondidum, sizeof(double));
	inputfile.read((char*)&rdeltaenum, sizeof(double));
	inputfile.read((char*)&romeganum, sizeof(double));
	condidum = (int)rcondidum;
	
	this->memory_create((int)rdeltaenum, (int)romeganum);

	inputfile.read((char*)(this->deltae), (this->deltaenum)*sizeof(double));
	for (deltaecounter = 0; deltaecounter < (this->deltaenum); deltaecounter++)
	{
		inputfile.seekg(condidum * (deltaecounter + 1));	
		for (omegacounter = 0; omegacounter < (this->omeganum); omegacounter++)
		{
			inputfile.read((char*)(this->omega + omegacounter), sizeof(double));
			inputfile.read((char*)(this->sigma[deltaecounter] + omegacounter), sizeof(double));	
		}
	}

	inputfile.close();

	return 0;
}

int conductivity::output_binary(const char* outputfilename)
{
	int deltaecounter;
	int omegacounter;
	int condidum;		//record lengh in bytes

	ofstream outputfile;
	double convert;

	outputfile.open(outputfilename, ios_base::out|ios_base::binary);
							//CONDUCTIVITY file record length calculation
	if ((this->deltaenum + 3) < (this->omeganum * 2))
        {
                condidum = (this->omeganum) * 2 * sizeof(double);
        }
        else
        {
                condidum = ((this->deltaenum) + 3) * sizeof(double);
        }
	
	outputfile.seekp(0);
	convert = condidum;
        outputfile.write ((char*)&convert, sizeof(double));
        convert = this->deltaenum;
        outputfile.write ((char*)&convert ,sizeof(double));
        convert = this->omeganum;
        outputfile.write ((char*)&convert, sizeof(double));
        outputfile.write((char*)(this->deltae), ((this->deltaenum) * sizeof(double)));

	for (deltaecounter = 0; deltaecounter < (this->deltaenum); deltaecounter++)
        {
                outputfile.seekp((deltaecounter + 1) * condidum);
                for (omegacounter = 0; omegacounter < (this->omeganum); omegacounter++)
                {
                        outputfile.write((char*)(this->omega + omegacounter), sizeof(double));
                        outputfile.write((char*)(this->sigma[deltaecounter] + omegacounter), sizeof(double));
                }
        }


	outputfile.close();
	return 0;
}

int conductivity::output_form(const char* outputfilename, const char* header, const char* valueheader)
{
	int deltaecounter, omegacounter;
	
	ofstream outputfile;
	
	outputfile.open(outputfilename, ios_base::out);
	
	outputfile << header << endl;
        for (deltaecounter = 0; deltaecounter < (this->deltaenum); deltaecounter++)
        {
                outputfile << endl << "deltae = " << this->deltae[deltaecounter] << " eV" << endl << endl;
                outputfile << "omega (eV)               " << valueheader << endl;
                for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
                {
                        outputfile  << this->omega[omegacounter] << "                           " << this->sigma[deltaecounter][omegacounter] << endl;
                }
        }


	outputfile.close();
}

int conductivity::output_form_Origin(const char* outputfilename, const char* xoriginstring, const char* xdimension, const char* yoriginstring, const char* ydimension, const char* columnheader, const char* columndimension)
{
	int deltaecounter, omegacounter;

	ofstream outputfile;

	outputfile.open(outputfilename, ios_base::out);
	
	outputfile << xoriginstring;
	if (strlen(xdimension) != 0)
	{
		outputfile << ", " << xdimension;
	}
	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		outputfile << "\t" << yoriginstring;
		if (strlen(ydimension) != 0)
		{
			outputfile << ", " << ydimension;
		}
	}
	outputfile << endl;

	for (deltaecounter = 0; deltaecounter < this-> deltaenum; deltaecounter++)
	{
		outputfile << "\t";
	}
	outputfile << endl;

	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		outputfile << "\t" << columnheader << " = " << this->deltae[deltaecounter] << " " << columndimension;
	}
	outputfile << endl;

	for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
	{
		outputfile << this->omega[omegacounter];
		for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
		{
			outputfile << "\t" << this->sigma[deltaecounter][omegacounter];
		}
		outputfile << endl;
	}

	outputfile.close();
}

int conductivity::output_form_transpose(const char* outputfilename, const char* header, const char* valueheader)
{
	int deltaecounter, omegacounter;
	
	ofstream outputfile;
	
	outputfile.open(outputfilename, ios_base::out);
	
	outputfile << header << endl;
        for (omegacounter = 0; omegacounter < (this->omeganum); omegacounter++)
        {
                outputfile << endl << "omega = " << this->omega[omegacounter] << " eV" << endl << endl;
                outputfile << "deltae (eV)               " << valueheader << endl;
                for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
                {
                        outputfile  << this->deltae[deltaecounter] << "                           " << this->sigma[deltaecounter][omegacounter] << endl;
                }
        }


	outputfile.close();
}

int conductivity::copy_frame(conductivity cond)
{
	int deltaecounter, omegacounter;

	this->memory_create(cond.deltaenum, cond.omeganum);

	for (deltaecounter = 0; deltaecounter < (this->deltaenum); deltaecounter++)
	{
		this->deltae[deltaecounter] = cond.deltae[deltaecounter];
		for (omegacounter = 0; omegacounter < (this->omeganum); omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = 0;
		}
	}

	for (omegacounter = 0; omegacounter < (this->omeganum); omegacounter++)
	{
		this->omega[omegacounter] = cond.omega[omegacounter];
	}
	return 0;	
}

int conductivity::check_dimensions(conductivity cond)
{
	if ((this->memoryexists == 0) || (cond.memoryexists == 0) || (this->deltaenum != cond.deltaenum) || (this->omeganum != cond.omeganum))
	{
		return 0;
	}
	else
	{
		return 1;
	}
}

conductivity conductivity::operator=(conductivity right)
{
	int deltaecounter, omegacounter;
	if (this->check_dimensions(right) == 0)
	{
		cerr << "Left and right part of assignment operator = for conductivity class have different dimensions, program terminated" << endl;
		exit (1);
	}
	else
	{
		for (omegacounter = 0; omegacounter < right.omeganum; omegacounter++)
		{
			this->omega[omegacounter] = right.omega[omegacounter];
		}
		for (deltaecounter = 0; deltaecounter < right.deltaenum; deltaecounter++)
		{
			this->deltae[deltaecounter] = right.deltae[deltaecounter];
			for (omegacounter = 0; omegacounter < right.omeganum; omegacounter++)
			{
				this->sigma[deltaecounter][omegacounter] = right.sigma[deltaecounter][omegacounter];
			}
		}
		return *this;
	}
}

conductivity operator+(conductivity cond1, conductivity cond2)
{
	conductivity result;
	int deltaecounter, omegacounter;
	
	if (cond1.check_dimensions(cond2)==0)
	{
		cerr <<"Operands of sum operator + for conductivity class have different dimensions, program terminated" << endl;
		exit(1);
	}
	else
	{	
		result.copy_frame(cond1);
		for (deltaecounter = 0; deltaecounter< cond1.deltaenum; deltaecounter++)
		{
			for (omegacounter = 0; omegacounter < cond1.omeganum; omegacounter++)
			{
				result.sigma[deltaecounter][omegacounter] = cond1.sigma[deltaecounter][omegacounter] + cond2.sigma[deltaecounter][omegacounter];
			}
		}
		return result;
	}
}

conductivity operator-(conductivity cond1, conductivity cond2)
{
	conductivity result;
	int deltaecounter, omegacounter;
	
	if (cond1.check_dimensions(cond2)==0)
	{
		cerr <<"Operands of sum operator + for conductivity class have different dimensions, program terminated" << endl;
		exit(1);
	}
	else
	{	
		result.copy_frame(cond1);
		for (deltaecounter = 0; deltaecounter< cond1.deltaenum; deltaecounter++)
		{
			for (omegacounter = 0; omegacounter < cond1.omeganum; omegacounter++)
			{
				result.sigma[deltaecounter][omegacounter] = cond1.sigma[deltaecounter][omegacounter] - cond2.sigma[deltaecounter][omegacounter];
			}
		}
		return result;
	}
}

conductivity operator*(conductivity cond1, conductivity cond2)
{
	conductivity result;
	int deltaecounter, omegacounter;
	
	if (cond1.check_dimensions(cond2)==0)
	{
		cerr <<"Operands of sum operator + for conductivity class have different dimensions, program terminated" << endl;
		exit(1);
	}
	else
	{	
		result.copy_frame(cond1);
		for (deltaecounter = 0; deltaecounter< cond1.deltaenum; deltaecounter++)
		{
			for (omegacounter = 0; omegacounter < cond1.omeganum; omegacounter++)
			{
				result.sigma[deltaecounter][omegacounter] = cond1.sigma[deltaecounter][omegacounter] * cond2.sigma[deltaecounter][omegacounter];
			}
		}
		return result;
	}
}

conductivity operator/(conductivity cond1, conductivity cond2)
{
	conductivity result;
	int deltaecounter, omegacounter;
	
	if (cond1.check_dimensions(cond2)==0)
	{
		cerr <<"Operands of sum operator + for conductivity class have different dimensions, program terminated" << endl;
		exit(1);
	}
	else
	{	
		result.copy_frame(cond1);
		for (deltaecounter = 0; deltaecounter< cond1.deltaenum; deltaecounter++)
		{
			for (omegacounter = 0; omegacounter < cond1.omeganum; omegacounter++)
			{
				result.sigma[deltaecounter][omegacounter] = cond1.sigma[deltaecounter][omegacounter] / cond2.sigma[deltaecounter][omegacounter];
			}
		}
		return result;
	}
}

conductivity mul_double(conductivity cond, double number)
{
	conductivity result;
	int deltaecounter, omegacounter;

	result.copy_frame(cond);
	for (deltaecounter = 0; deltaecounter < cond.deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < cond.omeganum; omegacounter++)
		{
			result.sigma[deltaecounter][omegacounter] = cond.sigma[deltaecounter][omegacounter] * number;
		}
	}
	return result;
}

conductivity div_double(conductivity cond, double number)
{
	conductivity result;
	int deltaecounter, omegacounter;

	result.copy_frame(cond);
	for (deltaecounter = 0; deltaecounter < cond.deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < cond.omeganum; omegacounter++)
		{
			result.sigma[deltaecounter][omegacounter] = cond.sigma[deltaecounter][omegacounter] / number;
		}
	}
	return result;
}


int conductivity::resigma_to_imeps(conductivity realcond)
{
	int deltaecounter, omegacounter;
	this->copy_frame(realcond);
	
	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < this-> omeganum; omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = ((double)1 / (double)13440) * ((realcond.sigma[deltaecounter][omegacounter]) / (this->omega[omegacounter]));
		}
	}

	return 0;
}

int conductivity::imsigma_to_reeps(conductivity imcond)
{
	int deltaecounter, omegacounter;
	this->copy_frame(imcond);
	
	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = 1 - ((double)1 / (double)13440) * ((imcond.sigma[deltaecounter][omegacounter]) / (this->omega[omegacounter]));
		}
	}

	return 0;
}

int conductivity::rerefrac_calc(conductivity reeps, conductivity imeps)
{
	int deltaecounter, omegacounter;
	this->copy_frame(reeps);

	if (reeps.check_dimensions(imeps) == 0)
	{
		cerr << "conductivty::rerefrac_calc: ERROR: reeps and imeps have different dimensions; program terminated" << endl;
		exit(1);
	}
	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = sqrt((sqrt(reeps.sigma[deltaecounter][omegacounter] * reeps.sigma[deltaecounter][omegacounter] + imeps.sigma[deltaecounter][omegacounter] * imeps.sigma[deltaecounter][omegacounter]) + reeps.sigma[deltaecounter][omegacounter]) / 2);
		}
	}

	return 0;
}

int conductivity::imrefrac_calc(conductivity reeps, conductivity imeps)
{
	int deltaecounter, omegacounter;
	this->copy_frame(reeps);

	if (reeps.check_dimensions(imeps) == 0)
	{
		cerr << "conductivty::imrefrac_calc: ERROR: reeps and imeps have different dimensions; program terminated" << endl;
		exit(1);
	}
	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = sqrt((sqrt(reeps.sigma[deltaecounter][omegacounter] * reeps.sigma[deltaecounter][omegacounter] + imeps.sigma[deltaecounter][omegacounter] * imeps.sigma[deltaecounter][omegacounter]) - reeps.sigma[deltaecounter][omegacounter]) / 2);
		}
	}

	return 0;
}

int conductivity::reflect_calc(conductivity rerefrac, conductivity imrefrac)
{
	int deltaecounter, omegacounter;
	this->copy_frame(rerefrac);
	
	if (rerefrac.check_dimensions(imrefrac) == 0)
	{
		cerr << "conductivity::rerefrac_calc: ERROR: rerefrac and imrefrac have different dimensions; program terminated" << endl;
		exit(1);
	}

	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = ((1 - rerefrac.sigma[deltaecounter][omegacounter]) * (1 - rerefrac.sigma[deltaecounter][omegacounter]) + imrefrac.sigma[deltaecounter][omegacounter] * imrefrac.sigma[deltaecounter][omegacounter]) / ((1 + rerefrac.sigma[deltaecounter][omegacounter]) * (1 + rerefrac.sigma[deltaecounter][omegacounter]) + imrefrac.sigma[deltaecounter][omegacounter] * imrefrac.sigma[deltaecounter][omegacounter]);
		}
	}
	return 0;
}

int conductivity::absorp_calc(conductivity imrefrac)
{
	int deltaecounter, omegacounter;
	this->copy_frame(imrefrac);

	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = 10135400 * imrefrac.sigma[deltaecounter][omegacounter] * this->omega[omegacounter];
		}
	}

	return 0;
}

int conductivity::square_root(conductivity inputcond)
{
	int deltaecounter, omegacounter;
	this->copy_frame(inputcond);

	for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
	{
		for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
		{
			this->sigma[deltaecounter][omegacounter] = sqrt(inputcond.sigma[deltaecounter][omegacounter]);
		}
	}

	return 0;
}

int conductivity::check_deltae(double deltaevalue)
{
	int deltaecounter;
	int position = -1;
	
	deltaecounter = 0;
	while ((deltaecounter < this->deltaenum) && (position < 0))
	{
		if (fabs(this->deltae[deltaecounter] - deltaevalue) < 0.0000000001)
			position = deltaecounter;
		deltaecounter++;
	}	
	return position;
}

int conductivity::check_omega(double omegavalue)
{
	int omegacounter;
	int position = -1;
	
	omegacounter = 0;
	while ((omegacounter < this->omeganum) && (position < 0))
	{
//		cout << omegacounter << "  " << this->omega[omegacounter] << endl;
		if (fabs(this->omega[omegacounter] - omegavalue) < 0.0000000001)
			position = omegacounter;
		omegacounter++;
	}	
	return position;
}

long int conductivity::get_arrays_size()
{
	long int size;
	size = (long int)this->deltaenum * (long int)this->omeganum * (long int)sizeof(double) + (long int)this->deltaenum * (long int)sizeof(double) + (long int)this->omeganum * (long int)sizeof(double);
	return size;
}

int conductivity::memory_remove()
{
	int deltaecounter;

	if (this->memoryexists != 0)
	{
		for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
		{
			delete [] this->sigma[deltaecounter];
		}
		delete [] this->sigma;

		delete [] this->omega;
		this->omeganum = 0;

		delete [] this->deltae;
		this->deltaenum = 0;

		this->memoryexists = 0;
	}
	return 0;
}

conductivity::~conductivity()
{
}

#endif
