#ifndef OMEGAINOUT_CPP
#define OMEGAINOUT_CPP

#include <iostream>
#include <fstream>
using namespace std;
#include "omegainout.h"

omegaclass::omegaclass()
{
	this->deltaenum = 0;
	this->omeganum = 0;
	this->deltae = NULL;
	this->omega = NULL;
	this->memoryexists = 0;
}

int omegaclass::memory_create(int newdeltaenum, int newomeganum)
{

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

	this->memoryexists = 1;

	return 0;
}

int omegaclass::input_binary(const char* inputfilename)
{
	ifstream inputfile;

	int idum;

	double ridum;
	double rdeltaenum;
	double romeganum;

	inputfile.open(inputfilename, ios_base::in|ios_base::binary);
	if (! inputfile)
	{
		return 1;
	}
	
	inputfile.seekg(0);
	inputfile.read((char*)&ridum, sizeof(double));
	inputfile.read((char*)&rdeltaenum, sizeof(double));
	inputfile.read((char*)&romeganum, sizeof(double));
	idum = (int)ridum;

	this->memory_create((int)rdeltaenum, (int)romeganum);

	inputfile.seekg(idum);
	inputfile.read((char*)this->deltae, (this->deltaenum) * sizeof(double));

	inputfile.seekg(2 * idum);
	inputfile.read((char*)this->omega, (this->omeganum) * sizeof(double));

	inputfile.close();

	return 0;
}

int omegaclass::output_form(const char* outputfilename, const char* header)
{
	ofstream outputfile;
	int deltaecounter, omegacounter;
	
	if (this->memoryexists != 0)
	{
		outputfile.open(outputfilename, ios_base::out);
		outputfile << header << endl << endl;

		outputfile << "number of deltae points: " << this->deltaenum << endl;
		outputfile << "number of omega points: " << this->omeganum << endl << endl;

		outputfile << "deltae values: " << endl;
		outputfile << "number of point               deltae, eV" << endl;
		for (deltaecounter = 0; deltaecounter < this->deltaenum; deltaecounter++)
		{
			outputfile << deltaecounter << "                               " << this->deltae[deltaecounter] << endl;
		}
		outputfile << endl;
	
		outputfile << "omega values: " << endl;
		outputfile << "number of point               omega, eV" << endl;
		for (omegacounter = 0; omegacounter < this->omeganum; omegacounter++)
		{
			outputfile << omegacounter << "                               " << this->omega[omegacounter] << endl;
		}
		outputfile << endl;

		outputfile.close();
	}
	return 0;
}

int omegaclass::output_binary(const char* outputfilename)
{
	ofstream outputfile;
	int idum;

	double convert;

	if (this->memoryexists != 0)
	{
		outputfile.open(outputfilename, ios_base::out|ios_base::binary);

		idum = 3;
		if (this->deltaenum > idum)
			idum = this->deltaenum;
		if (this->omeganum > idum)
			idum = this->omeganum;	
		
		idum = idum * sizeof(double);

		outputfile.seekp(0);
		convert = (double)idum;
		outputfile.write((char*)&convert, sizeof(double));
		convert = (double)this->deltaenum;
		outputfile.write((char*)&convert, sizeof(double));
		convert = (double)this->omeganum;
		outputfile.write((char*)&convert, sizeof(double));

		outputfile.seekp(idum);
		outputfile.write((char*)this->deltae, (this->deltaenum) * sizeof(double));

		outputfile.seekp(2 * idum);
		outputfile.write((char*)this->omega, (this->omeganum) * sizeof(double));

		outputfile.close();
	}

}

int omegaclass::memory_remove()
{
	if (this->memoryexists != 0)
	{
		delete [] this->omega;
		this->omeganum = 0;

		delete [] this->deltae;
		this->deltaenum = 0;

		this->memoryexists = 0;
	}
	return 0;
}

omegaclass::~omegaclass()
{
}

#endif
