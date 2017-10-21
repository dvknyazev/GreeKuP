#include <iostream>
#include <sstream>
#include <algorithm>
using namespace std;
#include <keys.h>
#include <math.h>
#include <omegainout.h>
#include <stdlib.h>

int main (int argc, char* argv [])
{
	string programname = "Omega_generate";

	omegaclass omega;
	int deltaenum = 0; 
	int omeganum = 0;

	double* deltae;

	int omegacounter, deltaecounter;

	int omegapresent = 0;	//defines whether -omega key is present
	string omeganame;	//name of file with omega values

	int freqstepfound = 0;
	double freqstep = 0;

	int freqminfound = 0;
	double freqminrequired = 0;
	double freqminreal = 0;

	int freqmaxfound = 0;
	double freqmax = 0;

	keyarray keyarr;
	int keynum;
	int keycounter;
	
	keyarr.find_parse_keys(argc, argv);
	keynum = keyarr.keynum;
	
	deltaenum = 0;
	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str()=="deltae")
		{
			deltaenum++;
		}
	}

	if (deltaenum == 0)
	{
		deltaenum = 4;
		deltae = new double [deltaenum];
		deltae[0] = 0.2;
		deltae[1] = 0.1;
		deltae[2] = 0.05;
		deltae[3] = 0.02;
	}
	else
	{
		deltae = new double [deltaenum];
	}
	

	deltaecounter = 0;

	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str()=="deltae")
		{
			keyarr.arr[keycounter].param >> deltae[deltaecounter];
			deltaecounter++;
		}
		else if (keyarr.arr[keycounter].name.str()=="omega")
		{
			omegapresent = 1;
			omeganame = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="freqstep")
		{
			freqstepfound = 1;
			keyarr.arr[keycounter].param >> freqstep;
		}
		else if (keyarr.arr[keycounter].name.str()=="freqmin")
		{
			freqminfound = 1;
			keyarr.arr[keycounter].param >> freqminrequired;
		}
		else if (keyarr.arr[keycounter].name.str()=="freqmax")
		{
			freqmaxfound = 1;
			keyarr.arr[keycounter].param >> freqmax;
		}
		else
		{
			cerr << programname << ": Unknown key -" <<  keyarr.arr[keycounter].name.str() << " found, program terminated" << endl;
			exit(1);
		}
	}

	if (omegapresent == 0)
	{
		omeganame = "OMEGA";
	}	

	if (freqstepfound == 0)
	{
		freqstep = 0.02;
	}

	if (freqminfound == 0)
	{
		freqminreal = freqstep;
	}
	else
	{
		freqminreal = ceil(freqminrequired / freqstep) * freqstep;
	}

	if (freqmaxfound == 0)
	{
		freqmax = 10;
	}

	
	omeganum = 0;
	while (freqminreal + freqstep * omeganum - freqmax < 1e-10)
		omeganum++;

	
	omega.memory_create(deltaenum, omeganum);
	
	sort (deltae, deltae + deltaenum, greater<double>());
	for (deltaecounter = 0; deltaecounter < deltaenum; deltaecounter++)
	{
		omega.deltae[deltaecounter] = deltae[deltaecounter];
	}

	for (omegacounter = 0; omegacounter < omeganum; omegacounter++)
	{
		omega.omega[omegacounter] = freqminreal + omegacounter * freqstep;
	}

	omega.output_binary(omeganame.c_str());

	omega.memory_remove();

	delete [] deltae;

	keyarr.memory_remove();

	return 0;
}
