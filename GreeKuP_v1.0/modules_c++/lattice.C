#include <iostream>
#include <fstream>
using namespace std;
#include <math.h>
#include "lattice.h"

//latticeprint - prints data from LATT_INI structure:
//	scale, direct lattice vectors, reciprocal lattice vectors, cell volume
int latticeprint(latt LATT_INI)
{       
        int a;
        int b;
	cout << "initial lattice parameters (assumed that cell shape and volume do not change):" << endl;
        cout << "scale = " << LATT_INI.SCALE << endl;
        cout << "direct lattice: " << endl;
        cout << "    vector number                        coordinates, A                vector length, A" << endl;
        for (b = 0; b < 3; b++)
        {
        	cout << "         A" << b << "                            ";
                for (a = 0; a < 3; a++)
                {
                	cout << "    " << LATT_INI.A[b][a];
                }
                cout << "                    " << LATT_INI.ANORM[b];
                cout << endl;
        }
        cout << "reciprocal lattice: " << endl;
        cout << "    vector number                        coordinates/2*pi, A^-1     vector length/2*pi, A^-1" << endl;
        for (b = 0; b < 3; b++)
        {
        	cout << "         B" << b << "                            ";
                for (a = 0; a < 3; a++)
                {
            		cout << "    " << LATT_INI.B[b][a];
                }
                cout << "                    " << LATT_INI.BNORM[b];
                cout << endl;
        }
        cout << "cell volume (direct space) = " << LATT_INI.OMEGA << " A^3" << endl << endl;
	return 0;
}

int latticeclean(latt* LATT_INI)
{
	int a;
	int b;
	LATT_INI->SCALE = 0;
	for (b = 0; b < 3; b++)
	{	
		for (a = 0; a < 3; a++)
		{
			LATT_INI->A[b][a] = 0;
		}
		LATT_INI->ANORM[b] = 0;
	}
	
	for (b = 0; b < 3; b++)
	{
		for (a = 0; a < 3; a++)
		{
			LATT_INI->B[b][a] = 0;
		} 
		LATT_INI->BNORM[b] = 0;
	}
	LATT_INI->OMEGA = 0;
	
	return 0;
}

int lattice_set_nan(latt* LATT_INI)
{
	int a;
	int b;
	LATT_INI->SCALE = NAN;
	for (b = 0; b < 3; b++)
	{	
		for (a = 0; a < 3; a++)
		{
			LATT_INI->A[b][a] = NAN;
		}
		LATT_INI->ANORM[b] = NAN;
	}
	
	for (b = 0; b < 3; b++)
	{
		for (a = 0; a < 3; a++)
		{
			LATT_INI->B[b][a] = NAN;
		} 
		LATT_INI->BNORM[b] = NAN;
	}
	LATT_INI->OMEGA = NAN;
	
	return 0;
}

int latticecopy(latt latt1, latt* latt2)
{
	int a;
	int b;

	latt2->SCALE = latt1.SCALE;
	for (b = 0; b < 3; b++)
	{	
		for (a = 0; a < 3; a++)
		{
			latt2->A[b][a] = latt1.A[b][a];
		}
		latt2->ANORM[b] = latt1.ANORM[b];
	}
	
	for (b = 0; b < 3; b++)
	{
		for (a = 0; a < 3; a++)
		{
			latt2->B[b][a] = latt1.B[b][a];
		} 
		latt2->BNORM[b] = latt1.BNORM[b];
	}
	latt2->OMEGA = latt1.OMEGA;
}

int latticeprintfile(latt LATT_INI, ofstream* outputfile)
{
	
        int a;
        int b;
	*outputfile << "initial lattice parameters (assumed that cell shape and volume do not change):" << endl;
        *outputfile << "scale = " << LATT_INI.SCALE << endl;
        *outputfile << "direct lattice: " << endl;
        *outputfile << "    vector number                        coordinates, A                vector length, A" << endl;
        for (b = 0; b < 3; b++)
        {
        	*outputfile << "         A" << b << "                            ";
                for (a = 0; a < 3; a++)
                {
                	*outputfile << "    " << LATT_INI.A[b][a];
                }
                *outputfile << "                    " << LATT_INI.ANORM[b];
                *outputfile << endl;
        }
        *outputfile << "reciprocal lattice: " << endl;
        *outputfile << "    vector number                        coordinates/2*pi, A^-1     vector length/2*pi, A^-1" << endl;
        for (b = 0; b < 3; b++)
        {
        	*outputfile << "         B" << b << "                            ";
                for (a = 0; a < 3; a++)
                {
            		*outputfile << "    " << LATT_INI.B[b][a];
                }
                *outputfile << "                    " << LATT_INI.BNORM[b];
                *outputfile << endl;
        }
        *outputfile << "cell volume (direct space) = " << LATT_INI.OMEGA << " A^3" << endl << endl;
	return 0;
}
