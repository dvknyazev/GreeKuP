#ifndef NABLAINOUT_CPP
#define NABLAINOUT_CPP
#include <iostream>
using namespace std;
#include <complex>
#include <fstream>
#include <math.h>
#include "nablainout.h"
#include "myinout.h"
#include "mymath.h"


int get_efermi_Te(int nb_tot, double* celtot, double* fertot, double* efermi, double* Te)
{
	double *y_efermi;       //used to calculate the Fermi energy; values of energy
        double *x_efermi;       //used to calculate the Fermi energy; ln((2 - weight) / weight)
        int npoints_efermi = 0; //the number of points for the Fermi energy calculation; weights must be between 1.5 and 0.5
        int start_efermi = 0;   //the number of the first point with the weight less than 1.5
        int end_efermi = 0; 	//the number of the last point with the weight larger than 0.5
	double upper_boundary_efermi = 0.9;
	double lower_boundary_efermi = 0.1;
	int counter_efermi;
	double xy_aver_efermi, x_aver_efermi, y_aver_efermi, x2_aver_efermi;

	start_efermi = 0;
	while (fertot[start_efermi] > upper_boundary_efermi)
	{
		start_efermi++;
	}
	end_efermi = start_efermi;
	while (fertot[end_efermi] > lower_boundary_efermi)
	{
		end_efermi++;
	}
	npoints_efermi = end_efermi - start_efermi;

	if (npoints_efermi < 2)
	{

		if (start_efermi - 4 < 0)
			start_efermi = 0;
		else
			start_efermi = start_efermi - 4;

		if (end_efermi + 4 > (nb_tot - 1))
			end_efermi = nb_tot;
		else
			end_efermi = end_efermi + 4;

		npoints_efermi = end_efermi - start_efermi;

		if (npoints_efermi < 2)
 		{
			return 0;
                }
		else
		{
		}
	}
	
	x_efermi = new double [npoints_efermi];
	y_efermi = new double [npoints_efermi];

	for (counter_efermi = 0; counter_efermi < npoints_efermi; counter_efermi++)
	{
		x_efermi[counter_efermi] = log ((1 - fertot[start_efermi + counter_efermi]) / fertot[start_efermi + counter_efermi]);
		y_efermi[counter_efermi] = celtot[start_efermi + counter_efermi];
                }

                x_aver_efermi = 0;
                y_aver_efermi = 0;
                x2_aver_efermi = 0;
                xy_aver_efermi = 0;
                for (counter_efermi = 0; counter_efermi < npoints_efermi; counter_efermi++)
                {
                        x_aver_efermi = x_aver_efermi + x_efermi[counter_efermi];
                        y_aver_efermi = y_aver_efermi + y_efermi[counter_efermi];
                        x2_aver_efermi = x2_aver_efermi + x_efermi[counter_efermi] * x_efermi[counter_efermi];
                        xy_aver_efermi = xy_aver_efermi + x_efermi[counter_efermi] * y_efermi[counter_efermi];
                }
                x_aver_efermi = x_aver_efermi / npoints_efermi;
                y_aver_efermi = y_aver_efermi / npoints_efermi;
                x2_aver_efermi = x2_aver_efermi / npoints_efermi;
                xy_aver_efermi = xy_aver_efermi / npoints_efermi;

                *Te = (xy_aver_efermi - x_aver_efermi * y_aver_efermi) / (x2_aver_efermi - x_aver_efermi * x_aver_efermi);
                *efermi = y_aver_efermi - *Te * x_aver_efermi;

                delete [] x_efermi;
                delete [] y_efermi;
	return 0;
}

nablaclass::nablaclass()
{
	this->memoryexists = 0;
	this->nkpts = 0;
	this->nb_tot = 0;
	this->vkpt = NULL;
	this->wtkpt = NULL;
	this->celtot = NULL;
	this->fertot = NULL;
	this->nabij = NULL;
	this->efermi = 0;
	this->Te = 0;
}

int nablaclass::memory_create(int nkpts, int nb_tot)
{
	int kcounter, bandcounter1, bandcounter2;

	if ((this->memoryexists != 0) && ((this->nkpts == nkpts) && (this->nb_tot == nb_tot)))
	{
		return 0;
	}
	else if ((this->memoryexists != 0) && ((this->nkpts != nkpts) || (this->nb_tot != nb_tot)))
	{
		this->memory_remove();
	}
	
	this->nkpts = nkpts;
	this->nb_tot = nb_tot;

	this->wtkpt = new double [nkpts];
	
	this->vkpt = new double* [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		vkpt[kcounter] = new double [3];
	}
	

	this->celtot = new complex<double>* [nkpts];
	this->fertot = new double* [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		this->celtot[kcounter] = new complex<double> [nb_tot];
		this->fertot[kcounter] = new double [nb_tot];
	}

	this->nabij = new complex<float>*** [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		this->nabij[kcounter] = new complex<float>** [nb_tot];
		for (bandcounter1 = 0; bandcounter1 < nb_tot; bandcounter1++)
		{
			this->nabij[kcounter][bandcounter1] = new complex<float>*[nb_tot];
			for (bandcounter2 = 0; bandcounter2 < nb_tot; bandcounter2++)
			{
				this->nabij[kcounter][bandcounter1][bandcounter2] = new complex<float> [3];
			}
		}
	}
	
	this->memoryexists = 1;

	return 0;	
}

int nablaclass::input_binary(const char* inputfilename)
{
	ifstream inputfile;

	long int idum;
	
	double rdum;
	double rispin;
	double rnkpts;
	double rnb_tot;
	latt LATT;

	int kcounter;
	int bandcounter1, bandcounter2;
	int a;

	double convert;


	double *y_efermi;       //used to calculate the Fermi energy; values of energy
        double *x_efermi;       //used to calculate the Fermi energy; ln((2 - weight) / weight)
        int npoints_efermi = 0; //the number of points for the Fermi energy calculation; weights must be between 1.5 and 0.5
        int start_efermi = 0;   //the number of the first point with the weight less than 1.5
        int end_efermi = 0; 	//the number of the last point with the weight larger than 0.5
	double upper_boundary_efermi = 0.9;
	double lower_boundary_efermi = 0.1;
	int counter_efermi;
	double xy_aver_efermi, x_aver_efermi, y_aver_efermi, x2_aver_efermi;

	inputfile.open(inputfilename, ios_base::in|ios_base::binary);
	if (! inputfile)
	{
		return 1;
	}
	
	inputfile.seekg(0);
	inputfile.read((char*)&rdum, sizeof(double));
	inputfile.read((char*)&rispin, sizeof(double));
	inputfile.read((char*)&rnkpts, sizeof(double));
	inputfile.read((char*)&rnb_tot, sizeof(double));
	inputfile.read((char*)&(this->latt_ini), sizeof(latt));
	idum = (int)rdum;
	
	this->memory_create((int)rnkpts, (int)rnb_tot);
	
	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		inputfile.seekg(idum * (1 + 4 * kcounter));
		inputfile.read((char*)this->vkpt[kcounter], 3 * sizeof(double));
		inputfile.read((char*)(this->wtkpt + kcounter), sizeof(double));
		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			inputfile.read((char*)(this->celtot[kcounter] + bandcounter1), sizeof(complex<double>));
			inputfile.read((char*)(this->fertot[kcounter] + bandcounter1), sizeof(double));
		}

		for (a = 0; a < 3; a++)
		{
			inputfile.seekg(idum * (1 + 4 * kcounter + 1 + a));
			for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
			{
				for (bandcounter2 = 0; bandcounter2 < this->nb_tot; bandcounter2++)
				{
					inputfile.read((char*)&convert, sizeof(double));
					(this->nabij[kcounter][bandcounter1][bandcounter2][a]).real((float)convert);
					inputfile.read((char*)&convert, sizeof(double));
					(this->nabij[kcounter][bandcounter1][bandcounter2][a]).imag((float)convert);
				}
			}
		}	
	}
	
	start_efermi = 0;
	while (this->fertot[0][start_efermi] > upper_boundary_efermi)
	{
		start_efermi++;
	}
	end_efermi = start_efermi;
	while (this->fertot[0][end_efermi] > lower_boundary_efermi)
	{
		end_efermi++;
	}
	npoints_efermi = end_efermi - start_efermi;

	if (npoints_efermi < 2)
	{
//		if (rank == 0)
//			cerr << programname << ": WARNING: There are " << npoints_efermi << " bands with fermi-weights between " << upper_boundary_efermi << " and " << lower_boundary_efermi << endl;

		if (start_efermi - 4 < 0)
			start_efermi = 0;
		else
			start_efermi = start_efermi - 4;

		if (end_efermi + 4 > (this->nb_tot - 1))
			end_efermi = this->nb_tot;
		else
			end_efermi = end_efermi + 4;

		npoints_efermi = end_efermi - start_efermi;

		if (npoints_efermi < 2)
 		{
//			if (rank == 0)
//				cerr << programname << ": ERROR: Only one band present, impossible to calculate chemical potential and temperature of electrons ; program terminated" << endl;
//                                MPI_Finalize();
//                                exit(1);
			return 0;
                }
		else
		{
//			if (rank == 0)
//				cerr << programname << ": WARNING: Bands with fermi-weights between " << FERTOT[0][0][start_efermi] << " and " << FERTOT[0][0][end_efermi-1] << " will be used instead" << endl;
		}
	}

	x_efermi = new double [npoints_efermi];
	y_efermi = new double [npoints_efermi];

	for (counter_efermi = 0; counter_efermi < npoints_efermi; counter_efermi++)
	{
		x_efermi[counter_efermi] = log ((1 - this->fertot[0][start_efermi + counter_efermi]) / this->fertot[0][start_efermi + counter_efermi]);
		y_efermi[counter_efermi] = real(this->celtot[0][start_efermi + counter_efermi]);
                }

                x_aver_efermi = 0;
                y_aver_efermi = 0;
                x2_aver_efermi = 0;
                xy_aver_efermi = 0;
                for (counter_efermi = 0; counter_efermi < npoints_efermi; counter_efermi++)
                {
                        x_aver_efermi = x_aver_efermi + x_efermi[counter_efermi];
                        y_aver_efermi = y_aver_efermi + y_efermi[counter_efermi];
                        x2_aver_efermi = x2_aver_efermi + x_efermi[counter_efermi] * x_efermi[counter_efermi];
                        xy_aver_efermi = xy_aver_efermi + x_efermi[counter_efermi] * y_efermi[counter_efermi];
                }
                x_aver_efermi = x_aver_efermi / npoints_efermi;
                y_aver_efermi = y_aver_efermi / npoints_efermi;
                x2_aver_efermi = x2_aver_efermi / npoints_efermi;
                xy_aver_efermi = xy_aver_efermi / npoints_efermi;

                this->Te = (xy_aver_efermi - x_aver_efermi * y_aver_efermi) / (x2_aver_efermi - x_aver_efermi * x_aver_efermi);
                this->efermi = y_aver_efermi - Te * x_aver_efermi;

                delete [] x_efermi;
                delete [] y_efermi;



	inputfile.close();

	return 0;
}

int nablaclass::output_form(const char* outputfilename, const char* header)
{
	ofstream outputfile;

	int kcounter;
	int a;
	int bandcounter1, bandcounter2;

	outputfile.open(outputfilename, ios_base::out);
	
	outputfile << header << endl << endl;
	outputfile << "number of k-points in the irreducible wedge of the Brillouin zone: " << this->nkpts << endl;
	outputfile << "number of bands: " << this->nb_tot << endl;
	latticeprintfile(this->latt_ini, &outputfile);

	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		outputfile << endl;
		outputfile << "k-point number " << kcounter << ":"<< endl << endl;
		outputfile << "k-point coordinates (in reciprocal vector units: B0, B1, B2): " << endl << "         ";
		for (a = 0; a < 3; a++)
		{
			outputfile << this->vkpt[kcounter][a] << " " ;
		}
		outputfile << endl;
		outputfile << "weight of the k-point: " << this->wtkpt[kcounter] << endl << endl;
		
		outputfile << "Band number     Energy, eV     fermi-weight" << endl;

		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{	
			outputfile << bandcounter1 << "        " << real(this->celtot[kcounter][bandcounter1]) << "             " << this->fertot[kcounter][bandcounter1] << endl;
		}
		outputfile << endl;
	
		outputfile << "Nabla matrix elements: " << endl;
		outputfile << "Band 1       Band 2                          matrix elements (x, y, z), A^-1" << endl;
		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			for (bandcounter2 = 0; bandcounter2 < this->nb_tot; bandcounter2++)
			{
				outputfile << bandcounter1 << "             " << bandcounter2;
				for (a = 0; a < 3; a++)
				{
					outputfile << "           " << this->nabij[kcounter][bandcounter1][bandcounter2][a];
				}				
				outputfile << endl;
			}
			outputfile << endl;
		}
		outputfile << endl << endl;
		
	}


	outputfile.close();

	return 0;
}

int nablaclass::memory_remove()
{
	int kcounter, bandcounter1, bandcounter2;
	if (this->memoryexists != 0)
	{
		delete [] this->wtkpt;
		
		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			delete [] this->vkpt[kcounter];
		}
		delete [] this->vkpt;

		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			delete [] this->celtot[kcounter];
			delete [] this->fertot[kcounter];
		}
		delete [] this->celtot;
		delete [] this->fertot;

		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
			{
				for (bandcounter2 = 0; bandcounter2 < this->nb_tot; bandcounter2++)
				{
					delete [] this->nabij[kcounter][bandcounter1][bandcounter2];
				}
				delete [] this->nabij[kcounter][bandcounter1];
			}
			delete [] this->nabij[kcounter];
		}
		delete [] this->nabij;

		this->nkpts = 0;
		this->nb_tot = 0;
		latticeclean(&(this->latt_ini));
		this->memoryexists = 0;
		
	}
	return 0;
}

nablaclass::~nablaclass()
{
}

nablaaverclass::nablaaverclass()
{
	this->memoryexists = 0;
	this->nkpts = 0;
	this->nb_tot = 0;
	this->vkpt = NULL;
	this->wtkpt = NULL;
	this->celtot = NULL;
	this->fertot = NULL;
	this->nabij = NULL;
	this->efermi = 0;
	this->Te = 0;
	latticeclean(&(this->latt_ini));
}

int nablaaverclass::memory_create (int nkpts, int nb_tot)
{
	int kcounter, bandcounter1, bandcounter2;

	if (this->memoryexists == 1)
	{
		this->memory_remove();
	}
	
	this->nkpts = nkpts;
	this->nb_tot = nb_tot;

	this->wtkpt = new double [nkpts];
	
	this->vkpt = new double* [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		vkpt[kcounter] = new double [3];
	}
	

	this->celtot = new double* [nkpts];
	this->fertot = new double* [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		this->celtot[kcounter] = new double [nb_tot];
		this->fertot[kcounter] = new double [nb_tot];
	}

	this->nabij = new double** [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		this->nabij[kcounter] = new double* [nb_tot];
		for (bandcounter1 = 0; bandcounter1 < nb_tot; bandcounter1++)
		{
			this->nabij[kcounter][bandcounter1] = new double[nb_tot];
		}
	}
	
	this->memoryexists = 1;

	return 0;	
}

int nablaaverclass::space_average(nablaclass* input)
{
	int kcounter, coordcounter, bandcounter1, bandcounter2;

	this->memory_create(input->nkpts, input->nb_tot);
	this->nkpts = input->nkpts;
	this->nb_tot = input->nb_tot;
	latticecopy(input->latt_ini, &this->latt_ini);
	this->Te = input->Te;
	this->efermi = input->efermi;
	
	for (kcounter = 0; kcounter < input->nkpts; kcounter++)
	{
		this->wtkpt[kcounter] = input->wtkpt[kcounter];
		for (coordcounter = 0; coordcounter < 3; coordcounter++)
		{
			this->vkpt[kcounter][coordcounter] = input->vkpt[kcounter][coordcounter];
		}

		for (bandcounter1 = 0; bandcounter1 < input->nb_tot; bandcounter1++)
		{
			this->celtot[kcounter][bandcounter1] = real(input->celtot[kcounter][bandcounter1]);
			this->fertot[kcounter][bandcounter1] = input->fertot[kcounter][bandcounter1];
		}
		
		for (bandcounter1 = 0; bandcounter1 < input->nb_tot; bandcounter1++)
		{
			for (bandcounter2 = 0; bandcounter2 < input->nb_tot; bandcounter2++)
			{
				this->nabij[kcounter][bandcounter1][bandcounter2] = 0;
				for (coordcounter = 0; coordcounter < 3; coordcounter++)
				{
					this->nabij[kcounter][bandcounter1][bandcounter2] = this->nabij[kcounter][bandcounter1][bandcounter2] + abs(input->nabij[kcounter][bandcounter1][bandcounter2][coordcounter]) * abs(input->nabij[kcounter][bandcounter1][bandcounter2][coordcounter]);
				}
				this->nabij[kcounter][bandcounter1][bandcounter2] = this->nabij[kcounter][bandcounter1][bandcounter2] / (double)3;
			}
		}
	}
}

int nablaaverclass::align(int aligntype)
{
	double alignenergy = 0;
	int kcounter, bandcounter;

	if (aligntype == 0)
	{
		return 0;
	}
	else if (aligntype == 1)
	{
		alignenergy = this->efermi;
	}
	else if (aligntype == 2)
	{
		alignenergy = this->celtot[0][0];
	}

	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		for (bandcounter = 0; bandcounter < this->nb_tot; bandcounter++)
		{
			this->celtot[kcounter][bandcounter] = this->celtot[kcounter][bandcounter] - alignenergy;
		}
		this->efermi = this->efermi - alignenergy;
	}


	return 0;
}

int nablaaverclass::input_binary(const char* inputfilename)
{
	ifstream inputfile;

	long int idum;
	
	double rdum;
	double rispin;
	double rnkpts;
	double rnb_tot;
	latt LATT;

	int kcounter;
	int bandcounter1, bandcounter2;
	int coordcounter;
	int a;

	double convert;


	inputfile.open(inputfilename, ios_base::in|ios_base::binary);
	if (! inputfile)
	{
		return 1;
	}
	
	inputfile.seekg(0);
	inputfile.read((char*)&rdum, sizeof(double));
	inputfile.read((char*)&rnkpts, sizeof(double));
	inputfile.read((char*)&rnb_tot, sizeof(double));
	inputfile.read((char*)&this->efermi, sizeof(double));
	inputfile.read((char*)&this->Te, sizeof(double));
	inputfile.read((char*)&(this->latt_ini), sizeof(latt));
	idum = (int)rdum;
	
	this->memory_create((int)rnkpts, (int)rnb_tot);
	
	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		inputfile.seekg(idum * (1 + 2 * kcounter));
		inputfile.read((char*)this->vkpt[kcounter], 3 * sizeof(double));

		inputfile.read((char*)(this->wtkpt + kcounter), sizeof(double));

		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			inputfile.read((char*)(this->celtot[kcounter] + bandcounter1), sizeof(double));
			inputfile.read((char*)(this->fertot[kcounter] + bandcounter1), sizeof(double));
		}

		inputfile.seekg(idum * (1 + 2 * kcounter + 1));
		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			for (bandcounter2 = 0; bandcounter2 < this->nb_tot; bandcounter2++)
			{
				inputfile.read((char*)(this->nabij[kcounter][bandcounter1] + bandcounter2), sizeof(double));
			}
		}
	}
	

	inputfile.close();

	return 0;
}

int nablaaverclass::output_binary(const char* outputfilename)
{
	int kcounter, bandcounter1, bandcounter2;
	ofstream outputfile;
	long int idum;
	int coordcounter;

	double convert;

	outputfile.open(outputfilename, ios_base::out|ios_base::binary);
	
	idum = 5 * sizeof(double) + sizeof(latt);
	if ((4 + 2 * this->nb_tot) * sizeof(double) > idum)
		idum = (4 + 2 * this->nb_tot) * sizeof(double);
	if (this->nb_tot * this->nb_tot * sizeof(double)> idum)
		idum = this->nb_tot * this->nb_tot * sizeof(double);

	outputfile.seekp(0);
	convert = idum;
	outputfile.write((char*)&convert, sizeof(double));
	convert = this->nkpts;
	outputfile.write((char*)&convert, sizeof(double));
	convert = this->nb_tot;
	outputfile.write((char*)&convert, sizeof(double));
	outputfile.write((char*)&this->efermi, sizeof(double));
	outputfile.write((char*)&this->Te, sizeof(double));
	outputfile.write((char*)&this->latt_ini, sizeof(latt));

	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		outputfile.seekp((1 + 2 * kcounter) * idum);
		for (coordcounter = 0; coordcounter < 3; coordcounter++)
		{
			outputfile.write((char*)(this->vkpt[kcounter] + coordcounter), sizeof(double));
		}
		outputfile.write((char*)(this->wtkpt + kcounter), sizeof(double));
		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			outputfile.write((char*)(this->celtot[kcounter] + bandcounter1),sizeof(double));
			outputfile.write((char*)(this->fertot[kcounter] + bandcounter1),sizeof(double));
		}

		outputfile.seekp((1 + 2 * kcounter + 1) * idum);
		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			for (bandcounter2 = 0; bandcounter2 < this->nb_tot; bandcounter2++)
			{
				outputfile.write((char*)(this->nabij[kcounter][bandcounter1] + bandcounter2),sizeof(double));
			}
		}
	}

	outputfile.close();
}

int nablaaverclass::output_form(const char* outputfilename, const char* header)
{
	ofstream outputfile;

	int kcounter;
	int a;
	int bandcounter1, bandcounter2;

		

	outputfile.open(outputfilename, ios_base::out);
	
	outputfile << header << endl << endl;
	outputfile << "number of k-points in the irreducible wedge of the Brillouin zone: " << this->nkpts << endl;
	outputfile << "number of bands/energy values: " << this->nb_tot << endl;
	outputfile << "Fermi energy: " << this->efermi << " eV" << endl;
	outputfile << "Te: " << this->Te * 11604.505 << " K" << endl;
	latticeprintfile(this->latt_ini, &outputfile);

	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		outputfile << endl;
		outputfile << "k-point number " << kcounter << ":"<< endl << endl;
		outputfile << "k-point coordinates (in reciprocal vector units: B0, B1, B2): " << endl << "         ";
		for (a = 0; a < 3; a++)
		{
			outputfile << this->vkpt[kcounter][a] << " " ;
		}
		outputfile << endl;
		outputfile << "weight of the k-point: " << this->wtkpt[kcounter] << endl << endl;
		
		outputfile << "Band/energy point number     Energy, eV     fermi-weight" << endl;

		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{	
			outputfile << bandcounter1 << "                     " << this->celtot[kcounter][bandcounter1] << "             " << this->fertot[kcounter][bandcounter1] << endl;
		}
		outputfile << endl;
	
		outputfile << "Nabla matrix elements: " << endl;
		outputfile << "Band/energy point 1       Band/energy point 2                          quantity" << endl;
		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			for (bandcounter2 = 0; bandcounter2 < this->nb_tot; bandcounter2++)
			{
				outputfile << bandcounter1 << "                          " << bandcounter2;
				outputfile << "                        " << this->nabij[kcounter][bandcounter1][bandcounter2];
				outputfile << endl;
			}
			outputfile << endl;
		}
		outputfile << endl << endl;
		
	}


	outputfile.close();

	return 0;
}

int nablaaverclass::output_section_number_Origin(int kpointnumber, int bandnumber, const char* outputname)
{
	int bandcounter;

	ofstream outputfile;
	
	outputfile.open(outputname, ios_base::out);

	outputfile << "\\g(e)\\-(2), eV\t|nabla(\\g(e)\\-(1),\\g(e)\\-(2))|\\+(2), A\\+(-2)\tf(\\g(e))\tnumber\\-(2)" << endl;
	outputfile << "\t\t\t" << endl;
	outputfile << "\tk-point" << kpointnumber << " \\g(e)\\-(1)=" << this->celtot[kpointnumber][bandnumber] << "eV number\\-(1)=" << bandnumber;
	outputfile << "\tk-point" << kpointnumber << " \\g(e)\\-(1)=" << this->celtot[kpointnumber][bandnumber] << "eV number\\-(1)=" << bandnumber;
	outputfile << "\tk-point" << kpointnumber << " \\g(e)\\-(1)=" << this->celtot[kpointnumber][bandnumber] << "eV number\\-(1)=" << bandnumber;
	outputfile << endl;

	for (bandcounter = 0; bandcounter < this->nb_tot; bandcounter++)
	{
		outputfile << this->celtot[kpointnumber][bandcounter] << "\t" << this->nabij[kpointnumber][bandnumber][bandcounter] << "\t" << this->fertot[kpointnumber][bandcounter] << "\t" << bandcounter << endl;
	}

	outputfile.close();

	return 0;
}

int nablaaverclass::uniform_mesh()
{
	int kcounter, bandcounter;
	double ediff = 0;

	ediff = this->celtot[0][1] - this->celtot[0][0];
	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		for (bandcounter = 0; bandcounter < this->nb_tot - 1; bandcounter++)
		{
			if (abs((this->celtot[kcounter][bandcounter+1] - this->celtot[kcounter][bandcounter]) - ediff) > 1e-10)
				return 0;
		}
	}
	
	return 1;
}

double nablaaverclass::get_omega_real(double omegarequired)
{
	double ediff = 0;
	int banddiff = 0;

	ediff = this->celtot[0][1] - this->celtot[0][0];
	banddiff = (int)round(omegarequired / ediff);

	return (banddiff * ediff);
}

int nablaaverclass::output_section_omega_Origin(int kpointnumber, double omegarequired, const char* outputname)
{
	int kcounter, bandcounter;
	double ediff = 0;
	double omegareal = 0;

	int banddiff = 0;

	ofstream outputfile;
	
	
	ediff = this->celtot[0][1] - this->celtot[0][0];
	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		for (bandcounter = 0; bandcounter < this->nb_tot - 1; bandcounter++)
		{
			if (abs((this->celtot[kcounter][bandcounter+1] - this->celtot[kcounter][bandcounter]) - ediff) > 1e-10)
				return 1;
		}
	}

	banddiff = (int)round(omegarequired / ediff);
	omegareal = banddiff * ediff;
	

	outputfile.open(outputname, ios_base::out);

	outputfile << "\\g(e), eV\t|nabla(\\g(e),\\g(e)+h\\g(w))|\\+(2), A\\+(-2)" << endl;
	outputfile << "\t" << endl;
	outputfile << "\tk-point" << kpointnumber << " \\g(w)=" << omegareal << "eV";
	outputfile << endl;

	if (banddiff >= 0)
	{
		bandcounter = 0;
		while ((bandcounter + banddiff) < this->nb_tot)
		{
			outputfile << this->celtot[kpointnumber][bandcounter] << "\t" << this->nabij[kpointnumber][bandcounter][bandcounter+banddiff] << endl;
			bandcounter++;
		}
	}
	else if (banddiff < 0)
	{
		bandcounter = (-1) * banddiff;
		while ((bandcounter + banddiff) < this->nb_tot)
		{
			outputfile << this->celtot[kpointnumber][bandcounter] << "\t" << this->nabij[kpointnumber][bandcounter][bandcounter+banddiff] << endl;
			bandcounter++;
		}
	}

	outputfile.close();

	return 0;
}

int nablaaverclass::memory_remove()
{
	int kcounter, bandcounter1, bandcounter2;
	if (this->memoryexists == 1)
	{
		delete [] this->wtkpt;
		
		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			delete [] this->vkpt[kcounter];
		}
		delete [] this->vkpt;

		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			delete [] this->celtot[kcounter];
			delete [] this->fertot[kcounter];
		}
		delete [] this->celtot;
		delete [] this->fertot;

		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
			{
				delete [] this->nabij[kcounter][bandcounter1];
			}
			delete [] this->nabij[kcounter];
		}
		delete [] this->nabij;

		this->nkpts = 0;
		this->nb_tot = 0;
		latticeclean(&(this->latt_ini));
		this->memoryexists = 0;
		
	}
	return 0;
}

nablaaverclass::~nablaaverclass()
{
}

bandsystemclass::bandsystemclass()
{
	this->memoryexists = 0;
	this->nkpts = 0;
	this->nb_tot = 0;
	this->vkpt = NULL;
	this->wtkpt = NULL;
	this->celtot = NULL;
	this->fertot = NULL;
	this->efermi = 0;
	this->Te = 0;
	latticeclean(&(this->latt_ini));
}


int bandsystemclass::memory_create(int nkpts, int nb_tot)
{
	int kcounter;

	if (this->memoryexists == 1)
	{
		this->memory_remove();
	}
	
	this->nkpts = nkpts;
	this->nb_tot = nb_tot;

	this->wtkpt = new double [nkpts];
	
	this->vkpt = new double* [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		vkpt[kcounter] = new double [3];
	}
	

	this->celtot = new double* [nkpts];
	this->fertot = new double* [nkpts];
	for (kcounter = 0; kcounter < nkpts; kcounter++)
	{
		this->celtot[kcounter] = new double [nb_tot];
		this->fertot[kcounter] = new double [nb_tot];
	}

	
	this->memoryexists = 1;

	return 0;
}

int bandsystemclass::input_wavecar(const char* wavecarname)
{
	ifstream wavecar;
	double rdum = 0;
	double rispin = 0;
	double rtag = 0;
	double rnkpts = 0;
	double rnb_tot = 0;
	double enmax = 0;
	double rnpl = 0;

	long int idum = 0;
	int nkpts = 0;
	int nb_tot = 0;

	int kcounter = 0;
	int bandcounter = 0;
	int coordcounter = 0;

	complex<double> convert (0,0);

	wavecar.open(wavecarname, ios_base::in|ios_base::binary);
	if (! wavecar)
		return 1;

	wavecar.read((char*)&rdum, sizeof(double));
	wavecar.read((char*)&rispin, sizeof(double));
	wavecar.read((char*)&rtag, sizeof(double));
	idum = (int)(floor(rdum + 0.5));

	wavecar.seekg(idum);
	wavecar.read((char*)&rnkpts, sizeof(double));
	wavecar.read((char*)&rnb_tot, sizeof(double));
	wavecar.read((char*)&enmax, sizeof(double));


	lattice_set_nan(&(this->latt_ini));
	for (coordcounter = 0; coordcounter < 3; coordcounter++)
	{
		wavecar.read((char*)(this->latt_ini.A[coordcounter]), 3*sizeof(double));
	}


	nkpts = (int)(floor(rnkpts + 0.5)); nb_tot = (int)(floor(rnb_tot + 0.5));

	this->memory_create(nkpts, nb_tot);

	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		wavecar.seekg((2 + kcounter * (this->nb_tot + 1)) * idum);
		wavecar.read((char*)&rnpl, sizeof(double));
		wavecar.read((char*)this->vkpt[kcounter],3*sizeof(double));
		this->wtkpt[kcounter] = NAN;		

		for (bandcounter = 0; bandcounter < this->nb_tot; bandcounter++)
		{
			wavecar.read((char*)&convert, sizeof(complex<double>));
			this->celtot[kcounter][bandcounter] = real(convert);
			wavecar.read((char*)(this->fertot[kcounter]+bandcounter), sizeof(double));
		}
	}


	get_efermi_Te(this->nb_tot, this->celtot[0], this->fertot[0], &(this->efermi), &(this->Te));

	return 0;
}

int bandsystemclass::input_nabla(const char* nablaname)
{

	ifstream inputfile;

	long int idum;
	
	double rdum;
	double rispin;
	double rnkpts;
	double rnb_tot;
	latt LATT;

	int nkpts = 0;
	int nb_tot = 0;

	int kcounter;
	int bandcounter;
	int a;

	complex<double> convert(0,0);


	inputfile.open(nablaname, ios_base::in|ios_base::binary);
	if (! inputfile)
	{
		return 1;
	}
	
	inputfile.seekg(0);
	inputfile.read((char*)&rdum, sizeof(double));
	inputfile.read((char*)&rispin, sizeof(double));
	inputfile.read((char*)&rnkpts, sizeof(double));
	inputfile.read((char*)&rnb_tot, sizeof(double));
	inputfile.read((char*)&(this->latt_ini), sizeof(latt));
	idum = (int)rdum;
	
	nkpts = (int)(floor(rnkpts + 0.5)); nb_tot = (int)(floor(rnb_tot + 0.5));
	this->memory_create(nkpts, nb_tot);
	
	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		inputfile.seekg(idum * (1 + 4 * kcounter));
		inputfile.read((char*)this->vkpt[kcounter], 3 * sizeof(double));
		inputfile.read((char*)(this->wtkpt + kcounter), sizeof(double));
		for (bandcounter = 0; bandcounter < this->nb_tot; bandcounter++)
		{
			inputfile.read((char*)(&convert), sizeof(complex<double>));
			this->celtot[kcounter][bandcounter] = real(convert);
			inputfile.read((char*)(this->fertot[kcounter] + bandcounter), sizeof(double));
		}

	}
	
	get_efermi_Te(this->nb_tot, this->celtot[0], this->fertot[0], &(this->efermi), &(this->Te));

	inputfile.close();

	return 0;

}

int bandsystemclass::input_spaceaver(const char* nablaaverfilename)
{
	ifstream inputfile;

	long int idum;
	
	double rdum;
	double rispin;
	double rnkpts;
	double rnb_tot;
	latt LATT;

	int nkpts = 0;
	int nb_tot = 0;

	int kcounter;
	int bandcounter1;
	int a;

	double convert;


	inputfile.open(nablaaverfilename, ios_base::in|ios_base::binary);
	if (! inputfile)
	{
		return 1;
	}
	
	inputfile.seekg(0);
	inputfile.read((char*)&rdum, sizeof(double));
	inputfile.read((char*)&rnkpts, sizeof(double));
	inputfile.read((char*)&rnb_tot, sizeof(double));
	inputfile.read((char*)&this->efermi, sizeof(double));
	inputfile.read((char*)&this->Te, sizeof(double));
	inputfile.read((char*)&(this->latt_ini), sizeof(latt));
	idum = (int)rdum;
	
	nkpts = (int)(floor(rnkpts + 0.5)); nb_tot = (int)(floor(rnb_tot + 0.5));
	this->memory_create(nkpts, nb_tot);
	
	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		inputfile.seekg(idum * (1 + 2 * kcounter));
		inputfile.read((char*)this->vkpt[kcounter], 3 * sizeof(double));

		inputfile.read((char*)(this->wtkpt + kcounter), sizeof(double));

		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{
			inputfile.read((char*)(this->celtot[kcounter] + bandcounter1), sizeof(double));
			inputfile.read((char*)(this->fertot[kcounter] + bandcounter1), sizeof(double));
		}

	}
	

	inputfile.close();

	return 0;
}

int bandsystemclass::align(int aligntype)
{

	double alignenergy = 0;
	int kcounter, bandcounter;

	if (aligntype == 0)
	{
		return 0;
	}
	else if (aligntype == 1)
	{
		alignenergy = this->efermi;
	}
	else if (aligntype == 2)
	{
		alignenergy = this->celtot[0][0];
	}
	else if (aligntype == 3)
	{
		alignenergy = this->celtot[0][0];
		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			if (alignenergy > this->celtot[kcounter][0])
				alignenergy = this->celtot[kcounter][0];
		}
	}

	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		for (bandcounter = 0; bandcounter < this->nb_tot; bandcounter++)
		{
			this->celtot[kcounter][bandcounter] = this->celtot[kcounter][bandcounter] - alignenergy;
		}
		this->efermi = this->efermi - alignenergy;
	}

	return 0;
}

int bandsystemclass::output_form(const char* outputfilename, const char* header)
{
	ofstream outputfile;

	int kcounter;
	int a;
	int bandcounter1, bandcounter2;

	outputfile.open(outputfilename, ios_base::out);
	
	outputfile << header << endl << endl;
	outputfile << "number of k-points in the irreducible wedge of the Brillouin zone: " << this->nkpts << endl;
	outputfile << "number of bands/energy values: " << this->nb_tot << endl;
	outputfile << "Fermi energy: " << this->efermi << " eV" << endl;
	outputfile << "Te: " << this->Te * 11604.505 << " K" << endl;
	latticeprintfile(this->latt_ini, &outputfile);

	for (kcounter = 0; kcounter < this->nkpts; kcounter++)
	{
		outputfile << endl;
		outputfile << "k-point number " << kcounter << ":"<< endl << endl;
		outputfile << "k-point coordinates (in reciprocal vector units: B0, B1, B2): " << endl << "         ";
		for (a = 0; a < 3; a++)
		{
			outputfile << this->vkpt[kcounter][a] << " " ;
		}
		outputfile << endl;
		outputfile << "weight of the k-point: " << this->wtkpt[kcounter] << endl << endl;
		
		outputfile << "Band/energy point number     Energy, eV     fermi-weight" << endl;

		for (bandcounter1 = 0; bandcounter1 < this->nb_tot; bandcounter1++)
		{	
			outputfile << bandcounter1 << "                     " << this->celtot[kcounter][bandcounter1] << "             " << this->fertot[kcounter][bandcounter1] << endl;
		}
	
		outputfile << endl << endl;
		
	}


	outputfile.close();

	return 0;
}

int bandsystemclass::make_dos_kpoint(conductivity* dos, int kpoint)
{
	int** enlower = NULL;
	int** enupper = NULL;
	int deltaecounter = 0;
	int energycounter = 0;
	int bandcounter = 0;

	enlower = new int* [dos->deltaenum];
	enupper = new int* [dos->deltaenum];
	for (deltaecounter = 0; deltaecounter < dos->deltaenum; deltaecounter++)
	{
		enlower[deltaecounter] = new int [dos->omeganum];
		enupper[deltaecounter] = new int [dos->omeganum];
	}


	for (deltaecounter = 0; deltaecounter < dos->deltaenum; deltaecounter++)
	{
		for (energycounter = 0; energycounter < dos->omeganum; energycounter++)
		{
			bandcounter = 0;
			while ((this->celtot[kpoint][bandcounter] - dos->omega[energycounter] < (-1) * dos->deltae[deltaecounter] * 10) && (bandcounter < this->nb_tot - 1))
				bandcounter++;


                        if (bandcounter > 0)
                                enlower[deltaecounter][energycounter] = bandcounter - 1;
                        else if (bandcounter == 0)
                                enlower[deltaecounter][energycounter] = 0;

                        while ((this->celtot[kpoint][bandcounter] - dos->omega[energycounter] < dos->deltae[deltaecounter] * 10) && (bandcounter < this->nb_tot - 1))
                                bandcounter++;

                        enupper[deltaecounter][energycounter] = bandcounter;
		}
	}

	for (deltaecounter = 0; deltaecounter < dos->deltaenum; deltaecounter++)
	{
		for (energycounter = 0; energycounter < dos->omeganum; energycounter++)
		{	
			dos->sigma[deltaecounter][energycounter] = 0;
			for (bandcounter = enlower[deltaecounter][energycounter]; bandcounter < enupper[deltaecounter][energycounter]; bandcounter++)
			{
				dos->sigma[deltaecounter][energycounter] += deltagaussian(dos->omega[energycounter] - this->celtot[kpoint][bandcounter],dos->deltae[deltaecounter]);
			}
		}
	}

	for (deltaecounter = 0; deltaecounter < dos->deltaenum; deltaecounter++)
	{
		delete [] enlower[deltaecounter];
		delete [] enupper[deltaecounter];
	}
	delete [] enlower;
	delete [] enupper;

	return 0;
}

int bandsystemclass::memory_remove()
{
	int kcounter;
	if (this->memoryexists == 1)
	{
		delete [] this->wtkpt;
		
		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			delete [] this->vkpt[kcounter];
		}
		delete [] this->vkpt;

		for (kcounter = 0; kcounter < this->nkpts; kcounter++)
		{
			delete [] this->celtot[kcounter];
			delete [] this->fertot[kcounter];
		}
		delete [] this->celtot;
		delete [] this->fertot;

		this->nkpts = 0;
		this->nb_tot = 0;
		latticeclean(&(this->latt_ini));
		this->memoryexists = 0;
		
	}
	return 0;
}

bandsystemclass::~bandsystemclass()
{

}

#endif
