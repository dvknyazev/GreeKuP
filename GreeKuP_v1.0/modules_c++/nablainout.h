#ifndef NABLAINOUT_H
#define NABLAINOUT_H
#include "lattice.h"
#include "myinout.h"

class nablaclass
{
	public:
		int nkpts;			//total number of k-points
		int nb_tot;			//total number of bands
		double** vkpt;			//coordinates of k-points
		double* wtkpt;			//weights of k-points
		complex<double>** celtot;	//eigenvalues
		double** fertot;		//fermi-weights
		complex<float>**** nabij;	//matrix elements
		latt latt_ini;
		int memoryexists;

		double efermi;
		double Te;

		nablaclass();
		int memory_create(int, int);

		int input_binary(const char*);

		int output_form(const char*, const char*);

		int memory_remove();
		~nablaclass();
};

class nablaaverclass
{
	public:
		int nkpts;
		int nb_tot;
		double** vkpt;
		double* wtkpt;
		double** celtot;
		double** fertot;
		double*** nabij;
		latt latt_ini;
		double efermi;
		double Te;
		int memoryexists;

		nablaaverclass();
		int memory_create(int, int);
		int space_average(nablaclass*);
		int align (int);
		int input_binary (const char *);
		int output_binary (const char *);
		int output_form (const char*, const char*);
		int output_section_number_Origin(int, int, const char*);

		int uniform_mesh();
		double get_omega_real(double);
		int output_section_omega_Origin(int, double, const char*);
		int memory_remove();
		~nablaaverclass();
};

class bandsystemclass
{
	public:
		int nkpts;
		int nb_tot;
		double** vkpt;
		double* wtkpt;
		double** celtot;
		double** fertot;
		latt latt_ini;
		double efermi;
		double Te;
		int memoryexists;

		bandsystemclass();
		int memory_create(int, int);

		int input_wavecar (const char*);
		int input_nabla (const char*);	
		int input_spaceaver (const char*);	

		int align (int);

		int make_dos_kpoint(conductivity*, int);

		int output_form (const char*, const char*);

		int memory_remove();
		~bandsystemclass();
};
#endif
