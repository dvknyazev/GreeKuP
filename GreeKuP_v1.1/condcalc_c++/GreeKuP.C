//version GreeKuP_v1.1

#include <iostream>
#include <cmath>
using namespace std;
#include <fstream>
#include <complex>
#include <sstream>
#include <math.h>
#include <time.h>
#include "mpi.h"
#include "myinout.h"
#include "keys.h"
#include "omegainout.h"
#include "mymath.h"
#include <stdlib.h>

struct latt //structure that contains data about direct and reciprocal lattice (same as defined in vasp, lattice.inc)
{
	double SCALE;
	double A[3][3]; //array of vectors of direct lattice A[vector number][coordinate(x,y,z) number]
	double B[3][3]; //array of vectors of reciprocal lattice B[vector number][coordinate(x,y,z) number]
	double ANORM[3]; //array of vector lengths (direct lattice) ANORM[vector number]
	double BNORM[3]; //array of vector lengths (reciprocal lattice) BNORM[vector number]
	double OMEGA; //cell volume (in direct space)
};

//complex<float> conj(complex<float> z) //conjugates the complex number
//{
//        return complex<float> (z.real(), -z.imag());
//}


int printtimehms(double timeseconds, ofstream* timefile)
{
	int hours, minutes;
	double seconds;
	hours = (int)timeseconds / 3600;
	minutes = (int)(timeseconds - hours * 3600) / 60;
	seconds = timeseconds - hours * 3600 - minutes * 60;
	*timefile << hours << ":" << minutes << ":" << seconds << endl;
	return 0;
}

double memory_to_GB (long int memory)
{	
	double GBs;
	GBs = (double)memory / 1000000000;
	GBs = (double)((int)(GBs * 1000 + 0.5)) / 1000;
	return GBs;
}
	


int main(int argc, char* argv[])
{
	char* programname = "GreeKuP";
	
	int allkptswrite = 0;	//defines wheteher scalar conductivity for each kpoint should be printed (in binary and formatted filea
				//default is zero; may be changed by 'k' key

	double RDUM; //record length in bytes
	double RISPIN; //number of spins
	double RTAG;  //precision of WAVECAR file: 45200 - single
	double RNKPTS;  //number of k-points in the irreducable wedge of the BZ (IBZ)
	double RNB_TOT; //total number of bands
	double ENMAX; //energy cutoff
	double RNPL; //number of plane waves (for a current k-point)
	double nablareadrdum;//record length of the input file with matrix elements
	
	long int IDUM; //record length (same as RDUM, but integer)
	int ISPIN; //number of spins (same as RISPIN, but integer)
	int ITAG; //precision of WAVECAR file (same as RTAG, but integer)
	int NKPTS; //number of k-points (same as RNKPTS, but integer)
	int NB_TOT; //total number of bands (same as RNB_TOT, but integer)
	int deltaenum; //number of deltae points (deltae - width of gaussian distribution? used instead of delta-function)
	int omeganum; //number of omega points (omega - frequency)
	int condidum; //CONDUCTIVITYFILE record length in bytes
	long int nablawriteidum;	//NABLA file record length in bytes
	long int nablareadidum;//record length of the input file with matrix elements (same as nablareadrdum, but integer)
	
	int ISP; //spin counter
	int K,K1; //k-point counter
	int J,J1,J2; //band counter
	int J1proc,J2proc; //band counter, number of bands with the reference point different for each processor
	int IREC; //record counter
	int I; //plane wave counter
	int a,a1; //coordinate counter
	int b; //basis vectors counter
	int deltaecounter; //counter of deltae points
	int omegacounter; //counter of omega points

	double convert; //used for converting from other types to double (necessary while writing in binary file CONDUCTIVITYFILE)
	double convert1;
	double convert_dbl_to_flt;	//used to convert double number to float

	latt LATT_INI; //structure with the initial lattice (it is assumed that cell shape and volume does not change during molecular dynamics)

	ifstream wavecar; //WAVECAR file
	ifstream planewave; //PLANEWAVE file
	ofstream conductivityfile; //CONDUCTIVITYFILE file
	ofstream conductivityform; //CONDUCTIVITYFORM file
	ofstream progress; //progress file
	ofstream nablawrite;	//output file with matrix elements (OPTICFULL format)
	ifstream nablaread;	//input file with matrix elements (OPTICFULL format)
	ofstream testoutput;

//	conductivity realcond;

	int wavecarkeypresent = 0; 	//defines whether -wavecar key is present
					//-wavecar parameter defines the name of WAVECAR file
	string wavecarname;		//defines the name of WAVECAR file

	int planewavekeypresent = 0;	//defines whether -planewave key is present
					//-planewave parameter defines the name of PLANEWAVE file
	string planewavename;		//defines the name of PLANEWAVE file

	int nablawritepresent = 0;	//defines whether -nablawrite key is present and whether matrix elements should be written
					//-nablawrite parameter define the name of the output file for matrix elements
	string nablawritename;		//defines the name of NABLA file

	int nablareadpresent = 0;	//defines whether -nablaread key is present
					//-nablaread parameter define the name of the input file with matrix elements
					//if -nablaread key is present the existence and values of -wavecar and -planewave keys
					//are ignored
					//if -nablaread key is present, nablawritepresent is set to 0 and matrix elements
					//will not be written
	string nablareadname;		//defines the name of input file with matrix elements

	int kpathpresent = 0;		//defines whether -kpath key is present
					//-kpath parameter defines the name of CONDUCTIVITYFILE for each k-point
					//default is CONDUCTIVITYFILE (all k-points)
	stringstream kpointfilenum;	//used to convert file (k-point) number to string
	string kpointfilename;		//used to make name of file for a k-point
	string kpathname;
	int kstart = 0;			//defines the number of first k-point CONDUCTIVITY file (default is zero)

	int zerocorrection = 2;
	double fermidiff;		//used to calculate difference between fermi weights during Kubo-Greenwood calculation
					//zerocorrection value is taken into account

	omegaclass omegaobject;		//used to read data from OMEGA file - description in omegainout.h
	int omegapresent = 0;		//defines whether -omega key is present
	string omeganame;		//the name of OMEGA file
	double* omega;
	double* deltae;
	int omegainputresult = 0;		//the result of OMEGA file opening: 0 - correct, 1 - incorrect

	int onsagerfound = 0;		//determines whether -onsager key is present

	int recoulesfound = 0;	//determines whether -recoules key is present
					//onsager coefficients will be calculated according to the formulas from
					//V.Recoules and J.-P.Crocombette, Phys. Rev. B 72, 104202 (2005)
					//(\eps_k-\mu) wiil be used instead of ((\eps_k+\eps_k')/2-\mu)a

	int dividebyfound = 0;		//determines whether the -divideby=divideby key is present
	int divideby = 0;		//determines what the Kubo-Greenwood formula should be divided by:
					//0 - \hbar\omega; 1 - |e_i-e_j|
					//default is 0
	double dividebyvalue = 0;	//temporary variable to store either \hbar\omega of e_i-e_j value

	int resigmapresent = 0;		//determines whether -resigma key is present
	string resigmaname = "";	//the name of the output ReCONDUCTIVITY file

	int resigmaformpresent = 0;	//determines whether -resigmaform key is present
	string resigmaformname = "";	//the name of the output ReCONDUCTIVITYFORM file

	int progressfound = 0;		//determines whether -progress key is found
	string progressname = "";	//the name of the output progress file

	time_t time1, time2; //time1 - used for starting of time counting; time2 - used for finishing of time counting;
	time_t time3, time4; //time3 - used for starting of time counting; time4 - used for finishing of time counting;
	time_t time_beg, time_end;
	int timesectioncounter; //used for calculating of time sections elapsed
	int timinginterval = 1; //interval between successive reports about calculation (measured in minutes)
	
	int nb_proc_max = 1;	//determines maximum number of bands per one processor
				//nb_proc_max = bands / procnumber, if bands % procnumber == 0
				//nb_proc_max = bands / procnumber + 1, if bands % procnumber != 0

	complex<float>*** nablamatrelem_proc = NULL;

	int size, rank;
	int proccounter;	//process counter
	int proccounter1;
	int nablapos;		//position of an element in nablamatrelemlin array
	MPI_Status status;
	int mpibuffsize, mpibuffsize1;
	double* bsendbuff, bsendbuff1;
	int syncsend, syncrecv;

	long int memory;	//used to calculate memory;

	keyarray keyarr;	//array of keys (fields and methods defined im keys.h)
	int keynum;		//number of keys
	int keycounter; 	//counter of keys	

	int* NPLWKP_TOT;	//total number of coefficients in plane wave array for each kpoint
				//(probably the number of pl. waves may vary from one k-point to another)
	int NPLWKP = 0;

	double **RIG;		//an array of plane wave G-vectors RIG[plane wave number][coordinate(B0,B1,B2) number]
	int **IG;		//an array of plane wave G-vectors IG[plane wave number][coordinate(B0,B1,B2) number]

	complex<float> **CPTWFP1;	//initializing of an array for wavefunction CPTWFP[band number][plane wave number]
					//last dimension equals NPLWKP_TOT(K)and may be different for different k-points
	complex<float> **CPTWFP2;


	float **planewavevector = NULL;	//an array of plane wave wave vectors 
					//planewavevector[plane wave number][coordinate(x,y,z) number]

	double *y_efermi;	//used to calculate the Fermi energy; values of energy
	double *x_efermi;	//used to calculate the Fermi energy; ln((2 - weight) / weight)
	int npoints_efermi = 0;	//the number of points for the Fermi energy calculation; weights must be between 1.5 and 0.5
	int start_efermi = 0;	//the number of the first point with the weight less than 1.5
	int end_efermi = 0;		//the number of the last point with the weight larger than 0.5
	double upper_boundary_efermi = 0.9;
	double lower_boundary_efermi = 0.1;
	int counter_efermi;
	double xy_aver_efermi, x_aver_efermi, y_aver_efermi, x2_aver_efermi; //used to calculate the Fermi energy and temperature
										//by means of the least square method
	double Te;		//electron temperature (measured in eV)
	double efermi;		//the Fermi energy (measured in eV)

//	float** L12;		//Onsager coefficient L12 (proportionality between electric current density and temperature gradient)
				//measured in (A / m^2) / (K / m) 
				//L12[deltae point number][omega point number]
//	float*** L12k;		//L12 coefficient for a particular k-point
				//L12k[k-point number][deltae point number][omega point number]
//	float** L21;		//Onsager coefficient L21 (proportionality between heat current density and electric field)
				//measured in (W / m^2) / (V / m)
				//L21[deltae point number][omega point number]
				//L21 = L12 * T
//	float*** L21k;		//L21 coefficient for a particular k-point
				//L21k[k-point number][deltae point number][omega point number]
//	float** L22;		//Onsager coefficient L22 (proportionality between heat current density and temperature gradient)
				//measured in (W / m^2) / (K / m) = W / (K * m)
				//L22[deltae point number][omega point number]
//	float*** L22k;		//L22 coefficient for a particular k-point
				//L22k[k-point number][deltae point number][omega point number]

	string L12name = "L12";
	string L12formname = "L12FORM";
	string L21name = "L21";
	string L21formname = "L21FORM";
	string L22name = "L22";
	string L22formname = "L22FORM";

	conductivity sigma_k_proc;
	conductivity L12_k_proc;
	conductivity L21_k_proc;
	conductivity L22_k_proc;
	conductivity* sigma_k;
	conductivity* L12_k;
	conductivity* L21_k;
	conductivity* L22_k;

	conductivity sigma;
	conductivity L12;
	conductivity L21;
	conductivity L22;
	
	double sigma_increment;	//used to store the increment of conductivity at each step of the Kubo-Greenwood formula
	double e_av_minus_efermi;	//used to store ((e_k + e_k') / 2 - efermi) during Onsager coefficients calculation
	

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &size);
	MPI_Comm_rank(MPI_COMM_WORLD, &rank);

	if (rank == 0)
		time(&time_beg);

	keyarr.find_parse_keys(argc, argv);
	keynum = keyarr.keynum;
	for (keycounter = 0; keycounter < keynum; keycounter++)
	{
		if (keyarr.arr[keycounter].name.str()=="k")
		{
			if (keyarr.arr[keycounter].compound==1)
			{
				cerr << programname << ": -k key must not have parameters; program terminated" << endl;
				exit(1);
			}
			allkptswrite = 1;
		}
		else if (keyarr.arr[keycounter].name.str()=="kpath")
		{
			kpathpresent = 1;
			kpathname = keyarr.arr[keycounter].param.str();	
		}
		else if (keyarr.arr[keycounter].name.str()=="kstart")
		{
			keyarr.arr[keycounter].param >> kstart;
		}
		else if (keyarr.arr[keycounter].name.str()=="wavecar")
		{
			wavecarkeypresent = 1;
			wavecarname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="planewave")
		{
			planewavekeypresent = 1;
			planewavename = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="nablawrite")
		{
			nablawritepresent = 1;
			nablawritename = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="nablaread")
		{
			nablareadpresent = 1;
			nablareadname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="zerocorr")
		{
			keyarr.arr[keycounter].param >> zerocorrection;
		}
		else if (keyarr.arr[keycounter].name.str()=="omega")
		{
			omegapresent = 1;
			omeganame = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="onsager")
		{
			onsagerfound = 1;
		}
		else if (keyarr.arr[keycounter].name.str()=="resigma")
		{
			resigmapresent = 1;
			resigmaname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="resigmaform")
		{
			resigmaformpresent = 1;
			resigmaformname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="progress")
		{
			progressfound = 1;
			progressname = keyarr.arr[keycounter].param.str();
		}
		else if (keyarr.arr[keycounter].name.str()=="recoules")
		{
			recoulesfound = 1;
		}
		else if (keyarr.arr[keycounter].name.str()=="divideby")
		{
			dividebyfound = 1;
			keyarr.arr[keycounter].param >> divideby;
		}
		else
		{
			cerr << programname << ": Unknown key " << *(argv + keycounter + 1) << "    Program terminated" << endl;
			exit(1);
		}
	}
	
	if (dividebyfound == 1)
	{
		if ((divideby != 0) && (divideby != 1))
		{
			cerr << programname << ": WARNING: No rule for divideby = " << divideby << " found; divideby = 0 will be used instead" << endl;
			divideby = 0;
		}
	}
	else
	{
		divideby = 0;
	}
	
	if (omegapresent != 0)
	{


		omegainputresult = omegaobject.input_binary(omeganame.c_str());

		if (omegainputresult != 0)
		{
			if (rank == 0)
				cerr << programname << ": ERROR: omega input file " << omeganame << " was not opened, program terminated" << endl;
			MPI_Finalize();
			exit(1);
		}
				//initializing of an array of widths of gaussian broadening for delta function deltae[deltae point number]
				//initializing of an array of frequency omega[omega point number]
				//deltae and frequency are measured in electronvolts (eV)
		deltaenum = omegaobject.deltaenum;
		omeganum = omegaobject.omeganum;

		deltae = new double [deltaenum];
		omega = new double [omeganum];
	
		for (deltaecounter = 0; deltaecounter < deltaenum; deltaecounter++)
			deltae[deltaecounter] = omegaobject.deltae[deltaecounter];

		for (omegacounter = 0; omegacounter < omeganum; omegacounter++)
			omega[omegacounter] = omegaobject.omega[omegacounter];

		omegaobject.memory_remove();
	}
	else
	{
		deltaenum = 4;
		omeganum = 500;

		deltae = new double [deltaenum];
		omega = new double [omeganum];

		deltae[0] = 0.2;
		deltae[1] = 0.1;
		deltae[2] = 0.05;
		deltae[3] = 0.02;

		for (omegacounter = 0; omegacounter < 500; omegacounter++)
		{
			omega[omegacounter] = 0.02 + (omegacounter) * 0.02;
		}
	}

	complex<float> wfproduct; //used for wavefunction scalar product calculation
	complex<float> norm; 


				//-wavecar and -planewave keys must be present simultaneously
	if (((wavecarkeypresent) && (! planewavekeypresent)) || ((planewavekeypresent) && (! wavecarkeypresent)))
	{
		if (rank == 0)
			cerr << programname << ": ERROR: -wavecar and -planewave keys should be present simultaneously; program terminated" << endl;
		MPI_Finalize();
		exit(1);
	}	
				//only one of keys -wavecar and -nablaread should be present
	if ((wavecarkeypresent) && (nablareadpresent))
	{
		if (rank == 0)
			cerr << programname << ": ERROR: -wavecar and -nablaread keys should not be present simultaneously; program terminated" << endl;
		MPI_Finalize();
		exit(1);
	}

	if ((! wavecarkeypresent) && (! nablareadpresent))
	{
		if (rank == 0)
			cerr << programname << ": ERROR: either -wavecar xor -nablaread key should be present; program terminated" << endl;
		MPI_Finalize();
		exit(1);
	}

	if ((nablareadpresent) && (nablawritepresent))
	{
		if (rank == 0)
			cerr << programname << ": WARNING: -nablaread and -nablawrite keys are present; nabla matrix elements will be read and then written, this is strange" << endl;
	}

	if ((zerocorrection > 2) || (zerocorrection < 0))
	{
		cerr << programname << ": WARNING: the value of zerocorrection " << zerocorrection << " makes no sense and will be replaced by 2 value" << endl;
		zerocorrection = 2;
	}

	if (resigmapresent == 0)
	{
		resigmaname = "ReCONDUCTIVITY";
	}

	if (resigmaformpresent == 0)
	{
		resigmaformname = resigmaname + "FORM";
	}

	//reading WAVECAR file
	if (progressfound == 0)
	{
		progressname = "progress";
	}
	progress.open (progressname.c_str(), ios_base::out);



	if (wavecarkeypresent)
	{
		wavecar.open(wavecarname.c_str(), ios_base::in|ios_base::binary);
	
		wavecar.read((char*)&RDUM,sizeof(double));
		wavecar.read((char*)&RISPIN,sizeof(double));
		wavecar.read((char*)&RTAG,sizeof(double));
		IDUM = (int)(floor(RDUM + .5)); ISPIN = (int)(floor(RISPIN + .5)); ITAG = (int)(floor(RTAG + .5));

		wavecar.seekg(IDUM);
		wavecar.read((char*)&RNKPTS,sizeof(double));
		wavecar.read((char*)&RNB_TOT,sizeof(double));
		wavecar.read((char*)&ENMAX,sizeof(double));
		NKPTS = (int)(floor(RNKPTS + .5)); NB_TOT = (int)(floor(RNB_TOT + .5));

		if (rank == 0)
		{
			progress << "Process #0: NKPTS = " << NKPTS << " was read from " << wavecarname.c_str() << endl << endl;
		}

		NPLWKP_TOT = new int [NKPTS]; //total number of coefficients in plane wave array for each kpoint
						//(probably the number of pl. waves may vary from one k-point to another)
	}
	else if (nablareadpresent)
	{
		nablaread.open(nablareadname.c_str(), ios_base::in|ios_base::binary);
		if (rank == 0)
		{
			progress << "Process #0: Starting reading nabla matrix elements from " << nablareadname << endl;
			time (&time3);
		}
		nablaread.read((char*)&nablareadrdum, sizeof(double));
		nablaread.read((char*)&RISPIN, sizeof(double));
		nablaread.read((char*)&RNKPTS, sizeof(double));
		nablaread.read((char*)&RNB_TOT, sizeof(double));
		nablaread.read((char*)&LATT_INI, sizeof(latt));
		nablareadidum = (int)(floor(nablareadrdum + .5)); ISPIN = (int)(floor(RISPIN + .5)); 
		NKPTS = (int)(floor(RNKPTS + .5)); NB_TOT = (int)(floor(RNB_TOT + .5));
	}

	double **VKPT = new double*[NKPTS]; //initializing of an array of k-points coordinates VKPT[k-point number][coordinate(B0,B1,B2) number]
	for (K = 0; K < NKPTS; K++)
	{
		VKPT[K] = new double [3];
	}
	
	complex<double> ***CELTOT = new complex<double>** [ISPIN]; //initializing of an array of eigenvalues CELTOT[spin number][k-point number][band number]
	double ***FERTOT = new double** [ISPIN]; //initializing of an array of fermi-weights FERTOT[spin number][k-point number][band number]
	for (ISP = 0; ISP < ISPIN; ISP++)
	{
		CELTOT[ISP] = new complex<double>* [NKPTS];
		FERTOT[ISP] = new double* [NKPTS];
		for (K = 0; K < NKPTS; K++)
		{
			CELTOT[ISP][K] = new complex<double> [NB_TOT];
			FERTOT[ISP][K] = new double [NB_TOT];
		}
	}
					//initializing of an array of symmetry weight factors WTKPT[k-point number]
	double *WTKPT = new double [NKPTS];
	

	if (NB_TOT % size == 0)
		nb_proc_max = NB_TOT / size;
	else if (NB_TOT % size != 0)
		nb_proc_max = NB_TOT / size + 1;


	nablamatrelem_proc = new complex<float>** [nb_proc_max];
	for (J = 0; J < nb_proc_max; J++)
	{
		nablamatrelem_proc[J] = new complex<float>* [NB_TOT];
		for (J1 = 0; J1 < NB_TOT; J1++)
		{
			nablamatrelem_proc[J][J1] = new complex<float> [3];
		}
	}

	if (rank == 0)
	{
		progress << "Process #0: " << memory_to_GB(nb_proc_max * NB_TOT * 3 * sizeof(complex<float>)) << " GB allocated for nabla matrix elements" << endl;
	}

	if (wavecarkeypresent)
	{
						//reading the number of plane waves for each k-point
						//since it does not depend on spin number, only loop on k-point number is performed (spin number = 1)
		NPLWKP = 0;
		for (K = 0; K < NKPTS; K++) 
		{
			IREC = 2 + K * (NB_TOT + 1);
			wavecar.seekg(IREC*IDUM);
			wavecar.read((char*)&RNPL,sizeof(double));
			NPLWKP_TOT[K] = (int)(floor(RNPL + .5));
			if (NPLWKP < NPLWKP_TOT[K])
				NPLWKP = NPLWKP_TOT[K];
		}
							//initializing of an array for wavefunction CPTWFP[spin number][k-point number][band number][plane wave number]
							//last dimension equals NPLWKP_TOT(K)and may be different for different k-points

		CPTWFP1 = new complex<float>* [nb_proc_max];
		CPTWFP2 = new complex<float>* [nb_proc_max];
		for (J = 0; J < nb_proc_max; J++)
		{
			CPTWFP1[J] = new complex<float> [NPLWKP];
			CPTWFP2[J] = new complex<float> [NPLWKP];
		}
	
		memory = 2 * (long int)nb_proc_max * (long int)sizeof(complex<float>) * (long int)NPLWKP;
	
		if (rank == 0)
		{
			progress << "Process #0: Memory for wavefunctions created" << endl;
			progress << memory_to_GB(memory) << " GB allocated" << endl << endl;
		}

		for (K = 0; K < NKPTS; K++)
		{
						//reading the header for k-point with number K
			IREC = 2 + K * (NB_TOT + 1);
			wavecar.seekg(IREC*IDUM);
			wavecar.read((char*)&RNPL,sizeof(double));
			wavecar.read((char*)VKPT[K],3*sizeof(double));

			for (J = 0; J < NB_TOT; J++)
			{
				wavecar.read((char*)(CELTOT[0][K]+J),sizeof(complex<double>));
				wavecar.read((char*)(FERTOT[0][K]+J),sizeof(double));
			}
		}


		// reading PLANEWAVE file
		planewave.open(planewavename.c_str(),ios_base::in|ios_base::binary);
		
		if (rank == 0)
		{
			time(&time1);
			progress << "Process #0: Starting PLANEWAVE reading" << endl;
		}
		planewave.seekg(2*sizeof(double)); //first two numbers in PLANEWAVE are RDUM and RNKPTS, same as in WAVECAR, skip them
		planewave.read((char*)&LATT_INI,sizeof(latt));
		

				//initializing of an array of plane wave G-vectors RIG[plane wave number][coordinate(B0,B1,B2) number]
				//G-vectors are measured in reciprocal lattice vector units
				//initializing of an array of plane wave G-vectors IG[plane wave number][coordinate(B0,B1,B2) number]
				//IG is same as RIG, but integer
		RIG = new double* [NPLWKP];
		IG = new int* [NPLWKP];
		for (I = 0; I < NPLWKP; I++)
		{
			RIG[I] = new double [3];
			IG[I] = new int [3];
		}
						
				//initializing of an array of plane wave wave vectors 
				//planewavevector[k-point number][plane wave number][coordinate(x,y,z) number]
		planewavevector = new float* [NPLWKP];
		for (I = 0; I < NPLWKP; I++)
		{
			planewavevector[I] = new float [3];
		}
		
		memory = (long int)NPLWKP * ((long int)3 * (long int)sizeof(float) + (long int)3 * (long int)sizeof(double) + (long int)3 * (long int)sizeof(int)); 
		if (rank == 0)
		{
			progress << "Process #0: Memory for G-vectors created" << endl;
			progress << "Process #0: " << memory_to_GB(memory) << " GB allocated" << endl << endl;
		}


		for (K = 0; K < NKPTS; K++)
		{
			IREC = 1 + 5 * K;	//reading symmetry weight factor for a k-point
			planewave.seekg(IREC * IDUM + 4 * sizeof(double));
			planewave.read((char*)(WTKPT+K),sizeof(double));
		}


	}
	else if (nablareadpresent)
	{
		for (K = 0; K < NKPTS; K++)
		{	
			nablaread.seekg((1 + K * 4) * nablareadidum);
			nablaread.read((char*)VKPT[K], 3 * sizeof(double));
			nablaread.read((char*)(WTKPT+K), sizeof(double));
			for (J = 0; J < NB_TOT; J++)
			{		
				nablaread.read((char*)(CELTOT[0][K]+J),sizeof(complex<double>));
				nablaread.read((char*)(FERTOT[0][K]+J),sizeof(double));
			}
		}
	}



	int *Jbeg = new int[size];
	int *Jfin = new int[size];
	for (proccounter = 0; proccounter < size; proccounter++)
	{
		if (proccounter == 0)
		{
			Jbeg[proccounter] = 0;
		}
		else
		{
			Jbeg[proccounter] = Jfin[(proccounter - 1)] + 1;
		}
		if (proccounter < (NB_TOT % size))
		{
			Jfin[proccounter] = Jbeg[proccounter] + (NB_TOT / size) - 1 + 1;
		}
		else
		{
			Jfin[proccounter] = Jbeg[proccounter] + (NB_TOT / size) - 1;
		}
	}



	//Here main loop over k-points starts
	for (K = 0; K < NKPTS; K++)
	{
		if (rank == 0)
		{
			progress << "Process #0: calculation for k-point #" << K << " started:" << endl << endl;
		}

		if (wavecarkeypresent)
		{
			if ((rank ==0)&& (K != 0))
			{
				time(&time1);
			}
			
							//reading body of PLANEWAVE file
			for (a = 0; a < 3; a++)
			{
				IREC = 1 + 5 * K + 1 + a; //skip the header for a k-point (it consists of the values of NPLWKP_TOT[K],VKPT[K][coordinate(B0,B1,B2) number])
							//these values are the same as in WAVECAR file
				planewave.seekg(IREC*IDUM);
				for (I = 0; I < (NPLWKP_TOT[K]); I++)
				{
					planewave.read((char*)(&(RIG[I][a])),sizeof(double));
					IG[I][a] = (int)(floor((RIG[I][a]) + .5));
				}
				for (I = NPLWKP_TOT[K]; I < NPLWKP; I++)
				{
					RIG[I][a] = 0;
					IG[I][a] = 0;
				}
			}
	
			if (K == (NKPTS - 1))
				planewave.close();
		
			if (rank == 0)
			{
				time(&time2);
				progress << "Process #0: PLANEWAVE file for k-point #" << K << " read" << endl << "Time elapsed:" << endl;
				printtimehms(difftime(time2,time1), &progress);
				progress << endl;
			}
	
		
		//calculation of the matrix elements of the nabla
	
								//initializing of an array of nabla matrix elements
								//nablamatrelem[k-point number][left band number][right band number][coordinate(x,y,z) number]
		}
	
		else if (nablareadpresent)
		{
	
		}
	
		if (wavecarkeypresent)
		{
							//calculation of plane wave wave vectors  
							//unlike VKPT and RIG planewavevector is given in direct space (measured in A^-1)
							//unlike VASP planewavevector contains values already multiplied by 2*pi factor (correct (k+G)-vectors)
			if (rank == 0)
			{
				progress << "Process #0: G-vectors conversion starting" << endl;
				time(&time1);
			}


			for (I = 0; I < (NPLWKP_TOT[K]); I++)
			{
				for (a = 0; a < 3; a++)
				{
					planewavevector[I][a] = 0;
					for (b = 0; b < 3; b++)		//converting from reciprocal to direct coordinates
					{
						planewavevector[I][a] = planewavevector[I][a] + float((VKPT[K][b] + RIG[I][b]) * LATT_INI.B[b][a] );
					}
					planewavevector[I][a] = planewavevector[I][a] * 2 * M_PI;
				}
			}
			

			for (I = NPLWKP_TOT[K]; I < NPLWKP; I++)
			{
				for (a = 0; a < 3; a++)
				{
					planewavevector[I][a] = 0;
				}
			}
	
					//removing RIG and IG arrays 
					//they are not necessary after G-vectors conversion performed
			if (K == (NKPTS - 1))
			{
				for (I = 0; I < NPLWKP; I++)
				{
					delete [] RIG[I];
					delete [] IG[I];
				}
				delete [] RIG;
				delete [] IG;

				memory = (long int)NPLWKP * ((long int)3 * (long int)sizeof(int) + (long int)3 * (long int)sizeof(double)); 
				if (rank == 0)
				{
					progress << "Process #0: Memory for G-vectors released" << endl;
					progress << "Process #0: " << memory_to_GB(memory) << " GB released" << endl << endl;
				}
			}
		
			if (rank == 0)
			{
				time(&time2);
				progress << "Process #0: G-vectors conversion finished" << endl << "Time elapsed:" << endl;
				printtimehms(difftime(time2,time1), &progress);
				progress << endl;
			}
	
						//calculation of nabla matrix elements
						//matrelem(J,J1) = sum_G (C_k+G(J))^**i(k+G)*(C_k+G(J1))
						//plane waves are orthonormal
			if (rank == 0)
			{
				progress << "Process #0: Nabla matrix elements calculation starting" << endl;
				time(&time1);
				timesectioncounter = 0;
			}
	
	
	

			//the information on wave functions is read
			//the informaltion is read for bands, which will be permanently stored on this node
			for (J1proc = 0; J1proc < (Jfin[rank] - Jbeg[rank] + 1); J1proc++)
			{
				IREC = 2 + 0 + K * (NB_TOT + 1) + 1 + Jbeg[rank] + J1proc;
				wavecar.seekg(IREC * IDUM);
				wavecar.read((char*)CPTWFP1[J1proc], (NPLWKP_TOT[K])*sizeof(complex<float>));
				for (I = NPLWKP_TOT[K]; I < NPLWKP; I++)
				{
					CPTWFP1[J1proc][I] = (0,0);
				}
			}
			for (J1proc = (Jfin[rank] - Jbeg[rank] + 1); J1proc < nb_proc_max; J1proc++)
			{
				for (I = NPLWKP_TOT[K]; I < NPLWKP; I++)
				{
					CPTWFP1[J1proc][I] = (0,0);
				}
			}
	
			for (proccounter = 0; proccounter < size; proccounter++)
			{

			
				//the information on wave functions is read
				//the information is read for bands, 
				//which will be temporarily stored on all nodes at current iteration
				for (J2proc = 0; J2proc < (Jfin[proccounter] - Jbeg[proccounter] + 1); J2proc++)
				{
					IREC = 2 + 0 + K * (NB_TOT + 1) + 1 + Jbeg[proccounter] + J2proc;
					wavecar.seekg(IREC * IDUM);
					wavecar.read((char*)CPTWFP2[J2proc], (NPLWKP_TOT[K])*sizeof(complex<float>));
					for (I = NPLWKP_TOT[K]; I < NPLWKP; I++)
					{
						CPTWFP2[J2proc][I] = (0,0);
					}
				}
				for (J2proc = (Jfin[proccounter] - Jbeg[proccounter] + 1); J2proc < nb_proc_max; J2proc++)
				{
					for (I = NPLWKP_TOT[K]; I < NPLWKP; I++)
					{
						CPTWFP2[J2proc][I] = (0,0);
					}
				}
			
				for (J = Jbeg[rank]; J < (Jfin[rank] + 1); J++)
				{
					for (J1 = Jbeg[proccounter]; J1 < (Jfin[proccounter] + 1); J1++)
					{
						for (a = 0; a < 3; a++)
						{
							nablamatrelem_proc[J-Jbeg[rank]][J1][a] = complex<float> (0,0);
							for (I = 0; I < NPLWKP_TOT[K]; I++)
							{
								nablamatrelem_proc[J-Jbeg[rank]][J1][a] = nablamatrelem_proc[J-Jbeg[rank]][J1][a] + conj(CPTWFP1[J-Jbeg[rank]][I]) * (complex<float> (0,planewavevector[I][a])) * (CPTWFP2[J1-Jbeg[proccounter]][I]);

							}
						}
					}
					if (rank == 0)
					{
						time(&time2);
						if ((((int)(difftime(time2, time1)) / (timinginterval * 60)) - timesectioncounter)!=0)
						{
							progress << "Process #0: " << (100 * ((Jfin[0] - Jbeg[0] + 1) * (Jbeg[proccounter] - Jbeg[0]) + (J - Jbeg[0] + 1) * (Jfin[proccounter] - Jbeg[proccounter] + 1))) / ((Jfin[0] - Jbeg[0] + 1) * NB_TOT * NKPTS) << " % acomplished" << endl;
							timesectioncounter = (int)(difftime(time2,time1)) / (timinginterval * 60);
						}
					
					}
				}
			
			}

			MPI_Barrier(MPI_COMM_WORLD);
//			progress << "rank = " << rank << endl;
	
			if (K == (NKPTS - 1))		
				wavecar.close();
	

			if (K == (NKPTS - 1))
			{
				//removing CPTWFP array
				for (J = 0; J < nb_proc_max; J++)
				{
					delete [] CPTWFP1[J];
					delete [] CPTWFP2[J];
				}
				delete [] CPTWFP1;
				delete [] CPTWFP2;
	
				memory = 2 * (long int)nb_proc_max * (long int)sizeof(complex<float>) * (long int)NPLWKP; 
				if (rank == 0)
				{
					progress << endl << "Process #0: Wave functions removed" << endl;
					progress << "Process #0: " << memory_to_GB(memory) << " GB released" << endl << endl;
				}
	
	
				//removing planewavevector array
				for (I = 0; I < NPLWKP; I++)
				{
					delete [] planewavevector[I];
				}
				delete [] planewavevector;
	
				memory = (long int)NPLWKP * ((long int)3 * (long int)sizeof(float)); 
				if (rank == 0)
				{
					progress << "Process #0: Memory for G-vectors released" << endl;
					progress << "Process #0: " << memory_to_GB(memory) << " GB released" << endl << endl;
				}
			}
	
			if (rank == 0)
			{
				time(&time2);
				progress << "Process #0: Nabla matrix elements calculation finished" << endl << "Time elapsed:" << endl;
				printtimehms(difftime(time2,time1), &progress);
				progress << endl;
			}
	
//			progress << "rank = " << rank << endl;
	
		}
		else if (nablareadpresent)
		{
			for (a = 0; a < 3; a++)
			{
				nablaread.seekg((1 + K * 4 + a + 1) * nablareadidum + Jbeg[rank] * NB_TOT * 2 * sizeof(double));
				for (J = Jbeg[rank]; J < Jfin[rank] + 1; J++)
				{
					for (J1 = 0; J1 < NB_TOT; J1++)
					{
						nablaread.read((char*)&convert_dbl_to_flt, sizeof(double));
						nablamatrelem_proc[J-Jbeg[rank]][J1][a].real((float)convert_dbl_to_flt);
						nablaread.read((char*)&convert_dbl_to_flt, sizeof(double));
						nablamatrelem_proc[J-Jbeg[rank]][J1][a].imag((float)convert_dbl_to_flt);
					}
				}
				for (J = Jfin[rank] + 1; J < nb_proc_max; J++)
				{
					for (J1 = 0; J1 < NB_TOT; J1++)
					{
						nablamatrelem_proc[J][J1][a] = (0, 0);
					}
				}
			}
	
			if (K == (NKPTS - 1))
				nablaread.close();
	
			if (rank == 0)
			{
				time(&time4);
				progress << "Process #0: Nabla matrix elements read for k-point #" << K << endl << "Time elapsed:" << endl;
				printtimehms(difftime(time4,time3), &progress);
				progress << endl;
			}
			
		}
	
//				progress << "proccounter = " << proccounter << " rank = " << rank << endl;
		if (nablawritepresent)
		{
			nablawriteidum = NB_TOT * NB_TOT * sizeof(complex<double>);	
			for (proccounter = 0; proccounter < size; proccounter++)
			{
//				progress << "proccounter = " << proccounter << " rank = " << rank << endl;

				if (proccounter == rank)
				{
					if ((rank == 0) && (K == 0))
					{
						nablawrite.open(nablawritename.c_str(), ios::out|ios::binary);
					}
					else
					{
						nablawrite.open(nablawritename.c_str(), ios::in|ios::out|ios::binary);
					}
					
					if (rank == 0)
					{
						//the header of the whole file is written by the zero process
						if (K == 0)
						{
							nablawrite.seekp(0);
							convert = nablawriteidum;
							nablawrite.write ((char*)&convert, sizeof(double));
							nablawrite.write ((char*)&RISPIN, sizeof(double));
							nablawrite.write ((char*)&RNKPTS, sizeof(double));
//							nablawrite.close();
//							nablawrite.open(nablawritename.c_str(), ios::in|ios::out|ios::binary);
//							nablawrite.seekp(3*sizeof(double));
							nablawrite.write ((char*)&RNB_TOT, sizeof(double));
							nablawrite.write ((char*)&LATT_INI, sizeof(latt));
						}


						//the header for each k-point is written by the zero process
						nablawrite.seekp((1 + K * 4) * nablawriteidum);
						nablawrite.write ((char*)VKPT[K], 3 * sizeof(double));
						nablawrite.write((char*)(WTKPT+K),sizeof(double));
						for (J = 0; J < NB_TOT; J++)
						{
							nablawrite.write((char*)(CELTOT[0][K]+J),sizeof(complex<double>));
							nablawrite.write((char*)(FERTOT[0][K]+J),sizeof(double));
						}
					}

					//the body of the file corresponding to processor #rank
					for (a = 0; a < 3; a++)
					{
						nablawrite.seekp((1 + K * 4 + a + 1) * nablawriteidum + Jbeg[rank] * NB_TOT * 2 * sizeof(double));
						for (J = Jbeg[rank]; J < Jfin[rank] + 1; J++)
						{
							for (J1 = 0; J1 < NB_TOT; J1++)
							{
								convert = real(nablamatrelem_proc[J-Jbeg[rank]][J1][a]);
								nablawrite.write((char*)&convert,sizeof(double));
								convert = imag(nablamatrelem_proc[J-Jbeg[rank]][J1][a]);
								nablawrite.write((char*)&convert,sizeof(double));
							}
						}
					}
					nablawrite.close();
				}
				MPI_Barrier(MPI_COMM_WORLD);
			}

		}
	
		if ((onsagerfound == 1) && (K == 0))
		{
			if (rank == 0)
			{
				progress << "Calculation of chemical potential and electron temperature starting" << endl;
				time(&time1);
			}
	
			start_efermi = 0;
			while (FERTOT[0][0][start_efermi] > upper_boundary_efermi)
			{
				start_efermi++;
			}
			end_efermi = start_efermi;
			while (FERTOT[0][0][end_efermi] > lower_boundary_efermi)
			{
				end_efermi++;
			}
			npoints_efermi = end_efermi - start_efermi;
	
			if (npoints_efermi < 2)
			{
				if (rank == 0)
					cerr << programname << ": WARNING: There are " << npoints_efermi << " bands with fermi-weights between " << upper_boundary_efermi << " and " << lower_boundary_efermi << endl;
	
				if (start_efermi - 4 < 0)
					start_efermi = 0;
				else
					start_efermi = start_efermi - 4;
	
				if (end_efermi + 4 > (NB_TOT - 1))
					end_efermi = NB_TOT;
				else
					end_efermi = end_efermi + 4;
	
				npoints_efermi = end_efermi - start_efermi;
	
				if (npoints_efermi < 2)
				{
					if (rank == 0)
						cerr << programname << ": ERROR: Only one band present, impossible to calculate chemical potential and temperature of electrons ; program terminated" << endl;
					MPI_Finalize();
					exit(1);
				}
				else
				{
					if (rank == 0)
						cerr << programname << ": WARNING: Bands with fermi-weights between " << FERTOT[0][0][start_efermi] << " and " << FERTOT[0][0][end_efermi-1] << " will be used instead" << endl;
				}
			}
			
			x_efermi = new double [npoints_efermi];
			y_efermi = new double [npoints_efermi];
	
			for (counter_efermi = 0; counter_efermi < npoints_efermi; counter_efermi++)
			{
					x_efermi[counter_efermi] = log ((1 - FERTOT[0][0][start_efermi + counter_efermi]) / FERTOT[0][0][start_efermi + counter_efermi]);
					y_efermi[counter_efermi] = real(CELTOT[0][0][start_efermi + counter_efermi]);
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
	
			Te = (xy_aver_efermi - x_aver_efermi * y_aver_efermi) / (x2_aver_efermi - x_aver_efermi * x_aver_efermi);
			efermi = y_aver_efermi - Te * x_aver_efermi;
		
			if (rank == 0)
			{
				progress << "Chemical potential and electron temperature calculated" << endl;
				progress << npoints_efermi << " bands with fermi-weights between " << FERTOT[0][0][start_efermi] << " and " << FERTOT[0][0][end_efermi - 1] << " used" << endl;
				progress << "Chemical potential (efermi) = " << efermi << " eV" << endl;
				progress << "Temperature of electrons = " << Te * 11604.505 << " K" << endl;
				time(&time2);
				progress << "Time elapsed:" << endl;
				printtimehms(difftime(time2,time1), &progress);
				progress << endl;
				
			}	
	
			delete [] x_efermi;
			delete [] y_efermi;
		}

			//now memory for conductivity and Onsager coefficients (if necessary) will be created
		if (K == 0)
		{
			sigma.memory_create(deltaenum,omeganum);
			for (deltaecounter = 0; deltaecounter < deltaenum; deltaecounter++)
			{
				sigma.deltae[deltaecounter] = deltae[deltaecounter];
			}
			for (omegacounter = 0; omegacounter < omeganum; omegacounter++)
			{
				sigma.omega[omegacounter] = omega[omegacounter];
			}
			
			sigma_k_proc.copy_frame(sigma);

			sigma_k = new conductivity [NKPTS];
			for (K1 = 0; K1 < NKPTS; K1++)
			{
				sigma_k[K1].copy_frame(sigma);
			}
			
			if (onsagerfound == 1)
			{
				L12.copy_frame(sigma);
				L21.copy_frame(sigma);
				L22.copy_frame(sigma);
				L12_k_proc.copy_frame(sigma);
				L21_k_proc.copy_frame(sigma);
				L22_k_proc.copy_frame(sigma);

				L12_k = new conductivity [NKPTS];
				L21_k = new conductivity [NKPTS];
				L22_k = new conductivity [NKPTS];
				for (K1 = 0; K1< NKPTS; K1++)
				{
					L12_k[K1].copy_frame(sigma);
					L21_k[K1].copy_frame(sigma);
					L22_k[K1].copy_frame(sigma);
				}
			}

			if (rank == 0)
			{
				memory = (2 + NKPTS) * sigma.get_arrays_size();
				if (onsagerfound == 1)
				{
					memory = memory * 4;
				}
				progress << "Process #0: objects for conductivity";
				if (onsagerfound == 1)
					progress << " and Onsager coefficients";
				progress <<" calculation created" << endl;
				progress << "Process #0: " <<memory_to_GB(memory) << " GB allocated for arrays in conductivity";
				if (onsagerfound == 1)
					progress << " and Onsager coefficients";				
				progress << " objects" << endl << endl;

			}

		}
		

		
					//conductivity calculation
		if (rank == 0)
		{
			progress << "Process #0: Conductivity ";
			if (onsagerfound == 1)
			{
				progress << "and Onsager coefficients ";
			}
			progress << "calculation starting for k-point #" << K << endl;
			time(&time1);
			timesectioncounter = 0;
		}



		for (deltaecounter = 0; deltaecounter < deltaenum; deltaecounter++)
		{
			for (omegacounter = 0; omegacounter < omeganum; omegacounter++)
			{
				sigma_k_proc.sigma[deltaecounter][omegacounter] = 0;
				if (onsagerfound == 1)
				{
					L12_k_proc.sigma[deltaecounter][omegacounter] = 0;
					L21_k_proc.sigma[deltaecounter][omegacounter] = 0;
					L22_k_proc.sigma[deltaecounter][omegacounter] = 0;
				}
				for (J = Jbeg[rank]; J < Jfin[rank] + 1; J++)
				{
					for (J1 = 0; J1 < NB_TOT; J1++)
					{
						fermidiff = (FERTOT[0][K][J] - FERTOT[0][K][J1]);
						if ((real(CELTOT[0][K][J1])-real(CELTOT[0][K][J])) < 0)
						{
							if (zerocorrection == 0)
								fermidiff = fermidiff * (-1);
							else if (zerocorrection == 1)
								fermidiff = 0;
							else if (zerocorrection == 2)
								fermidiff = fermidiff;
						}

						if (divideby == 0)
							dividebyvalue = omega[omegacounter];
						else if (divideby == 1)
						{
							dividebyvalue = fabs(real(CELTOT[0][K][J1])-real(CELTOT[0][K][J]));
							if (dividebyvalue < 1e-10)
							{
								dividebyvalue = omega[omegacounter];
							}
						}
						
						sigma_increment = 0;
						for (a = 0; a < 3; a++)
						{
							sigma_increment += real((nablamatrelem_proc[J-Jbeg[rank]][J1][a] * conj(nablamatrelem_proc[J-Jbeg[rank]][J1][a])));
						}
					
						sigma_increment *= fermidiff * (deltagaussian((real(CELTOT[0][K][J1]) - real(CELTOT[0][K][J]) - omega[omegacounter]), deltae[deltaecounter])) / dividebyvalue;
						sigma_k_proc.sigma[deltaecounter][omegacounter] += sigma_increment;

						if (onsagerfound == 1)
						{
							if (recoulesfound == 1)
							{
								L12_k_proc.sigma[deltaecounter][omegacounter] += (-1) * sigma_increment * (real(CELTOT[0][K][J]) - efermi);
								L21_k_proc.sigma[deltaecounter][omegacounter] += (-1) * sigma_increment * (real(CELTOT[0][K][J1]) - efermi);
								L22_k_proc.sigma[deltaecounter][omegacounter] += sigma_increment * (real(CELTOT[0][K][J]) - efermi) * (real(CELTOT[0][K][J1]) - efermi);
							}
							else 
							{
								e_av_minus_efermi = (real(CELTOT[0][K][J]) + real(CELTOT[0][K][J1])) / 2 - efermi;
								L12_k_proc.sigma[deltaecounter][omegacounter] += (-1) * sigma_increment * e_av_minus_efermi;
								L21_k_proc.sigma[deltaecounter][omegacounter] += (-1) * sigma_increment * e_av_minus_efermi;
								L22_k_proc.sigma[deltaecounter][omegacounter] += sigma_increment * e_av_minus_efermi * e_av_minus_efermi;
							}
						}
					}
				}
				

				sigma_k_proc.sigma[deltaecounter][omegacounter] *= ((double)1 / (double)3) * 8.880338 * 100000000 / LATT_INI.OMEGA;
				
				if (onsagerfound == 1)
				{
					
					L12_k_proc.sigma[deltaecounter][omegacounter] *= ((double)1 / (Te * 11604.505)) * ((double)1 / (double)3) * 8.880338 * 100000000 / LATT_INI.OMEGA;
					L21_k_proc.sigma[deltaecounter][omegacounter] *= ((double)1 / (double)3) * 8.880338 * 100000000 / LATT_INI.OMEGA;
					L22_k_proc.sigma[deltaecounter][omegacounter] *= ((double)1 / (Te * 11604.505)) * ((double)1 / (double)3) * 8.880338 * 100000000 / LATT_INI.OMEGA;
				}

				if (rank == 0)
				{
					time(&time2);
		                        if ((((int)(difftime(time2, time1)) / (timinginterval * 60)) - timesectioncounter)!=0)
		      	                {
	                	                progress << "Process #0: "<< (100 * (deltaecounter * omeganum + omegacounter  + 1)) / (deltaenum * omeganum) << " % acomplished" << endl;
	                        	        timesectioncounter = (int)(difftime(time2,time1)) / (timinginterval * 60);
	                      		}
				}
			}
	
   			
			MPI_Barrier(MPI_COMM_WORLD);

			MPI_Reduce(sigma_k_proc.sigma[deltaecounter], sigma_k[K].sigma[deltaecounter], omeganum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			if (onsagerfound == 1)
			{
				MPI_Reduce(L12_k_proc.sigma[deltaecounter], L12_k[K].sigma[deltaecounter], omeganum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(L21_k_proc.sigma[deltaecounter], L21_k[K].sigma[deltaecounter], omeganum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
				MPI_Reduce(L22_k_proc.sigma[deltaecounter], L22_k[K].sigma[deltaecounter], omeganum, MPI_DOUBLE, MPI_SUM, 0, MPI_COMM_WORLD);
			}
		}


		if (allkptswrite == 1)
		{
			if (rank == 0)
			{
				kpointfilenum.str("");
				kpointfilenum << (K + kstart);
				if (kpathpresent)
				{
					kpointfilename = kpathname;
				}
				else
				{
					kpointfilename = resigmaname;
				}
				kpointfilename = kpointfilename + "_kpoint" + kpointfilenum.str();
				sigma_k[K].output_binary(kpointfilename.c_str());
				
				if (kpathpresent)
				{
					kpointfilename = kpathname + "FORM";
				}
				else
				{
					kpointfilename = resigmaformname;
				}
				kpointfilename = kpointfilename + "_kpoint" + kpointfilenum.str();
				sigma_k[K].output_form(kpointfilename.c_str(), "Conductivity calculated", "conductivity (1 / (Ohm * m))");
				if (onsagerfound == 1)
				{
					kpointfilename = L12name + "_kpoint" + kpointfilenum.str();
					L12_k[K].output_binary(kpointfilename.c_str());
					kpointfilename = L12formname + "_kpoint" + kpointfilenum.str();
					L12_k[K].output_form(kpointfilename.c_str(), "L12 Onsager coefficient calculated: ", "L12, (A / m^2) / (K / m)");
					kpointfilename = L21name + "_kpoint" + kpointfilenum.str();
					L21_k[K].output_binary(kpointfilename.c_str());
					kpointfilename = L21formname + "_kpoint" + kpointfilenum.str();
					L21_k[K].output_form(kpointfilename.c_str(), "L21 Onsager coefficient calculated: ", "L21, (W / m^2) / (V / m)");
					kpointfilename = L22name + "_kpoint" + kpointfilenum.str();
					L22_k[K].output_binary(kpointfilename.c_str());
					kpointfilename = L22formname + "_kpoint" + kpointfilenum.str();
					L22_k[K].output_form(kpointfilename.c_str(), "L22 Onsager coefficient calculated: ", "L22, W  / (m * K)");
				}
			}
		}

		if (rank == 0)
		{
			time(&time2);
			progress << "Process #0: Conductivity";
			if (onsagerfound == 1)
			{
				progress << " and Onsager coefficients";
			}
			progress << " calculated for k-point #"<< K << endl << "Time elapsed:" << endl;
			printtimehms(difftime(time2,time1), &progress);
			progress << endl;
		}

	}

		//weighted sum over k-points
	for (K = 0; K < NKPTS; K++)
	{
		sigma = sigma + mul_double(sigma_k[K], WTKPT[K]);
	}	

	if (onsagerfound == 1)
	{
		for (K = 0; K < NKPTS; K++)
		{
			L12 = L12 + mul_double(L12_k[K], WTKPT[K]);
			L21 = L21 + mul_double(L21_k[K], WTKPT[K]);
			L22 = L22 + mul_double(L22_k[K], WTKPT[K]);
		}
	}
	
	if (rank == 0)
	{
		progress << "Process #0: Conductivity output starting" << endl;
		time(&time1);
		timesectioncounter = 0;
		
		sigma.output_binary(resigmaname.c_str());

		sigma.output_form(resigmaformname.c_str(), "Conductivity calculated: ", "conductivity (1 / (Ohm * m))");
		
		if (onsagerfound == 1)
		{

			L12.output_binary(L12name.c_str());
			L12.output_form(L12formname.c_str(), "L12 Onsager coefficient calculated: ", "L12, (A / m^2) / (K / m)");
			L21.output_binary(L21name.c_str());
			L21.output_form(L21formname.c_str(), "L21 Onsager coefficient calculated: ", "L21, (W / m^2) / (V / m)");
			L22.output_binary(L22name.c_str());
			L22.output_form(L22formname.c_str(), "L22 Onsager coefficient calculated: ", "L22, W / (m * K)");

			L12 = mul_double(L12, Te * 11604.505);
			L12.output_binary("L12Te");
			L12.output_form("L12TeFORM", "L12 * Te Onsager coefficient calculated: ", "L12 * Te, (A / m)");
		}

	
		time(&time2);
		progress << "Process #0: Conductivity";
		if (onsagerfound == 1)
		{
			progress << " and Onsager coefficients";
		}
		progress << " printed" << endl << "Time elapsed:" << endl;
		printtimehms(difftime(time2,time1), &progress);
		progress << endl;
	}


	//now memory for conductivity and Onsager coefficients will be removed
	memory = (2 + NKPTS) * sigma.get_arrays_size();
	sigma.memory_remove();
	sigma_k_proc.memory_remove();
	if (onsagerfound == 1)
	{
		L12.memory_remove();
		L21.memory_remove();
		L22.memory_remove();
		L12_k_proc.memory_remove();
		L21_k_proc.memory_remove();
		L22_k_proc.memory_remove();
	}
	for (K1 = 0; K1 < NKPTS; K1++)
	{
		sigma_k[K1].memory_remove();
		if (onsagerfound == 1)
		{
			L12_k[K1].memory_remove();
			L21_k[K1].memory_remove();
			L22_k[K1].memory_remove();
		}
	}
	delete [] sigma_k;
	if (onsagerfound == 1)
	{
		delete [] L12_k;
		delete [] L21_k;
		delete [] L22_k;
	}

	if (rank == 0)
	{
		if (onsagerfound == 1)
		{
			memory = memory * 4;
		}
		progress << "Process #0: objects for conductivity";
		if (onsagerfound == 1)
			progress << " and Onsager coefficients";
		progress <<" calculation removed" << endl;
		progress << "Process #0: " << memory_to_GB(memory) << " GB released for arrays in conductivity";
		if (onsagerfound == 1)
			progress << " and Onsager coefficients";				
		progress << "  objects" << endl << endl;

	}



			//removing nablamatrelem array
	for (J = 0; J < nb_proc_max; J++)
	{
		for (J1 = 0; J1 < NB_TOT; J1++)
		{
			delete [] nablamatrelem_proc[J][J1];
		}
		delete [] nablamatrelem_proc[J];
	}
	delete [] nablamatrelem_proc;
	if (rank == 0)
	{
		progress << "Process #0: Memory for nabla matrix elements removed" << endl;
		progress << "Process #0: " << memory_to_GB(nb_proc_max * NB_TOT * 3 * sizeof(complex<float>)) << " GB released" << endl;
	}



	delete [] WTKPT;

					//removing CELTOT and FERTOT arrays
	for (ISP = 0; ISP < ISPIN; ISP++) 
	{
		for (K = 0; K < NKPTS; K++)
		{
			delete [] CELTOT[ISP][K];
			delete [] FERTOT[ISP][K];
		}
		delete [] CELTOT [ISP];
		delete [] FERTOT [ISP];
	}
	delete [] CELTOT;
	delete [] FERTOT;

	for (K = 0; K < NKPTS; K++) 
	{
		delete [] VKPT[K];
	}
	delete [] VKPT;

	if (wavecarkeypresent)
	{
		delete [] NPLWKP_TOT;
	}

	delete [] deltae;
	delete [] omega;
	
	keyarr.memory_remove();

	if (rank == 0)
	{
		time(&time_end);
		progress << "Program finished" << endl;
		progress << "Total time elapsed:" << endl;
		printtimehms(difftime(time_end, time_beg), &progress);
		progress << "Total time elapsed (seconds): " << endl << difftime(time_end, time_beg) << endl;
		progress << endl;
	}
	
	progress.close();

	MPI_Finalize();

	return 0;
}
