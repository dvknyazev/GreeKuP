#ifndef MYINOUT_H
#define MYINOUT_H
class conductivity
{
	friend conductivity operator+(conductivity, conductivity);	//overloaded operator + for conductivity class
								//checks the dimensions of operands, if they coincide, 
								//sums sigma values of operands
	
	friend conductivity operator-(conductivity, conductivity);	//overloaded operator - for conductivity class
								//checks the dimensions of operands, if they coincide, 
								//substracts sigma value of 2nd operand from the 1st one

	friend conductivity operator*(conductivity, conductivity);	//overloaded operator + for conductivity class
								//checks the dimensions of operands, if they coincide, 
								//multiplies sigma values of operands

	friend conductivity operator/(conductivity, conductivity);	//overloaded operator + for conductivity class
								//checks the dimensions of operands, if they coincide, 
								//divides sigma values of operands

	friend conductivity mul_double(conductivity, double);		//multiplicates conductivity argument by double number
	friend conductivity div_double(conductivity, double);		//divides conductivity argument by double number
	public:
		int deltaenum;		//number of elements in deltae array
		int omeganum;		//number of elements in omega array
		double* deltae;		//an array of deltae values deltae[deltae point number]
		double* omega;		//an array of omega values omega[omega point number]
		double** sigma;		//an array of conductivity values sigma[deltae point number][omega point number]
		int memoryexists;	//used to show if memory is allocated (1) or not (0)
	
		conductivity();
		int memory_create(int, int);	//memory_create(int newdeltaenum, int newomeganum)
						//if memory does not exists, creates it (memory size newdeltaenum*newomeganum)
						//if size of memory differs from newdeltenum*newomeganum, recreates it
						//sets deltaenum to newdeltaenum, omeganum to newomeganum, memoryexists to 1	
		int input_binary(const char*);	//input_binary (char* inputfilename)
						//reads data from CONDUCTIVITY file
						//if memory does not exists (or has another size) creates it
		int output_binary(const char*);	//output_binary (const char* outputfilename)
						//writes data to CONDUCTIVITY file
		int output_form(const char*, const char*, const char*);	//output_form (const char* outputfile, const char* header, const char* valueheader)
							//writes data to CONDUCTIVITYFORM file
		int output_form_transpose(const char*, const char*, const char*); //output_form_transpose (const char* outputfile, const char* header, const char* valueheader)
							//writes transposed data to formatted file
				
		int copy_frame(conductivity);	//copy_frame(conductivity cond)
						//creates (or recreates) memory according to the size of cond
						//copies deltae and omega values from cond; sets sigma values to zero
		int check_dimensions(conductivity);	//int check_dimensions (conductivity cond)
						//checks dimensions of object, called by and cond
						//if (one of) dimensions do not coincide or memory in one of objects does not exist,
						//returns 0
						//else returns 1
		conductivity operator=(conductivity);	//overloaded assignment operator = for conductivity
							//checks dimensions of left and right operand
							//then copies right operand to the left

		
		int resigma_to_imeps(conductivity);	//resigma_to_imepsilon(conductivity realcond)
							//calculates imaginary part of dielectric function
							//according to the formula Im eps = Re sigma / (omega * 13440)
							//sigma is measured in Ohm^-1*m^-1; omega - in eV; epsilon - dimensionless;
						//creates enough memory and copies deltae and omega values first (calls copy_frame)
							//writes data to the object, called by
							//takes real sigma values from realcond


		int imsigma_to_reeps(conductivity);	//imsigma_to_reeps(conductivity imcond)
							//calculates real part of dielectric function
							//according to the formula Re eps = 1 - (Im sigma / (omega * 13440))
							//sigma is measured in Ohm^-1*m^-1; omega - in eV; epsilon - dimensionless;
						//creates enough memory and copies deltae and omega values first (calls copy_frame)
							//writes data to the object, called by
							//takes im sigma values from imcond

		int rerefrac_calc(conductivity, conductivity);	//rerefrac_calc(conductivity reeps, conductivity imeps)
							//calculates the real part of the index of refraction
							//according to the formula 
							//n = (1 / sqrt(2)) * sqrt(sqrt(reeps^2 + imeps^2) + reeps)
							//n, reeps, imeps - dimensionless
							//creates enough memory and copies deltae and omega values first 
							//(calls copy_frame)
							//writes data to the object, called by

		int imrefrac_calc(conductivity, conductivity);	//imrefrac_calc(conductivity reeps, conductivity imeps)
							//calculates the imaginary part of the index of refraction
							//according to the formula 
							//k = (1 / sqrt(2)) * sqrt(sqrt(reeps^2 + imeps^2) - reeps)
							//k, reeps, imeps - dimensionless
							//creates enough memory and copies deltae and omega values first 
							//(calls copy_frame)
							//writes data to the object, called by
		
		int reflect_calc(conductivity, conductivity);	//reflect_calc(conductivity rerefrac, conductivity imrefrac)
							//calculates the reflectivity according to the formula
							//r = (((1 - n)^2 + k^2) / ((1 + n)^2 + k^2))
							//Fresnel formula in the case of the normal incidence
							//r, n, k - dimensionless
							//creates enough memory and copies deltae and omega values first
							//(calls copy_frame)
							//writes data to the object, called by

		int absorp_calc(conductivity);		//absorb_calc(conductivity imrefrac)
							//calculates the absorption coefficient according to the formula
							//alpha = 2 * imrefrac * omega / c 
							//claculation formula: imrefrac * omega * 10135400
							//where omega is measured -n eV, alpha in m^-1
							//creates enough memory and cpies deltae and omega values first
							//(calls copy_frame)
							//writes data to the object, called by

		int square_root(conductivity);
							//sqrt(conductivity inputcond)
							//calculates square root of conductivity
							//(necessary while dispersion calculation)

		int check_deltae(double);		//int check_deltae(double deltaevalue)
							//runs deltae array and searches for the first value
							//that differs from the deltaevalue less then 1E-10
							//if such value is found, 
							//the position of the value in the deltae array is returned
							//else -1 is returned

		int check_omega(double);		//int check_omega(double omegavalue)
							//runs omega array and searches for the first value
							//that differs from the omegavalue less then 1E-10
							//if such value is found, 
							//the position of the value in the omega array is returned
							//else -1 is returned
		int output_form_Origin(const char*, const char*, const char*, const char*, const char*, const char*, const char*); 
							//int output_form_Origin(const char* outputfilename, const char* xoriginstring, const char* xdimension, const char* yoriginstring, const char* ydimension, const char* columnheader, const char* columndimension);

		long int get_arrays_size();

		int memory_remove();	//removes memory for sigma, omega, deltae; sets omeganum, deltaenum, memoryexists to zero
		~conductivity();
};

#endif
