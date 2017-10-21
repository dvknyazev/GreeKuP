#ifndef OMEGAINOUT_H
#define OMEGAINOUT_H

class omegaclass
{
	public:
		int deltaenum;		//number of deltae points
		int omeganum;		//number of omega points
		double* deltae;		//an array of deltae points deltae[deltae point number]
		double* omega;		//an array of omega points omega[omega point number]
		int memoryexists;	//used to show if memory is allocated or now

		omegaclass();

		int memory_create(int, int);

		int input_binary(const char*);

		int output_form(const char*, const char*);

		int output_binary(const char*);

		int memory_remove();	//removes memory (if exists), sets memoryexists to 0
		~omegaclass();
};

#endif
