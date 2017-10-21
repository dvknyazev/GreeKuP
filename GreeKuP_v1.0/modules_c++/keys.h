#ifndef KEYS_H
#define KEYS_H

class key
{
	public:
		stringstream name;		//name of the key
		int compound;			//determines if key is compound (i.e -name=parameter)
		stringstream param; 		//parameter of the key (begins after = sign, ends at the end of the key)

		key();				//cleans name and param, sets compound to zero
		int fill(const char*);		//parses the key and fills name, compound and param fields
		~key();
};

class keyarray
{
	public:
		int keynum;			//number of the keys
		key* arr;			//array of keys
		int memoryexists;		//determines whether memory exists (memoryexists = 1) or not (memoryexists = 0)

		keyarray();			//cleans keynum
		int memory_create(int);		//int memory_create(int newkewnum)
						//if memory does not exist, creates it
						//if exists, recreates it
						//the dimension of (re)created memory is newkeynum
		int find_parse_keys(int, char**);
						//int find_parse_keys(int argc, char** argv)
						//first element of argv is program name
						//then finds all subsequent arguments, which start from '-'
						//(the search stops when first argument, that does not start from '-' is found)
						//then creates enough memory and invokes key::fill method for each key
		int print_keys();		//prints information about keys:
						//total number of keys
						//name, compound, param for each keys 
		int memory_remove();		//removes the memory if it exists and sets keynum and memoryexists to zero
		~keyarray();
};

#endif
