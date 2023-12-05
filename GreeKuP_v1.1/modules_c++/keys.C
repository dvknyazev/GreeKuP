#ifndef KEYS_CPP
#define KEYS_CPP
#include <iostream>
#include <sstream>
using namespace std;
#include "keys.h"

key::key()
{
	this->name.str("");
	this->compound = 0;
	this->param.str("");
}

int key::fill(const char* argum)
{
	int keyparser;

	if (*argum!='-')
	{	
		return 1;
	}	
	keyparser = 1;
	while ((*(argum + keyparser)!='=')&&(*(argum + keyparser)!='\0'))
	{
		(this->name) << *(argum + keyparser);
		keyparser++;
	}

	if (*(argum + keyparser)=='\0')
	{
		this->compound = 0;
		return 0;
	}
	else
	{
		keyparser++;
		this->compound = 1;
		while (*(argum + keyparser)!='\0')
		{
			(this->param) << *(argum + keyparser);
			keyparser++;
		}

		return 0;
	}
}

key::~key()
{
}

keyarray::keyarray()
{
	this->keynum = 0;
	this->memoryexists = 0;
}

int keyarray::memory_create(int newkeynum)
{
	if (this->memoryexists != 0)
	{
		this->memory_remove();
	}
	this->arr = new key [newkeynum];
	this->keynum = newkeynum;
	this->memoryexists = 1;
	
	return 0;
}

int keyarray::find_parse_keys(int argvsize, char** argvargum)
{
	int keycounter = 0;
	while (((keycounter + 1) < argvsize)&&(*(*(argvargum + keycounter + 1))=='-'))
	{	
		keycounter++;
	}
	this->memory_create(keycounter);
	for (keycounter = 0; keycounter < this->keynum; keycounter++)
	{
		this->arr[keycounter].fill(*(argvargum + 1 + keycounter));
	}
	return 0;
}

int keyarray::print_keys()
{
	int keycounter;
	cout << this->keynum << " keys found" << endl;
	for (keycounter = 0; keycounter < this->keynum; keycounter++)
	{
		cout << endl;
		cout << "Key name: " << this->arr[keycounter].name.str() << endl;
		if (this->arr[keycounter].compound)
		{
			cout << "Compound" << endl;
		}
		else
		{
			cout << "Simple" << endl;
		}
		cout << "Key parameter: " << this->arr[keycounter].param.str() << endl;

	}
	return 0;
}

int keyarray::memory_remove()
{
	if (this->memoryexists != 0)
	{
		delete [] this->arr;
		this->keynum = 0;
		this->memoryexists = 0;
	}
	return 0;
}

keyarray::~keyarray()
{
}

#endif
