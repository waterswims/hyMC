#ifndef _MKLRAND
#define _MKLRAND

#include "mkl.h"

class mkl_randbase
{
protected:
	int arr_size;
	int curr;
	VSLStreamStatePtr stream;
public:
	mkl_randbase(){}
	~mkl_randbase(){vslDeleteStream(&stream);}
	virtual void fill(){}
	virtual void change_seed(){}
	void save(const char* name) {vslSaveStreamF(stream, name);}
	void load(const char* name) {vslDeleteStream(&stream); vslLoadStreamF(&stream, name);}
};

class mkl_drand: public mkl_randbase
{
private:
	double *randarr;

public:
	mkl_drand(int size, int seed=1);
	~mkl_drand();
	double gen();
	void fill();
	void change_seed(int seed);
};

class mkl_irand: public mkl_randbase
{
private:
	int *randarr;

public:
	mkl_irand(int size, int seed=1);
	~mkl_irand();
	int gen();
	void fill();
	void change_seed(int seed);
};

class mkl_lnrand: public mkl_randbase
{
private:
	double *randarr;
	double lmean, lsd;

public:
	mkl_lnrand(double m, double sd, int size, int seed=1);
	~mkl_lnrand();
	double gen();
	void fill();
	void change_seed(int seed);
};

#endif
