#include "../include/mklrand.hpp"

mkl_drand::mkl_drand(int size, int seed)
{
	arr_size = size;
	curr = 0;
	randarr = (double*)malloc(arr_size*sizeof(double));
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this -> fill();
}

mkl_drand::~mkl_drand()
{
	free(randarr);
}

double mkl_drand::gen()
{
	double out = randarr[curr];
	curr++;
	if(curr >= arr_size)
	{
		this->fill();
	}
	return out;
}

void mkl_drand::fill()
{
	curr = 0;
	vdRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, arr_size, randarr, 0, 1);
}

void mkl_drand::change_seed(int seed)
{
	vslDeleteStream(&stream);
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this->fill();
}

mkl_irand::mkl_irand(int size, int seed)
{
	arr_size = size;
	curr = 0;
	randarr = (int*)malloc(arr_size*sizeof(int));

	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this -> fill();
}

mkl_irand::~mkl_irand()
{
	free(randarr);
}

int mkl_irand::gen()
{
	int out = randarr[curr];
	curr++;
	if(curr >= arr_size)
	{
		this->fill();
	}
	return out;
}

void mkl_irand::fill()
{
	curr = 0;
	viRngUniform(VSL_RNG_METHOD_UNIFORM_STD, stream, arr_size, randarr, 0, 2);
}

void mkl_irand::change_seed(int seed)
{
	vslDeleteStream(&stream);
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this->fill();
}

mkl_lnrand::mkl_lnrand(double m, double sd, int size, int seed)
{
	lmean = m;
	lsd = sd;
	arr_size = size;
	curr = 0;
	randarr = (double*)malloc(arr_size*sizeof(double));
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this -> fill();
}

mkl_lnrand::~mkl_lnrand()
{
	free(randarr);
}

double mkl_lnrand::gen()
{
	double out = randarr[curr];
	curr++;
	if(curr >= arr_size)
	{
		this->fill();
	}
	return out;
}

void mkl_lnrand::fill()
{
	curr = 0;
	vdRngLognormal(VSL_RNG_METHOD_LOGNORMAL_BOXMULLER2, stream, arr_size, randarr, lmean, lsd, 0, 1);
}

void mkl_lnrand::change_seed(int seed)
{
	vslDeleteStream(&stream);
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this->fill();
}

mkl_nrand::mkl_nrand(double m, double sdin, int size, int seed)
{
	mean = m;
	sd = sdin;
	arr_size = size;
	curr = 0;
	randarr = (double*)malloc(arr_size*sizeof(double));
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this -> fill();
}

mkl_nrand::~mkl_nrand()
{
	free(randarr);
}

double mkl_nrand::gen()
{
	double out = randarr[curr];
	curr++;
	if(curr >= arr_size)
	{
		this->fill();
	}
	return out;
}

void mkl_nrand::fill()
{
	curr = 0;
	vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, arr_size, randarr, mean, sd);
}

void mkl_nrand::change_seed(int seed)
{
	vslDeleteStream(&stream);
	vslNewStream(&stream, VSL_BRNG_SFMT19937, seed);
	this->fill();
}
