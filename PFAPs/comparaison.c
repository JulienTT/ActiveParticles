#include <stdio.h>
#include <stdlib.h>
#include <string.h>
#include <math.h>
#include <time.h>
#include "mt19937-64.c"
#include "mkl_vsl.h"

// COMPILATION ON TOWER
//   -----------
//   gcc -O3 -DMKL_ILP64 -m64 -I/opt/intel/mkl/include AOUP_current_mkl.c -o RUN_mkl -Wl,--start-group /opt/intel/mkl/lib/intel64/libmkl_intel_ilp64.a /opt/intel/mkl/lib/intel64/libmkl_gnu_thread.a /opt/intel/mkl/lib/intel64/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
//
//   COMPILATION ON CASSOULET :
//   gcc -O3 -DMKL_ILP64 -m64 -I/usr/include/mkl/ AOUP_current_mkl.c -o RUN_mkl -Wl,--start-group /usr/lib/x86_64-linux-gnu/libmkl_intel_ilp64.a /usr/lib/x86_64-linux-gnu/libmkl_gnu_thread.a /usr/lib/x86_64-linux-gnu/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl
//
//   COMPILATION ON CONFIT :
//   gcc -O3 -DMKL_ILP64 -m64 -I/usr/include/mkl/ AOUP_current_mkl.c -o RUN_mkl -Wl,--start-group /usr/lib/x86_64-linux-gnu/libmkl_intel_ilp64.a /usr/lib/x86_64-linux-gnu/libmkl_gnu_thread.a /usr/lib/x86_64-linux-gnu/libmkl_core.a -Wl,--end-group -lgomp -lpthread -lm -ldl


void gauss_stl(long N, double sigma) {
	int i;
	clock_t start, end;
	init_genrand64(2018);
	double *array;
	array = (double*) malloc(sizeof(double)*N);
	start = clock();

	for (i = 0 ; i < N ; i++) {
		array[i] = gasdevMT();
	}

	end = clock();
	double time_taken = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);

	printf("mt19937: %lg sec\n",time_taken);

	free(array);
}

// void pointer_try(long N, double sigma) {
// 	int i;
// 	clock_t start, end;
// 	init_genrand64(2018);
// 	array = (double*) malloc(sizeof(double)*N);
// 	array_ini = (double*) malloc(sizeof(double)*10);
// 	start = clock();
//
// 	for (i = 0 ; i < N ; i++) {
// 		array[i] = gasdevMT();
// 	}
// 	array = array_ini
//
// 	end = clock();
// 	double time_taken = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);
//
// 	printf("mt19937: %lg sec\n",time_taken)
//
// 	free(array);
// }

void gauss_mkl1(const size_t N, const double sigma, const size_t batch_size) {
	VSLStreamStatePtr stream;
	clock_t start, end;
	vslNewStream(&stream, VSL_BRNG_MT19937, 2018);

	const size_t n_batch = N / batch_size;
	double *array;
	array = (double*) malloc(sizeof(double)*batch_size*n_batch);

	start = clock();

	for ( size_t i = 0 ; i < n_batch ; i++) {
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, batch_size,
					  &array[i * batch_size], 0.0, sigma);
	}

	end = clock();
	double time_taken = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);

	printf("MKL ICDF: %lg sec\n",time_taken);

	free(array);
	vslDeleteStream(&stream);
}

void gauss_mkl2(const size_t N, const double sigma, const size_t batch_size) {
	VSLStreamStatePtr stream;
	clock_t start, end;
	vslNewStream(&stream, VSL_BRNG_MT19937, 2018);

	const size_t n_batch = N / batch_size;
	double *array;
	array = (double*) malloc(sizeof(double)*batch_size*n_batch);

	start = clock();

	for ( size_t i = 0 ; i < n_batch ; i++) {
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER, stream, batch_size,
					  &array[i * batch_size], 0.0, sigma);
	}

	end = clock();
	double time_taken = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);

	printf("MKL Box-Muller: %lg sec\n",time_taken);

	free(array);
	vslDeleteStream(&stream);
}

void gauss_mkl3(const size_t N, const double sigma, const size_t batch_size) {
	VSLStreamStatePtr stream;
	clock_t start, end;
	vslNewStream(&stream, VSL_BRNG_MT19937, 2018);

	const size_t n_batch = N / batch_size;
	double *array;
	array = (double*) malloc(sizeof(double)*batch_size*n_batch);

	start = clock();

	for ( size_t i = 0 ; i < n_batch ; i++) {
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, batch_size,
					  &array[i * batch_size], 0.0, sigma);
	}

	end = clock();
	double time_taken = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);

	printf("MKL Box-Muller 2: %lg sec\n",time_taken);

	free(array);
	vslDeleteStream(&stream);
}

void gauss_mkl4(const size_t N, const double sigma, const size_t batch_size) {
	VSLStreamStatePtr stream;
	clock_t start, end;
	vslNewStream(&stream, VSL_BRNG_SFMT19937, 2018);

	const size_t n_batch = N / batch_size;
	double *array;
	array = (double*) malloc(sizeof(double)*batch_size*n_batch);

	start = clock();

	for ( size_t i = 0 ; i < n_batch ; i++) {
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_ICDF, stream, batch_size,
					  &array[i * batch_size], 0.0, sigma);
	}

	end = clock();
	double time_taken = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);

	printf("MKL SFMT19937 ICDF: %lg sec\n",time_taken);

	free(array);
	vslDeleteStream(&stream);
}

void gauss_mkl5(const size_t N, const double sigma, const size_t batch_size) {
	VSLStreamStatePtr stream;
	clock_t start, end;
	vslNewStream(&stream, VSL_BRNG_SFMT19937, 2018);

	const size_t n_batch = N / batch_size;
	double *array;
	array = (double*) malloc(sizeof(double)*batch_size*n_batch);

	start = clock();

	for ( size_t i = 0 ; i < n_batch ; i++) {
		vdRngGaussian(VSL_RNG_METHOD_GAUSSIAN_BOXMULLER2, stream, batch_size,
					  &array[i * batch_size], 0.0, sigma);
	}

	end = clock();
	double time_taken = ((double)(end - start)) / ((double)CLOCKS_PER_SEC);

	printf("MKL SFMT19937 Box-Muller 2: %lg sec\n",time_taken);
	// printf("essai de nombre: %lg sec\n",array[10]);

	free(array);
	vslDeleteStream(&stream);
}

int main() {
	size_t N = 100000000;
	double sigma = 1.0;
	int batch_size = 1000;

	gauss_stl(N, sigma);
	// gauss_gsl1(N, sigma);
	// gauss_gsl2(N, sigma);
	// gauss_gsl3(N, sigma);
	gauss_mkl1(N, sigma, batch_size);
	gauss_mkl2(N, sigma, batch_size);
	gauss_mkl3(N, sigma, batch_size);
	gauss_mkl4(N, sigma, batch_size);
	gauss_mkl5(N, sigma, batch_size);

	return 0;
}
