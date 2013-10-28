
#include <stdio.h>
#include <sys/time.h>
#include <time.h>

#include <mkl_cblas.h>


#ifndef SIZE_M
#define SIZE_M 2000
#endif 

#ifndef SIZE_N 
#define SIZE_N 2000
#endif 

#ifndef SIZE_K
#define SIZE_K 2000
#endif


#define MYTIMEVAL( tv_ )			\
    ((tv_.tv_sec)+(tv_.tv_usec)*1.0e-6)

#define MYTIMESTAMP( time_ )				\
    {							\
	static struct timeval tv;			\
	gettimeofday( &tv, NULL );			\
	time_=MYTIMEVAL(tv);				\
    }

//
// C = A*B
// (M x N) * (N x K) -> (M x K)
// 
void mmult( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K])
{
    int i, j, k;
    double sum;
  
    for (i = 0; i < SIZE_M; i++) {
	for (j = 0; j < SIZE_K; j++) {
	    sum = 0.0;
	    for (k = 0; k < SIZE_N; k++) {
		sum += A[i][k] * B[k][j];
	    }
	    C[i][j] = sum;
	}
    }
 
}


double A[SIZE_M][SIZE_N];
double B[SIZE_N][SIZE_K];
double C[SIZE_M][SIZE_K];

int main(int argc, char* argv[])
{
    int i, j;
    double sum;
    double tstart, tstop;
    double nflop;
    double tmmult, tdgemm;

    for( i=0; i<SIZE_M; i++ ) {
	for( j=0; j<SIZE_N; j++ ) {
	    A[i][j]=(double)(i)+(double)(j);
	}
    }
    for( i=0; i<SIZE_N; i++ ) {
	for( j=0; j<SIZE_K; j++ ) {
	    B[i][j]=(double)(i)+(double)(j);
	}
    }

    nflop = 2.0*(double)SIZE_M*(double)SIZE_N*(double)SIZE_K;

    
    MYTIMESTAMP(tstart);
    mmult( A, B, C);
    MYTIMESTAMP(tstop);
    
    tmmult = tstop-tstart;

/*
    sum=0.0;
    for( i=0; i<M && i<K; i++ ) {
	sum+=C[i][i];
    }
    fprintf(stderr, "Trace mmult: %f\n", sum);
*/  
  
    MYTIMESTAMP(tstart);
    cblas_dgemm(CblasRowMajor,  
		CblasNoTrans, CblasNoTrans, SIZE_M, SIZE_N, SIZE_K,  
		1.0, (const double*)A, SIZE_M, (const double*)B, 
		SIZE_N, 0.0, (double*)C, SIZE_K);  
    MYTIMESTAMP(tstop);

    tdgemm = tstop-tstart;

/*
    sum=0.0;
    for( i=0; i<M && i<K; i++ ) {
	sum+=C[i][i];
    }
    fprintf(stderr, "Trace dgemm: %f\n", sum);
*/

    fprintf(stderr, "#M,N,K,tmmult,tdgemm,gflops_mmult,gflops_dgemm\n");
    fprintf(stderr, "%d,%d,%d,%f,%f,%f,%f\n",
	    SIZE_M, SIZE_N, SIZE_K, tmmult, tdgemm, 
	    1.0e-6*nflop/tmmult,
	    1.0e-6*nflop/tdgemm);

    return 0;
}
