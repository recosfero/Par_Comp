#include <stdio.h>
#include <sys/time.h>
#include <time.h>

#include <mkl_cblas.h>


#ifndef SIZE_M
#define SIZE_M 1000
#endif 

#ifndef SIZE_N 
#define SIZE_N 1000
#endif 

#ifndef SIZE_K
#define SIZE_K 1000
#endif

#define max(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a > _b ? _a : _b; })

#define min(a,b) \
  ({ __typeof__ (a) _a = (a); \
  __typeof__ (b) _b = (b); \
  _a < _b ? _a : _b; })

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
void mmult( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K],int BLOCK)
{
    int i, j, k, jj, kk;
    
    
    for (jj=0; jj< SIZE_K; jj += BLOCK){
      for (kk=0; kk < SIZE_N ; kk+=BLOCK){
	for (i = 0; i < SIZE_M ; i++) {
	  for (j = jj; j < min(jj+BLOCK-1,SIZE_K); j++) {
	    
	    for (k = kk; k < min(SIZE_N,kk+BLOCK-1); k++) {
		C[k][i] += C[k][i] +  A[j][i] * B[k][j];
	    }
	   
	  }
	}
      }
    }
 
}

double measure_mult(double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K],int BLOCK) {
    double tstart, tstop;
    MYTIMESTAMP(tstart);

    mmult( A, B, C,BLOCK);
    MYTIMESTAMP(tstop);

    return tstop-tstart;
}

double A[SIZE_M][SIZE_N];
double B[SIZE_N][SIZE_K];
double C[SIZE_M][SIZE_K];

int main(int argc, char* argv[])
{
    int i, j, k;
    double sum;
    double tstart, tstop;
    double nflop;
    double tmmult, tdgemm;
    double times[5];
    //int BLOCKSIZE = 8;

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

    
    //MYTIMESTAMP(tstart);
    
    //mmult( A, B, C,BLOCKSIZE);
    //MYTIMESTAMP(tstop);
    k = 0;
    for( i=2;i<=32;i=i*2){
    	tmmult = measure_mult(A,B,C,i);
    	times[k] = tmmult;
	k++;
    }
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

    fprintf(stderr, "#M,N,K,tmmult_2,tmmult_4,tmmult_8,tmmult_16,tmmult_32,tdgemm,gflops_mmult_2,gflops_mmult_4,gflops_mmult_8,gflops_mmult_16,gflops_mmult_32,gflops_dgemm\n");
    fprintf(stderr, "%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f\n",
	    SIZE_M, SIZE_N, SIZE_K, times[0], times[1], times[2], times[3], times[4], tdgemm, 
	    1.0e-6*nflop/times[0],
	    1.0e-6*nflop/times[1],
            1.0e-6*nflop/times[2],
            1.0e-6*nflop/times[3],
            1.0e-6*nflop/times[4],
	    1.0e-6*nflop/tdgemm);

    return 0;
}
