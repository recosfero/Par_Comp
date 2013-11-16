
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


//
// C = A*B
// (M x N) * (N x K) -> (M x K)
// 
void mmult_unroll4( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K])
{
    int i, j, k;
    double sum;
  
    for (i = 0; i < SIZE_M; i++) {
	    for (j = 0; j < SIZE_K; j++) {
	      sum = 0.0;
	      for (k = 0; k < SIZE_N; k+=4) {
		      sum += A[i][k] * B[k][j];
          sum += A[i][k+1] * B[k+1][j];
          sum += A[i][k+2] * B[k+2][j];
          sum += A[i][k+3] * B[k+3][j];
	      }
	      C[i][j] = sum;
	    }
    }
 
}

//
void mmult_unroll8( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K])
{
    int i, j, k;
    double sum;
    
    for (i = 0; i < SIZE_M; i++) {
	    for (j = 0; j < SIZE_K; j++) {
            sum = 0.0;
            for (k = 0; k < SIZE_N; k+=8) {
                sum += A[i][k] * B[k][j];
                sum += A[i][k+1] * B[k+1][j];
                sum += A[i][k+2] * B[k+2][j];
                sum += A[i][k+3] * B[k+3][j];
                sum += A[i][k+4] * B[k+4][j];
                sum += A[i][k+5] * B[k+5][j];
                sum += A[i][k+6] * B[k+6][j];
                sum += A[i][k+7] * B[k+7][j];
            }
            C[i][j] = sum;
	    }
    }
    
}

void mmult_unroll16( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K])
{
    int i, j, k;
    double sum;
    
    for (i = 0; i < SIZE_M; i++) {
	    for (j = 0; j < SIZE_K; j++) {
            sum = 0.0;
            for (k = 0; k < SIZE_N; k+=16) {
                sum += A[i][k] * B[k][j];
                sum += A[i][k+1] * B[k+1][j];
                sum += A[i][k+2] * B[k+2][j];
                sum += A[i][k+3] * B[k+3][j];
                sum += A[i][k+4] * B[k+4][j];
                sum += A[i][k+5] * B[k+5][j];
                sum += A[i][k+6] * B[k+6][j];
                sum += A[i][k+7] * B[k+7][j];
                sum += A[i][k+8] * B[k+8][j];
                sum += A[i][k+9] * B[k+9][j];
                sum += A[i][k+10] * B[k+10][j];
                sum += A[i][k+11] * B[k+11][j];
                sum += A[i][k+12] * B[k+12][j];
                sum += A[i][k+13] * B[k+13][j];
                sum += A[i][k+14] * B[k+14][j];
                sum += A[i][k+15] * B[k+15][j];
            }
            
            //1000 % 16 = 8;

            
            sum += A[i][k] * B[k][j];
            sum += A[i][k+1] * B[k+1][j];
            sum += A[i][k+2] * B[k+2][j];
            sum += A[i][k+3] * B[k+3][j];
            sum += A[i][k+4] * B[k+4][j];
            sum += A[i][k+5] * B[k+5][j];
            sum += A[i][k+6] * B[k+6][j];
            sum += A[i][k+7] * B[k+7][j];
            sum += A[i][k+8] * B[k+8][j];

            
            C[i][j] = sum;
	    }
    }
    
}

void mmult_unroll4( double A[SIZE_M][SIZE_N])


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
    
    mainUnloop4(nflops);
    mainUnloop8(nflops);
    mainUnloop16(nflops);

    return 0;
}

int mainUnloop4(int nflops) {
    double tstart, tstop, tmmult;
    
    
    MYTIMESTAMP(tstart);
    mmult_unroll4( A, B, C);
    MYTIMESTAMP(tstop);
    
    tmmult = tstop-tstart;
    
    printf("Unroll 4: %fsec, %fgflops", tmmult, nflops/tmmult*1.0e-6);
    
    return 0;
}

int mainUnloop8(int nflops) {
    double tstart, tstop, tmmult;
    
    
    MYTIMESTAMP(tstart);
    mmult_unroll8( A, B, C);
    MYTIMESTAMP(tstop);
    
    tmmult = tstop-tstart;
    
    printf("Unroll 8: %fsec, %fgflops", tmmult, nflops/tmmult*1.0e-6);
    
    return 0;
}

int mainUnloop16(int nflops) {
    double tstart, tstop, tmmult;
    
    
    MYTIMESTAMP(tstart);
    mmult_unroll16( A, B, C);
    MYTIMESTAMP(tstop);
    
    tmmult = tstop-tstart;
    
    printf("Unroll 16: %fsec, %fgflops", tmmult, nflops/tmmult*1.0e-6);
    
    return 0;
}
