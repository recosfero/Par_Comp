#include <stdio.h>
#include <sys/time.h>
#include <time.h>
#include <math.h>

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


#ifndef CHECK
#define CHECK 0
#endif

#define MYTIMEVAL( tv_ )			\
    ((tv_.tv_sec)+(tv_.tv_usec)*1.0e-6)

#define MYTIMESTAMP( time_ )				\
    {							\
	static struct timeval tv;			\
	gettimeofday( &tv, NULL );			\
	time_=MYTIMEVAL(tv);				\
    }

double A[SIZE_M][SIZE_N];
double B[SIZE_N][SIZE_K];
double C[SIZE_M][SIZE_K];

/*Prototypes*/
void executeUnloopMethod(double nflops, void (*unrollMethod) (double (*)[SIZE_N], double (*)[SIZE_N], double (*)[SIZE_N]), char *name);
void mmult(double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K]);
void mmult_unroll4( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K]);
void mmult_unroll8( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K]);
void mmult_unroll16( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K]);
void mmult_unroll24(double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K]);
void mmult_unroll28( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K]);
void mmult_unroll222( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K]);



int main(int argc, char* argv[])
{
    int i, j;
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

    MYTIMESTAMP(tstart);
    cblas_dgemm(CblasRowMajor,  
		CblasNoTrans, CblasNoTrans, SIZE_M, SIZE_N, SIZE_K,  
		1.0, (const double*)A, SIZE_M, (const double*)B, 
		SIZE_N, 0.0, (double*)C, SIZE_K);  
    MYTIMESTAMP(tstop);

    tdgemm = tstop-tstart;

    fprintf(stderr, "#M,N,K,tmmult,tdgemm,gflops_mmult,gflops_dgemm\n");
    fprintf(stderr, "%d,%d,%d,%f,%f,%f,%f\n",
	    SIZE_M, SIZE_N, SIZE_K, tmmult, tdgemm, 
	    1.0e-6*nflop/tmmult,
	    1.0e-6*nflop/tdgemm);
    
    executeUnloopMethod(nflop, mmult_unroll4, "Unloop 4");
    executeUnloopMethod(nflop, mmult_unroll8, "Unloop 8");
    executeUnloopMethod(nflop, mmult_unroll16, "Unloop 26");
    executeUnloopMethod(nflop, mmult_unroll24, "Unloop 24");
    executeUnloopMethod(nflop, mmult_unroll28, "Unloop 28");
    executeUnloopMethod(nflop, mmult_unroll222, "Unloop 222");
    return 0;
}

void compare(double D[SIZE_M][SIZE_N])
{
  int i,j;
  double x = 10e-3;

  double E[SIZE_M][SIZE_N];

  mmult(A,B,E);

  for (i = 0; i < SIZE_M; i++) {
	    for (j = 0; j < SIZE_K; j++) {
          double diff = D[i][j] - E[i][j];
          if (fabs(diff) > x){
            fprintf(stderr,"FAILED:%d, %d \n", i,j);
            return;
           }
          
	    }
    }
    fprintf(stderr,"PASSED\n");
}

void executeUnloopMethod(double nflops, void (*unrollMethod) (double (*)[SIZE_N], double (*)[SIZE_N], double (*)[SIZE_N]), char *name)  {
    double tstart, tstop, tmmult;
    
    MYTIMESTAMP(tstart);
    (*unrollMethod)( A, B, C);
    MYTIMESTAMP(tstop);
    
    tmmult = tstop-tstart;
    
    fprintf(stderr, "%s: %fsec, %fgflops\n", name, tmmult, 1.0e-6*nflops/tmmult);
  
    if(CHECK == 1){
      fprintf(stderr, "%s\t", name);
      compare(C);
    }
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


void mmult_unroll222( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K])
{
    int i, j, k;
    double sum, sum1, sum2, sum3;
  
    for (i = 0; i < SIZE_M; i+=2) {
	    for (j = 0; j < SIZE_K; j+=2) {
	      sum = 0.0;
        sum1 = 0.0;
        sum2 = 0.0;
        sum3 = 0.0;
	      for (k = 0; k < SIZE_N; k+=2) {
		      sum += A[i][k] * B[k][j];
          sum += A[i][k+1] * B[k+1][j];
		      
          sum1 += A[i][k] * B[k][j+1];
          sum1 += A[i][k+1] * B[k+1][j+1];
		      
          sum2 += A[i+1][k] * B[k][j];
          sum2 += A[i+1][k+1] * B[k+1][j];
		      
          sum3 += A[i+1][k] * B[k][j+1];
          sum3 += A[i+1][k+1] * B[k+1][j+1];
	      }
	      
        C[i][j] = sum;
	      C[i][j+1] = sum1;
        C[i+1][j] = sum2;
	      C[i+1][j+1] = sum3;
	    }
    }
}

void mmult_unroll24( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K])
{
    int i, j, k;
    double sum, sum1;
  
    for (i = 0; i < SIZE_M; i++) {
	    for (j = 0; j < SIZE_K; j+=2) {
	      sum = 0.0;
        sum1 = 0.0;
	      for (k = 0; k < SIZE_N; k+=4) {
		      sum += A[i][k] * B[k][j];
          sum += A[i][k+1] * B[k+1][j];
          sum += A[i][k+2] * B[k+2][j];
          sum += A[i][k+3] * B[k+3][j];
		      
          sum1 += A[i][k] * B[k][j+1];
          sum1 += A[i][k+1] * B[k+1][j+1];
          sum1 += A[i][k+2] * B[k+2][j+1];
          sum1 += A[i][k+3] * B[k+3][j+1];
	      }
	      C[i][j] = sum;
	      C[i][j+1] = sum1;
	    }
    }
 
}

void mmult_unroll28( double A[SIZE_M][SIZE_N], double B[SIZE_N][SIZE_K], double C[SIZE_M][SIZE_K])
{
    int i, j, k;
    double sum, sum1;
    
    for (i = 0; i < SIZE_M; i++) {
	    for (j = 0; j < SIZE_K; j+=2) {
            sum = 0.0;
            sum1 = 0.0;
            for (k = 0; k < SIZE_N; k+=8) {
                sum += A[i][k] * B[k][j];
                sum += A[i][k+1] * B[k+1][j];
                sum += A[i][k+2] * B[k+2][j];
                sum += A[i][k+3] * B[k+3][j];
                sum += A[i][k+4] * B[k+4][j];
                sum += A[i][k+5] * B[k+5][j];
                sum += A[i][k+6] * B[k+6][j];
                sum += A[i][k+7] * B[k+7][j];
                sum1 += A[i][k] * B[k][j+1];
                sum1 += A[i][k+1] * B[k+1][j+1];
                sum1 += A[i][k+2] * B[k+2][j+1];
                sum1 += A[i][k+3] * B[k+3][j+1];
                sum1 += A[i][k+4] * B[k+4][j+1];
                sum1 += A[i][k+5] * B[k+5][j+1];
                sum1 += A[i][k+6] * B[k+6][j+1];
                sum1 += A[i][k+7] * B[k+7][j+1];
            }
            C[i][j] = sum;
            C[i][j+1] = sum1;
	    }
    }

     
}

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
            
            for (k = 0; k < SIZE_N - 8; k+=16) {
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

            
            C[i][j] = sum;
	    }
    }
}

