#include "mpi.h"
#include <stdio.h>
#include <stdlib.h>

void write_pgm(FILE *f, unsigned char *v, int sizex, int sizey) 
{
  int i, j;

  fprintf(f, "P2\n");
  fprintf(f, "%d %d\n", sizex, sizey);
  fprintf(f, "%d\n", 255);
  
  for(i=0; i<sizey; i++ ) {
    for(j=0; j<sizex; j++ ) {
      fprintf(f," %d", (unsigned int)v[i*sizex+j]);
    }
    fprintf(f,"\n");
  }
}


#define setpixel(x_,y_)						\
  {								\
    int xx,yy;							\
    xx=(x_)%sizex;						\
    yy=(y_)%sizey;						\
    xx=(xx+sizex)%sizex;					\
    yy=(yy+sizey)%sizey;					\
    v[yy*sizex+xx]=1;						\
  }


void smooth(unsigned char *v, int sizex, int sizey) {
  int x, y;

  for( x=1; x<sizex-1; x++ ) {
    for( y=1; y<sizey-1; y++ ) {
      v[y*sizex+x]=
	( 0.40 * v[y*sizex+x] +
	  0.15 * v[(y-1)*sizex+x]+
	  0.15 * v[(y+1)*sizex+x]+
	  0.15 * v[y*sizex+(x+1)]+
	  0.15 * v[y*sizex+(x-1)] );
    }
  }
  
}


void draw_circle(unsigned char *v, int sizex, int sizey,
		 int x0, int y0, int r) 
{
  int f=1-r;
  int ddF_x=1;
  int ddF_y=-2*r;
  int x=0;
  int y=r;

  setpixel(x0-r, y0);
  setpixel(x0+r, y0);
  setpixel(x0  , y0-r);
  setpixel(x0  , y0+r);

  while(x<y) {
    if(f>=0) {
      y--;
      ddF_y+=2;
      f+=ddF_y;
    }
    x++;
    ddF_x+=2;
    f+=ddF_x;
    setpixel(x0+x, y0+y);
    setpixel(x0-x, y0+y);
    setpixel(x0+x, y0-y);
    setpixel(x0-x, y0-y);
    setpixel(x0+y, y0+x);
    setpixel(x0-y, y0+x);
    setpixel(x0+y, y0-x);
    setpixel(x0-y, y0-x);
  }
}


int main(int argc, char* argv[])
{
  unsigned char* array;
  unsigned char* subarray;
  int i,j,k,l, sizex, sizey, targetrank;
  int niter;
  FILE *outfile;

//MPI initialisation
  int myrank, myvalue, numprocs, P,Q, subsizex, subsizey;
  int* coords;
  P = 2;
  Q = 2;
  int ndims = 2;
  int dims[2] = { P , Q  };
  int periods[2] = { 1,1 };
  int reorder = 1;
  MPI_Status status;
  MPI_Comm cart_comm;
  

  MPI_Init(&argc, &argv);
  MPI_Comm_size(MPI_COMM_WORLD,&numprocs);
  MPI_Comm_rank(MPI_COMM_WORLD,&myrank);
  
  //MPI_Dims_create(P,Q,dims);
  MPI_Cart_create(MPI_COMM_WORLD, ndims,dims,periods,reorder, &cart_comm);

  sizex = 1000;
  sizey = 1000;
  niter = 20;
  subsizex = sizex /P;
  subsizey = sizey /Q;
  subarray=malloc(subsizex*subsizey*sizeof(unsigned char));
  if(myrank==0)
  { 
  	unsigned char* zerosub;
	zerosub=malloc(subsizex*subsizey*sizeof(unsigned char));
  	array=malloc(sizex*sizey*sizeof(unsigned char));
  	for(i=0; i<sizex*sizey; i++) 
    		array[i]=255;
  

  	draw_circle(array, sizex, sizey, 0, 0, 40);
  	draw_circle(array, sizex, sizey, 0, 0, 30);
  	draw_circle(array, sizex, sizey, 100, 100, 10);
  	draw_circle(array, sizex, sizey, 100, 100, 20);
  	draw_circle(array, sizex, sizey, 100, 100, 30);
  	draw_circle(array, sizex, sizey, 100, 100, 40);
  	draw_circle(array, sizex, sizey, 100, 100, 50);
	int up,down,left,right;
	MPI_Cart_shift(cart_comm,0,1,&left,&right);
	MPI_Cart_shift(cart_comm,1,1,&up,&down);
	printf("P:%d My neighbors are r: %d d:%d 1:%d u:%d\n",myrank,right,down,left,up);
		
	coords = malloc(2*sizeof(int));
	for(i=0;i<P;i++)
	{
		for(j=0;j<Q;j++)
		{
			for(k=0;k<subsizex;k++)
			{
				for(l=0;l<subsizey;l++)
				{
					subarray[l*subsizex+k] = array[Q*sizex+l*subsizex+P*subsizey+k];				

				}

			}
		coords[0] = i;
		coords[1] = j;
		MPI_Cart_rank(cart_comm, coords, &targetrank);
		printf("The target rank is %d and the array has size %d\n", targetrank,subsizex*subsizey);
		if(targetrank!=0){
			MPI_Send(subarray,subsizex*subsizey,MPI_UNSIGNED_CHAR,targetrank,123,cart_comm);
		}
		else{
			memcpy(zerosub,subarray,subsizex*subsizey);		
		}
		}

	}
	subarray = zerosub;
	
  }
  else
  {
	MPI_Recv(subarray,subsizex*subsizey,MPI_UNSIGNED_CHAR,0,123,cart_comm,&status);

  }

 // MPI_Bcast(array,sizex*sizey,MPI_UNSIGNED_CHAR,0,MPI_COMM_WORLD);
  
  for( i=0; i<niter; i++ ) {
    smooth(subarray, subsizex, subsizey);
  }
  
  if(myrank==0) {

  	for(i=0;i<P;i++)
	{
		for(j=0;j<Q;j++)
		{
			for(k=0;k<subsizex;k++)
			{
				for(l=0;l<subsizey;l++)
				{
					array[Q*sizex+l*subsizex+P*subsizey+k]	= subarray[l*subsizex+k];				

				}

			}
		coords[0] = i;
		coords[1] = j;
		MPI_Cart_rank(cart_comm, coords, &targetrank);
		printf("The target rank is %d and the array has size %d\n", targetrank,subsizex*subsizey);
		if(targetrank!=0){
			MPI_Recv(subarray,subsizex*subsizey,MPI_UNSIGNED_CHAR,targetrank,123,cart_comm,&status);
		}
		}

	}


 	 outfile = fopen("testimg.pgm", "w");
  	write_pgm(outfile, array, sizex, sizey);

  	fclose(outfile);
  }
  else {

	MPI_Send(subarray,subsizex*subsizey,MPI_UNSIGNED_CHAR,0,123,cart_comm);

  }

  MPI_Finalize();
  return 0;
}
