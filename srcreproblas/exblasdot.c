#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <vector>
#include "exblas/exdot.h"


double dotproduct(double *a, double *b, int nrowsI)
{
  //locals
  int i = 0;
  double local_dot = 0.;
  //*c = 0.; //hard reset
  
  //perform local dot product
  for( i = 0; i < nrowsI; i++)
    local_dot += (a[i] * b[i]);

  //completed successfully!
  return local_dot;

}


int main(int argc,char **argv)
{
    int n;
	double *a=NULL,*b=NULL,resdot;


	n=1000;
	a=(double *)malloc(n*sizeof(double));
    b=(double *)malloc(n*sizeof(double));
    FILE *fpRead=fopen("res0.txt","r");
	if(fpRead==NULL)
	{
		return 0;
	}
    for(i = 0; i < n; i++)
    {
        fscanf(fpRead,"%lf ",&a[i]);
    }
	fpRead=fopen("res00.txt","r");
	if(fpRead==NULL)
	{
		return 0;
	}
    for(i = 0; i < n; i++)
    {
        fscanf(fpRead,"%lf ",&b[i]);
    }
	
	// ReproAllReduce -- Begin
    // 因为前面是inf and -nan 所以这里需要2*
	
    //std::vector<int64_t> h_superacc(2*exblas::BIN_COUNT);
    //算的出来可以用下面
    std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    exblas::cpu::exdot (n, a, b, &h_superacc[0]);
    int imin=exblas::IMIN, imax=exblas::IMAX;
    exblas::cpu::Normalize(&h_superacc[0], imin, imax);
	/*
    if (myId == 0) 
	{
        MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    } 
	else 
	{
        MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }
	*/
	/*
    if (myId == 0) 
	{
        beta = exblas::cpu::Round( &h_superacc[0] );
    }
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	*/
	rou1 = exblas::cpu::Round( &h_superacc[0] );
    // ReproAllReduce -- End
	
	//mpfrDot(dimL, sol2L, sol2L, &beta, commun);
	//dotproduct(res0, res00, &rou1, n_dist, commun); //dot product
	/*
    if(myId==0)
		printf ("dot result: %.10e\n", rou1);
	*/
	
	//resdot=dotproduct(a, b, n);
	printf("recursive dot:%.10e\n",rou1);

    return 0;
}
