#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <string.h>
#include <math.h>
#include <mkl_blas.h>
#include <mpi.h>
#include <hb_io.h>
#include <vector>

#include "reloj.h"
#include "ScalarVectors.h"
#include "SparseProduct.h"
#include "ToolsMPI.h"
#include "matrix.h"
#include "common.h"

#include "exblas/exdot.h"

// ================================================================================

#define DIRECT_ERROR 1
#define PRECOND 1
#define VECTOR_OUTPUT 1


int dotproduct(double *a, double *b, double *c, int nrowsI, MPI_Comm *commun)
{
  //locals
  int i = 0;
  double local_dot = 0.;
  //*c = 0.; //hard reset
  
  //perform local dot product
  for( i = 0; i < nrowsI; i++)
    local_dot += (a[i] * b[i]);

  //reduceall to single values
  /* int MPI_Allreduce( */
  /*   void *sendbuf, */
  /*   void *recvbuf, */
  /*   int count, */
  /*   MPI_Datatype datatype, */
  /*   MPI_Op op, */
  /*   MPI_Comm comm */
  /* ); */

  MPI_Allreduce(
		&local_dot,
		c,
		1,
		MPI_DOUBLE,
		MPI_SUM,
		(*commun));

  //completed successfully!
  return 0;

}


void ConjugateGradient (SparseMatrix mat, double *x, double *b, int *sizes, int *dspls, double *rbuf, int myId) 
{
	//size是子矩阵的列维度（也就是原矩阵的维度），sizeR是子矩阵的行维度
    int size = mat.dim2, sizeR = mat.dim1; 
    int IONE = 1; 
    double DONE = 1.0, DMONE = -1.0, DZERO = 0.0;
    int n, n_dist, iter, maxiter, nProcs;
    double beta, tol, rho, alpha1, umbral;
    double *res = NULL, *z = NULL, *d = NULL, *y = NULL;
    double *aux = NULL;
    double t1, t2, t3, t4;
	
	double *r_tld=NULL;
	double *p1=NULL;
	double *p_hat=NULL;
	double *p0=NULL;
	double *v1=NULL;
	double *v0=NULL;
	double rv;
	double *s_hat=NULL;
	double *t=NULL;
	double rho_1;
	double sn;
	double alphan;
	double omegan;
	double *s=NULL;
	double ts,tt;
	double *res00=NULL;
	double *res0=NULL;
	double rou1;
	double rou0;
	double *x0=NULL;
	double *x1=NULL;
	double omega,alpha0;
	double omega1;
	double nb;
#if PRECOND
    int i, *posd = NULL;
    double *diags = NULL;
#endif
    
	MPI_Comm *commun = (MPI_Comm *)malloc(1*sizeof(MPI_Comm));
    commun[0] = MPI_COMM_WORLD;
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs);
    n = size; 
	n_dist = sizeR; 
	maxiter = 16 * size; 
	//maxiter=16000;
	umbral = 1.0e-8;
    CreateDoubles (&res, n_dist); 
	CreateDoubles (&z, n_dist); 
    CreateDoubles (&d, n_dist);

	CreateDoubles (&r_tld,n_dist);
	CreateDoubles (&p1,n_dist);
	CreateDoubles (&p_hat,n_dist);
	CreateDoubles (&p0,n_dist);
	CreateDoubles (&v1,n_dist);
	CreateDoubles(&v0,n_dist);
	CreateDoubles (&s_hat,n_dist);
	CreateDoubles (&t,n_dist);
	CreateDoubles (&s,n_dist);
	CreateDoubles (&res00,n_dist);
	CreateDoubles (&res0,n_dist);
	CreateDoubles (&x0,n_dist);
	CreateDoubles (&x1,n_dist);
#ifdef DIRECT_ERROR
    // init exact solution
    double *res_err = NULL, *x_exact = NULL;
	CreateDoubles (&x_exact, n_dist);
	CreateDoubles (&res_err, n_dist);
	//初始精确解为1
    InitDoubles(x_exact, n_dist, DONE, DZERO);
#endif // DIRECT_ERROR 

#if PRECOND
    CreateDoubles (&y, n_dist);
    CreateInts (&posd, n_dist);
    CreateDoubles (&diags, n_dist);
    GetDiagonalSparseMatrix2 (mat, dspls[myId], diags, posd);
#pragma omp parallel for
    for (i=0; i<n_dist; i++)
        diags[i] = DONE / diags[i];
#endif
    CreateDoubles (&aux, n); 

#if VECTOR_OUTPUT
    // write to file for testing purpose
    FILE *fp;
    if (myId == 0) {
        char name[50];
        sprintf(name, "%d.txt", nProcs);
        fp = fopen(name,"w");
    }
#endif

	
	iter=0;
	int flag=0;
	double beom;
	
	
  /*
  bnrm2 = norm( b );
  if  ( bnrm2 == 0.0 ), bnrm2 = 1.0; end

  r = b - A*x;
  error = norm( r ) / bnrm2;
  if ( error < tol ) return, end
  */
    MPI_Allgatherv (x, sizeR, MPI_DOUBLE, aux, sizes, dspls, MPI_DOUBLE, MPI_COMM_WORLD);
	InitDoubles (z, sizeR, DZERO, DZERO);
    ProdSparseMatrixVectorByRows (mat, 0, aux, z);            			// z = A * x   
    dcopy (&n_dist, b, &IONE, res0, &IONE);                          		// res = b
    daxpy (&n_dist, &DMONE, z, &IONE, res0, &IONE);                      // res -= z   (res=b-Ax)	
	dcopy (&n_dist, res0, &IONE, res00, &IONE);                          		// r00 = r0
	dcopy (&n_dist, res0, &IONE, res, &IONE);                          		// r1 = r0
/*    
    if(myId==0)
    {
	for(int i=0;i<n_dist;i++)
	{
	    printf("%lf\n",res[i]);
	}	
    }
*/	


	
	double vAux[2];
    std::vector<int64_t> h_superacc(2 * exblas::BIN_COUNT);
    std::vector<int64_t> h_superacc_tol(exblas::BIN_COUNT);
	int imin=exblas::IMIN, imax=exblas::IMAX;
	
		exblas::cpu::exdot (n_dist, res, res, &h_superacc[0]);     //r'*r
		// ReproAllReduce -- Begin
		imin=exblas::IMIN, imax=exblas::IMAX;
		exblas::cpu::Normalize(&h_superacc[0], imin, imax);
		if (myId == 0) {
			MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		} else {
			MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (myId == 0) {
			tol = exblas::cpu::Round( &h_superacc[0] );
		}
		MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// ReproAllReduce -- End
		tol=sqrt(tol);
	
	/*
     dotproduct(res, res, &tol, n_dist, commun); //dot product
     tol = sqrt(tol); //root
	 */
	 /*
     dotproduct(b, b, &nb, n_dist, commun); //dot product
     nb = sqrt(nb); //root
	 */
	 
#ifdef DIRECT_ERROR
    // compute direct error
    double direct_err;
	dcopy (&n_dist, x_exact, &IONE, res_err, &IONE);                        // res_err = x_exact
	daxpy (&n_dist, &DMONE, x, &IONE, res_err, &IONE);                      // res_err -= x

    // compute inf norm
    direct_err = norm_inf(n_dist, res_err);
    MPI_Allreduce(MPI_IN_PLACE, &direct_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

//    // compute euclidean norm
//    // direct_err = res_err' * res_err
//    exblas::cpu::exdot (n_dist, res_err, res_err, &h_superacc[0]);
//    // ReproAllReduce -- Begin
//    imin=exblas::IMIN, imax=exblas::IMAX;
//    exblas::cpu::Normalize(&h_superacc[0], imin, imax);
//    if (myId == 0) {
//        MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//    } else {
//        MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//    }
//    if (myId == 0) {
//        direct_err = exblas::cpu::Round( &h_superacc[0] );
//    }
//    MPI_Bcast(&direct_err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//    // ReproAllReduce -- End
//
//    direct_err = sqrt(direct_err);
#endif // DIRECT_ERROR
		
	
	MPI_Barrier(MPI_COMM_WORLD);
    if (myId == 0) 
        reloj (&t1, &t2);
	while((tol>umbral) && (iter<maxiter))
	//while(((tol/nb)>umbral) && (iter<maxiter))
	{
		
	if (myId == 0) 
#ifdef DIRECT_ERROR
            //printf ("%d \t %a \t %a \n", iter, tol, direct_err);
            printf ("%d \t %20.10e \t %20.10e \n", iter, tol, direct_err);
#else        
            printf ("%d \t %20.10e \n", iter, tol);
#endif // DIRECT_ERROR
/*
if(iter==5)
{
#if VECTOR_OUTPUT
//#if 1
    // print aux
    MPI_Allgatherv (x, n_dist, MPI_DOUBLE, aux, sizes, dspls, MPI_DOUBLE, MPI_COMM_WORLD);
    if (myId == 0) {
        //fprintf(fp, "%d ", iter);
        fprintf(fp, "%d\n", iter);
        for (int ip = 0; ip < n; ip++)
            //fprintf(fp, "%20.10e ", aux[ip]);
            fprintf(fp, "%20.10e\n", aux[ip]);
        fprintf(fp, "\n");

        fclose(fp);
    }
#endif
}
*/		
		
		exblas::cpu::exdot (n_dist, res0, res00, &h_superacc[0]);     //r_tld'*r
		// ReproAllReduce -- Begin
		imin=exblas::IMIN, imax=exblas::IMAX;
		exblas::cpu::Normalize(&h_superacc[0], imin, imax);
		if (myId == 0) {
			MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		} else {
			MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (myId == 0) {
			rou1 = exblas::cpu::Round( &h_superacc[0] );
		}
		MPI_Bcast(&rou1, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// ReproAllReduce -- End
		
	    //dotproduct(res0, res00, &rou1, n_dist, commun); //dot product
		
		/*
		if(myId==0)
		{
			printf("%lf\n",rho);
		}
		*/
		if(rou1==0)
		//if(rou1<1e-16)
		{
			//x=x1;
			dcopy (&n_dist, x1, &IONE, x, &IONE);            //x=x1;
			printf("rou1=0,iterative failed!\n");
			break;
		}
		
		if(iter==0)
		{
			dcopy (&n_dist, res0, &IONE, p1, &IONE);            //p=r;	
		}		
		else
		{
			//omega is w0
			beta=(rou1/rou0)*(alpha0/omega);
			//p = r + beta*( p - omega*v );       //mkl去写这个
			beom = -beta*omega;
		//	dscal (&n_dist, &beta, p, &IONE);                                // p = beta * p
			dcopy(&n_dist,res0,&IONE,p1,&IONE);
			daxpy (&n_dist, &beta, p0, &IONE, p1, &IONE);                       // p += r
			daxpy (&n_dist, &beom, v0, &IONE, p1, &IONE);                      // p -= beta*omega * v
		}
/*
		else
		{
			dcopy (&n_dist, res0, &IONE, p1, &IONE);            //p=r;
		}
*/		
		VvecDoubles (DONE, diags, p1, DZERO, p_hat, n_dist);                    // p_hat = D^-1 * p1 
		MPI_Allgatherv (p_hat, sizeR, MPI_DOUBLE, aux, sizes, dspls, MPI_DOUBLE, MPI_COMM_WORLD);
		InitDoubles (v1, sizeR, DZERO, DZERO);
		ProdSparseMatrixVectorByRows (mat, 0, aux, v1);            			// v = A * p1 
		
		exblas::cpu::exdot (n_dist, res00, v1, &h_superacc[0]);     //r00'*v1
		// ReproAllReduce -- Begin
		imin=exblas::IMIN, imax=exblas::IMAX;
		exblas::cpu::Normalize(&h_superacc[0], imin, imax);
		if (myId == 0) {
			MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		} else {
			MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (myId == 0) {
			rv = exblas::cpu::Round( &h_superacc[0] );
		}
		MPI_Bcast(&rv, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// ReproAllReduce -- End
		
		//dotproduct(v1, res00, &rv, n_dist, commun); //dot product
		alpha1=rou1/rv;

		//s = r - alpha*v;
		alphan=-1.0*alpha1;
		dcopy (&n_dist, res0, &IONE, s, &IONE);            //s=r;
		daxpy (&n_dist, &alphan, v1, &IONE, s, &IONE);                      	// s = r - alpha*v 

/*		
		exblas::cpu::exdot (n_dist, s, s, &h_superacc[0]);     //s'*s
		// ReproAllReduce -- Begin
		imin=exblas::IMIN, imax=exblas::IMAX;
		exblas::cpu::Normalize(&h_superacc[0], imin, imax);
		if (myId == 0) 
		{
			MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		else 
		{
			MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (myId == 0) 
		{
			//我感觉这个Round有问题，舍入到0去了!
			sn = exblas::cpu::Round( &h_superacc[0] );
		}
		MPI_Bcast(&sn, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);

		sn = sqrt (sn);

		
		
		//dotproduct(s, s, &sn, n_dist, commun); //dot product
        //sn = sqrt(sn); //root
		


		if(sn<=umbral)
		{
			//x1=x0+alpha1*p1;
            //x=x1;
            		if(myId==0)
			{
            			printf("sn=%lf  norm s < umbral!\n",sn);
			}
			dcopy (&n_dist, x0, &IONE, x1, &IONE);            //x1=x0;
			daxpy (&n_dist, &alpha1, p_hat, &IONE, x1, &IONE);            	// x1=x0+alpha1*p1;
			dcopy (&n_dist, x1, &IONE, x, &IONE);            //x=x1;
			
			//daxpy (&n_dist, &alpha, p_hat, &IONE, x, &IONE);            	// x += alpha * p_hat;
			break;
		}
*/

		/*
		//得有一个归约算norm(s)
		if ( norm(s) < tol ),                          % early convergence check
			x = x + alpha*p_hat;
			resid = norm( s ) / bnrm2;
			break;
		end
		*/
		
		VvecDoubles (DONE, diags, s, DZERO, s_hat, n_dist);                    // s_hat = D^-1 * s
		MPI_Allgatherv (s_hat, sizeR, MPI_DOUBLE, aux, sizes, dspls, MPI_DOUBLE, MPI_COMM_WORLD);
		InitDoubles (t, sizeR, DZERO, DZERO);
		ProdSparseMatrixVectorByRows (mat, 0, aux, t);            			// t = A * s_hat 

		
		// ReproAllReduce -- Begin
        // ts = t' * s 
        exblas::cpu::exdot (n_dist, t, s, &h_superacc[0]);
        imin=exblas::IMIN, imax=exblas::IMAX;
        exblas::cpu::Normalize(&h_superacc[0], imin, imax);

        // compute tolerance
        //     tt = t' * t
        exblas::cpu::exdot (n_dist, t, t, &h_superacc_tol[0]);
        imin=exblas::IMIN, imax=exblas::IMAX;
        exblas::cpu::Normalize(&h_superacc_tol[0], imin, imax);

        // merge two superaccs into one for reduction
        for (int i = 0; i < exblas::BIN_COUNT; i++) 
		{
            h_superacc[exblas::BIN_COUNT + i] = h_superacc_tol[i]; 
        } 

        if (myId == 0) 
		{
            MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], 2*exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        } 
		else 
		{
            MPI_Reduce (&h_superacc[0], NULL, 2*exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
        }
        if (myId == 0) 
		{
            // split them back
            for (int i = 0; i < exblas::BIN_COUNT; i++) 
			{
                h_superacc_tol[i] = h_superacc[exblas::BIN_COUNT + i]; 
            } 
            vAux[0] = exblas::cpu::Round( &h_superacc[0] );
            vAux[1] = exblas::cpu::Round( &h_superacc_tol[0] );
        }
        MPI_Bcast(vAux, 2, MPI_DOUBLE, 0, MPI_COMM_WORLD);
	    ts = vAux[0]; 
	    tt  = vAux[1]; 
        // ReproAllReduce -- End
	
/*	
		exblas::cpu::exdot (n_dist, t, s, &h_superacc[0]);     //t'*s
		// ReproAllReduce -- Begin
		imin=exblas::IMIN, imax=exblas::IMAX;
		exblas::cpu::Normalize(&h_superacc[0], imin, imax);
		if (myId == 0) {
			MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		} else {
			MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (myId == 0) {
			ts = exblas::cpu::Round( &h_superacc[0] );
		}
		MPI_Bcast(&ts, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		
		exblas::cpu::exdot (n_dist, t, t, &h_superacc[0]);     //t'*t
		// ReproAllReduce -- Begin
		imin=exblas::IMIN, imax=exblas::IMAX;
		exblas::cpu::Normalize(&h_superacc[0], imin, imax);
		if (myId == 0) {
			MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		} else {
			MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (myId == 0) {
			tt = exblas::cpu::Round( &h_superacc[0] );
		}
		MPI_Bcast(&tt, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
*/
		/*
		dotproduct(t, s, &ts, n_dist, commun); //dot product
		dotproduct(t, t, &tt, n_dist, commun); //dot product
		*/
		//omega1 is w1
		omega1 = ts/tt;
		/*
		if(myId==0)
		{
			printf("%lf\n",omega);	
		}
		*/
		//x1=x0+alpha1*p1+w1*s;
		dcopy (&n_dist, x0, &IONE, x1, &IONE);            //x1=x0;
		daxpy (&n_dist, &alpha1, p_hat, &IONE, x1, &IONE);
		daxpy (&n_dist, &omega1, s_hat, &IONE, x1, &IONE);
		//daxpy (&n_dist, &alpha, p_hat, &IONE, x, &IONE);                      	// x = x + alpha*p_hat 
		//daxpy (&n_dist, &omega, s_hat, &IONE, x, &IONE);                      	// x = x + alpha*p_hat + omega*s_hat
/* 
    if(myId==0)
    {
	for(int i=0;i<n_dist;i++)
	{
	    printf("%lf\n",x[i]);
	}	
    }
*/		

		//r = s - omega*t;
		omegan=-1.0*omega1;
		dcopy (&n_dist, s, &IONE, res, &IONE);            //r=s;
		daxpy (&n_dist, &omegan, t, &IONE, res, &IONE);                      	// r = s - omega*t 
		
		
		exblas::cpu::exdot (n_dist, res, res, &h_superacc[0]);     //r'*r
		// ReproAllReduce -- Begin
		imin=exblas::IMIN, imax=exblas::IMAX;
		exblas::cpu::Normalize(&h_superacc[0], imin, imax);
		if (myId == 0) {
			MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		} else {
			MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
		}
		if (myId == 0) {
			tol = exblas::cpu::Round( &h_superacc[0] );
		}
		MPI_Bcast(&tol, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
		// ReproAllReduce -- End
		tol = sqrt (tol);
		
		/*
		dotproduct(res, res, &tol, n_dist, commun); //dot product
        tol = sqrt(tol); //root
		*/
		
#ifdef DIRECT_ERROR
        // compute direct error
        dcopy (&n_dist, x_exact, &IONE, res_err, &IONE);                        // res_err = x_exact
        daxpy (&n_dist, &DMONE, x, &IONE, res_err, &IONE);                      // res_err -= x

        // compute inf norm
        direct_err = norm_inf(n_dist, res_err);
        MPI_Allreduce(MPI_IN_PLACE, &direct_err, 1, MPI_DOUBLE, MPI_MAX, MPI_COMM_WORLD);

//        // compute euclidean norm
//        // direct_err = res_err' * res_err
//        exblas::cpu::exdot (n_dist, res_err, res_err, &h_superacc[0]);
//        // ReproAllReduce -- Begin
//        imin=exblas::IMIN, imax=exblas::IMAX;
//        exblas::cpu::Normalize(&h_superacc[0], imin, imax);
//        if (myId == 0) {
//            MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//        } else {
//            MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
//        }
//        if (myId == 0) {
//            direct_err = exblas::cpu::Round( &h_superacc[0] );
//        }
//        MPI_Bcast(&direct_err, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
//        // ReproAllReduce -- End
//
//        direct_err = sqrt(direct_err);
#endif // DIRECT_ERROR
		
/*
    if(myId==0)
    {
	    printf("%lf\n",tol);	
    }
*/		
		rou0=rou1;
		alpha0=alpha1;
		omega=omega1;
		dcopy (&n_dist, res, &IONE, res0, &IONE);            //r0=r1;
		dcopy (&n_dist, p1, &IONE, p0, &IONE);            //p0=p1;
		dcopy (&n_dist, v1, &IONE, v0, &IONE);            //v0=v1;
		dcopy (&n_dist, x1, &IONE, x0, &IONE);            //x0=x1;
		dcopy (&n_dist, x1, &IONE, x, &IONE);            //x=x1;
		iter++;
		if(omega1==0)
		{
			dcopy (&n_dist, x1, &IONE, x, &IONE);            //x=x1;
			printf("w=0, back\n");
			break;
		}
	}
	
	MPI_Barrier(MPI_COMM_WORLD);
    if (myId == 0) 
        reloj (&t3, &t4);
	

	
	
	


#if VECTOR_OUTPUT
//#if 1
    // print aux
    MPI_Allgatherv (x, n_dist, MPI_DOUBLE, aux, sizes, dspls, MPI_DOUBLE, MPI_COMM_WORLD);
    if (myId == 0) {
        //fprintf(fp, "%d ", iter);
        fprintf(fp, "%d\n", iter);
        for (int ip = 0; ip < n; ip++)
            //fprintf(fp, "%20.10e ", aux[ip]);
            fprintf(fp, "%20.10e\n", aux[ip]);
        fprintf(fp, "\n");

        fclose(fp);
    }
#endif

    if (myId == 0) 
	{
        printf ("Size: %d \n", n);
        printf ("Iter: %d \n", iter);
        //printf ("Tol: %a \n", tol);
	printf("flag: %d \n",flag);
        printf ("Tol: %20.10e \n", tol);
        printf ("Time_loop: %20.10e\n", (t3-t1));
        printf ("Time_iter: %20.10e\n", (t3-t1)/iter);
    }

    RemoveDoubles (&aux); RemoveDoubles (&res); RemoveDoubles (&z); RemoveDoubles (&d);
	RemoveDoubles (&res0);
	RemoveDoubles (&res00);
	
	RemoveDoubles (&p0);
	RemoveDoubles (&p1);
	RemoveDoubles (&p_hat);
	RemoveDoubles (&v1);
	RemoveDoubles (&v0);
	RemoveDoubles (&s);
	RemoveDoubles (&s_hat);
	RemoveDoubles (&t);
	RemoveDoubles (&x0);
	RemoveDoubles (&x1);

#if PRECOND
    RemoveDoubles (&diags); RemoveInts (&posd); RemoveDoubles(&y);
#endif
}

/*********************************************************************************/

int main (int argc, char **argv) {
    int dim; 
    double *vec = NULL, *sol1 = NULL, *sol2 = NULL;
    int index = 0, indexL = 0;
    SparseMatrix mat  = {0, 0, NULL, NULL, NULL}, sym = {0, 0, NULL, NULL, NULL};

    int root = 0, myId, nProcs;
    int dimL, dspL, *vdimL = NULL, *vdspL = NULL;
    SparseMatrix matL = {0, 0, NULL, NULL, NULL};
    double *vecL = NULL, *sol1L = NULL, *sol2L = NULL, *rbuf = NULL;
    int mat_from_file, nodes, size_param, stencil_points;
    if (argc == 3) 
	{
        mat_from_file = atoi(argv[2]);
    } 
	else 
	{
        mat_from_file = atoi(argv[2]);
        nodes = atoi(argv[3]);
        size_param = atoi(argv[4]);
        stencil_points = atoi(argv[5]);
    }

    /***************************************/

    MPI_Init (&argc, &argv);

    // Definition of the variables nProcs and myId
    MPI_Comm_size(MPI_COMM_WORLD, &nProcs); MPI_Comm_rank(MPI_COMM_WORLD, &myId);
    root = nProcs-1;
    root = 0;

    /***************************************/

    CreateInts (&vdimL, nProcs); CreateInts (&vdspL, nProcs); 
    if(mat_from_file) 
	{
        if (myId == root) 
		{
            // Creating the matrix
            ReadMatrixHB (argv[1], &sym);
            DesymmetrizeSparseMatrices (sym, 0, &mat, 0);
            dim = mat.dim1;
        }

        // Distributing the matrix
        dim = DistributeMatrix (mat, index, &matL, indexL, vdimL, vdspL, root, MPI_COMM_WORLD);
        dimL = vdimL[myId]; dspL = vdspL[myId];
    }
	//应该是这一段生成病态矩阵
    else 
	{
		//这个变量是原矩阵维度
        dim = size_param * size_param * size_param;
        int divL, rstL, i;
        divL = (dim / nProcs); 
		rstL = (dim % nProcs);
		//每个进程分配的维数，前面的进程数会多1个空间，大于rst的进程数少一个
        for (i=0; i<nProcs; i++) 
			vdimL[i] = divL + (i < rstL);
        vdspL[0] = 0; 
		//每个进程分配的维数的起始位置
		for (i=1; i<nProcs; i++) 
			vdspL[i] = vdspL[i-1] + vdimL[i-1];
		//这两个变量是每个进程的维数及在原矩阵的起始位置
        dimL = vdimL[myId]; 
		dspL = vdspL[myId];
		//这个变量是什么含义？
        int band_width = size_param * (size_param + 1) + 1;
        band_width = 100 * nodes;
		//这个变量是什么含义？非零元个数？
        long nnz_here = ((long) (stencil_points + 2 * band_width)) * dimL;
        printf ("dimL: %d, nodes: %d, size_param: %d, band_width: %d, stencil_points: %d, nnz_here: %ld\n",
                dimL, nodes, size_param, band_width, stencil_points, nnz_here);
		//分配每个进程稀疏子矩阵matL的空间
        allocate_matrix(dimL, dim, nnz_here, &matL);
		//把矩阵的值生成并填充进去
        generate_Poisson3D_filled(&matL, size_param, stencil_points, band_width, dspL, dimL, dim);

        // To generate ill-conditioned matrices
        //double factor = 1.0e6;
        //这个可以
        //double factor = 1.0;
        //double factor = 1.0e1;
        double factor = 1.0e2;
        //目前只能到1e2
        //double factor = 1.0e3;
        ScaleFirstRowCol(matL, dspL, dimL, myId, root, factor);
    }
    MPI_Barrier(MPI_COMM_WORLD);

    // Creating the vectors
    if (myId == root) 
	{
        CreateDoubles (&vec , dim);
        CreateDoubles (&sol1, dim);
        CreateDoubles (&sol2, dim);
        CreateDoubles (&rbuf, nProcs);
        InitRandDoubles (vec, dim, -1.0, 1.0);
        InitDoubles (sol1, dim, 0.0, 0.0);
        InitDoubles (sol2, dim, 0.0, 0.0);
        InitDoubles (rbuf , nProcs, 0.0, 0.0);
    } 
	else 
	{
        CreateDoubles (&vec , dim);
        CreateDoubles (&sol2, dim);
        InitDoubles (vec , dim, 0.0, 0.0);
        InitDoubles (sol2, dim, 0.0, 0.0);
    }
    CreateDoubles (&vecL , dimL);
    CreateDoubles (&sol1L, dimL);
    CreateDoubles (&sol2L, dimL);
    InitDoubles (vecL , dimL, 0.0, 0.0);
    InitDoubles (sol1L, dimL, 0.0, 0.0);
    InitDoubles (sol2L, dimL, 0.0, 0.0);

    /***************************************/

    int i;
    double beta;
    if (myId == root) 
	{
		//vec是精确解，进程0设置为1
        InitDoubles (vec, dim, 1.0, 0.0);
        InitDoubles (sol1, dim, 0.0, 0.0);
        InitDoubles (sol2, dim, 0.0, 0.0);
        //  ProdSparseMatrixVectorByRows (mat, 0, vec, sol1);
    }
	//这个for是算右端项吗？
    int k=0;
    int *vptrM = matL.vptr;
    for (int i=0; i < matL.dim1; i++) 
	{
        for(int j=vptrM[i]; j<vptrM[i+1]; j++) 
		{
            sol1L[k] += matL.vval[j];
        }
        k++;
    }
/*
    if(myId==root)
    {
	for(int i=0;i<dimL;i++)
	{
	    printf("%lf\n",sol1L[i]);
	}	
    }
*/
	//0进程把sol2（初值为0）分发给每个进程的sol2L
    MPI_Scatterv (sol2, vdimL, vdspL, MPI_DOUBLE, sol2L, dimL, MPI_DOUBLE, root, MPI_COMM_WORLD);
/*
    //0
    if(myId==root)
    {
	for(int i=0;i<dimL;i++)
	{
	    printf("%lf\n",sol2L[i]);
	}	
    }
*/
	//matL是分配给每个进程的A，sol2L是每个进程的初始值，sol1L是每个进程的右端项
    ConjugateGradient (matL, sol2L, sol1L, vdimL, vdspL, rbuf, myId);
/*
    //inf and -nan why?
    if(myId==root)
    {
	for(int i=0;i<dimL;i++)
	{
	    printf("%lf\n",sol2L[i]);
	}	
    }
*/
	
	//sol2L应该是每个处理器得到的解
    // Error computation
    for (i=0; i<dimL; i++) sol2L[i] -= 1.0;
/*
    if(myId==root)
    {
	for(int i=0;i<dimL;i++)
	{
	    printf("%lf\n",sol2L[i]);
	}	
    }
*/

    // ReproAllReduce -- Begin
    // 因为前面是inf and -nan 所以这里需要2*
    std::vector<int64_t> h_superacc(2*exblas::BIN_COUNT);
    //算的出来可以用下面
    //std::vector<int64_t> h_superacc(exblas::BIN_COUNT);
    exblas::cpu::exdot (dimL, sol2L, sol2L, &h_superacc[0]);
    int imin=exblas::IMIN, imax=exblas::IMAX;
    exblas::cpu::Normalize(&h_superacc[0], imin, imax);
    if (myId == 0) 
	{
        MPI_Reduce (MPI_IN_PLACE, &h_superacc[0], exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    } 
	else 
	{
        MPI_Reduce (&h_superacc[0], NULL, exblas::BIN_COUNT, MPI_LONG, MPI_SUM, 0, MPI_COMM_WORLD);
    }
    if (myId == 0) 
	{
        beta = exblas::cpu::Round( &h_superacc[0] );
    }
    MPI_Bcast(&beta, 1, MPI_DOUBLE, 0, MPI_COMM_WORLD);
    // ReproAllReduce -- End

    beta = sqrt(beta);
	//%a以16进制格式输出近似解与真解的误差的平方和开方
    if (myId == 0) 
        //printf ("Error: %a\n", beta);
        printf ("Error: %.10e\n", beta);

    /***************************************/
    // Freeing memory
    RemoveDoubles (&sol2L); 
	RemoveDoubles (&sol1L); 
	RemoveDoubles (&vecL);
    RemoveInts (&vdspL); RemoveInts (&vdimL); 
    if (myId == root) 
	{
        RemoveDoubles (&sol2); 
		RemoveDoubles (&sol1); 
		RemoveDoubles (&vec);
        RemoveSparseMatrix (&mat); 
		RemoveSparseMatrix (&sym);
    } 
	else 
	{
        RemoveDoubles (&sol2); 
		RemoveDoubles (&vec);
    }

    MPI_Finalize ();

    return 0;
}

