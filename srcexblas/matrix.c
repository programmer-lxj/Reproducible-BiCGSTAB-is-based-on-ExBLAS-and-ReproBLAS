#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <assert.h>

//#include "global.h"
//#include "debug.h"

#include "matrix.h"
//#include "cg_aux_conLabel.h"

// finite-difference method for a 3D Poisson's equation with a 7, 19 or 27 point stencil
void generate_Poisson3D(ptr_SparseMatrix A, const int p, const int stencil_points, int dspL, int dimL, int dim)
{
	int p2 = p * p, i, j=0; //, pos=0;
	int pos = 0;
	int *vptr = A->vptr;
	const int    *stenc_c;
	const double *stenc_v;

	const int    stenc_c7[]  = { -p2,  -p,  -1,   0,   1,   p,  p2};
	const double stenc_v7[]  = { -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0};

	const double r = 1.0;
	const int    stenc_c19[] =
	{
		       -p2-p,          -p2-1,  -p2+0, -p2+1,          -p2+p,
		 -p-1,    -p,    -p+1,    -1,      0,     1,     p-1,     p,     p+1,
		        p2-p,           p2-1,   p2+0,  p2+1,           p2+p
	};
	const double stenc_v19[] =
	{
		     -(1+r),      -(1+r),    -(8*r-4),   -(1+r),      -(1+r),
		-2, -(6-2*r), -2, -(6-2*r), -(-32-16*r), -(6-2*r), -2, -(6-2*r), -2,
		     -(1+r),      -(1+r),    -(8*r-4),   -(1+r),      -(1+r)
	};

	const int    stenc_c27[] =
	{
		-p2-p-1, -p2-p, -p2-p+1, -p2-1,  -p2+0, -p2+1, -p2+p-1, -p2+p, -p2+p+1,
		   -p-1,    -p,    -p+1,    -1,      0,     1,     p-1,     p,     p+1,
		 p2-p-1,  p2-p,  p2-p+1,  p2-1,   p2+0,  p2+1,  p2+p-1,  p2+p,  p2+p+1
	};
	const double stenc_v27[] =
	{
		  -(2+r),  -(8-10*r),    -(2+r),  -(8-10*r),   -(100*r-40),  -(8-10*r),    -(2+r),  -(8-10*r),    -(2+r),
		-(20-2*r), -(80-20*r), -(20-2*r), -(80-20*r), -(-400-200*r), -(80-20*r), -(20-2*r), -(80-20*r), -(20-2*r),
		   -(2+r),  -(8-10*r),    -(2+r),  -(8-10*r),   -(100*r-40),  -(8-10*r),    -(2+r),  -(8-10*r),    -(2+r)
	};

	if( stencil_points == 7 )
	{
		stenc_c = stenc_c7;
		stenc_v = stenc_v7;
	}
	else if( stencil_points == 19 )
	{
		stenc_c = stenc_c19;
		stenc_v = stenc_v19;
	}
	else if( stencil_points == 27 )
	{
		stenc_c = stenc_c27;
		stenc_v = stenc_v27;
	}
	else
		// this should be impossible, but silences compiler warnings
		return;

	// to compute the nnz, we just need to know that each stencil point at distance |d| from the diagonal
	// will be excluded from the matrix on d lines, otherwise each stencil point is on each line

	// let's only do the part here.
	//for(j=start_row; j<end_row; j++)
	printf("Generate matrix ---- dim: %d, dimL: %d, dspL: %d\n", dim, dimL, dspL);
	for(j=0; j<dimL; j++) //M
	{
		vptr[j] = pos;
	//	printf("A->vptr[%d]: %d \n", j, A->vptr[j]);
		for(i=0; i<stencil_points; i++){
			int val = j + dspL + stenc_c[i];
			//long val = j + dspL + stenc_c[i];
		//	printf("j: %d dspL: %d, stenc_c[%d]: %d, dim: %d \n", j, dspL, i, stenc_c[i], dim);
			if( val >= 0 && val < dim )
			{
				A->vpos[pos] = val;
				A->vval[pos] = stenc_v[i];
		//		printf("j: %d, stenc_c[%d]: %d, A->vpos[%d]: %d, A->vval[%d]: %f\n", j, i, stenc_c[i], pos, A->vpos[pos], pos, A->vval[pos]);
				pos++;
			}
		}
	}

	// point to just beyond last element
	//MA->vptr[j] = pos;
	vptr[j] = pos;
	//A->elemc = pos;
//	A->vptr[j-start_row] = pos;
}

// finite-difference method for a 3D Poisson's equation with a 7, 19 or 27 point stencil
void generate_Poisson3D_filled(ptr_SparseMatrix A, const int p, const int stencil_points, int band_width, int dspL, int dimL, int dim)
{
	int p2 = p * p, i, j=0; //, pos=0;
	int pos = 0;
	int *vptr = A->vptr;
	const double value = 0.1;
	const int    *stenc_c;
	const double *stenc_v;

	const int    stenc_c7[]  = { -p2,  -p,  -1,   0,   1,   p,  p2};
	const double stenc_v7[]  = { -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0};

	const double r = 1.0;
	const int    stenc_c19[] =
	{
		       -p2-p,          -p2-1,  -p2+0, -p2+1,          -p2+p,
		 -p-1,    -p,    -p+1,    -1,      0,     1,     p-1,     p,     p+1,
		        p2-p,           p2-1,   p2+0,  p2+1,           p2+p
	};
	const double stenc_v19[] =
	{
		     -(1+r),      -(1+r),    -(8*r-4),   -(1+r),      -(1+r),
		-2, -(6-2*r), -2, -(6-2*r), -(-32-16*r), -(6-2*r), -2, -(6-2*r), -2,
		     -(1+r),      -(1+r),    -(8*r-4),   -(1+r),      -(1+r)
	};

	const int    stenc_c27[] =
	{
		-p2-p-1, -p2-p, -p2-p+1, -p2-1,  -p2+0, -p2+1, -p2+p-1, -p2+p, -p2+p+1,
		   -p-1,    -p,    -p+1,    -1,      0,     1,     p-1,     p,     p+1,
		 p2-p-1,  p2-p,  p2-p+1,  p2-1,   p2+0,  p2+1,  p2+p-1,  p2+p,  p2+p+1
	};
	const double stenc_v27[] =
	{
		  -(2+r),  -(8-10*r),    -(2+r),  -(8-10*r),   -(100*r-40),  -(8-10*r),    -(2+r),  -(8-10*r),    -(2+r),
		-(20-2*r), -(80-20*r), -(20-2*r), -(80-20*r), -(-400-200*r), -(80-20*r), -(20-2*r), -(80-20*r), -(20-2*r),
		   -(2+r),  -(8-10*r),    -(2+r),  -(8-10*r),   -(100*r-40),  -(8-10*r),    -(2+r),  -(8-10*r),    -(2+r)
	};

	if( stencil_points == 7 )
	{
		stenc_c = stenc_c7;
		stenc_v = stenc_v7;
	}
	else if( stencil_points == 19 )
	{
		stenc_c = stenc_c19;
		stenc_v = stenc_v19;
	}
	else if( stencil_points == 27 )
	{
		stenc_c = stenc_c27;
		stenc_v = stenc_v27;
	}
	else
		// this should be impossible, but silences compiler warnings
		return;

	// to compute the nnz, we just need to know that each stencil point at distance |d| from the diagonal
	// will be excluded from the matrix on d lines, otherwise each stencil point is on each line

	// let's only do the part here.
	//for(j=start_row; j<end_row; j++)
	printf("Generate matrix ---- dim: %d, dimL: %d, dspL: %d, band_width: %d \n", dim, dimL, dspL, band_width);
	for(j=0; j<dimL; j++) //M
	{
		int iii;
		//long iii;
		int jjj = j + dspL;
		int prv = 0;
		vptr[j] = pos;
		//printf("j: %d, dspL: %d, pos: %d, dim: %d, band_width: %d\n", j, dspL, pos, dim, band_width);
		//printf(stderr, "A->vptr[%d]: %ld \n", j, vptr[j]);
		for(i=0; i<stencil_points; i++){
			int	val = j + dspL + stenc_c[i];
			//long	val = j + dspL + stenc_c[i];
			//printf("j: %d, i: %d, val: %d, dspL: %d, pos: %d, stenc_c[%d]: %d, dim: %d, band_width: %d\n", j, i, val, dspL, pos, i, stenc_c[i], dim, band_width);
			//printf("j: %d, val: %d, dspL: %d, pos: %ld , stenc_c[%d]: %d, dim: %d, band_width: %d \n", j, val, dspL, pos, i, stenc_c[i], dim, band_width);
//			if( j + dspL + stenc_c[i] >= 0 && j + dspL + stenc_c[i] < dim )
			if( val >= 0 && val < dim )
			{
				// Analyzing if val is into the band
//				if( val >= (j-band_width) && val <= (j+band_width) ) {
				if( val >= (jjj-band_width) && val <= (jjj+band_width) ) {
					// Adding the elements in the band which are previous to val
//					int kk1 = ((j-band_width) <    0)?     0: (j-band_width);
					int kk1 = ((jjj-band_width) <    0)?     0: (jjj-band_width);
					if (prv != 0) kk1 = prv;
					//printf("Entra en if kk1: %d, val: %d j: %d, band_width: %d, j-band_width: %d, j+band_width: %d \n", kk1, val, j, band_width, j-band_width, j+band_width);
					for (iii=kk1; iii<val; iii++) {
						A->vpos[pos] = iii;
						A->vval[pos] = value;
						pos++;
					}
					prv = val + 1;
					// Choosing the correct value to val 
//					if ( val == j ) {
					if ( val == jjj ) {
						//printf("Diagonal val: %d, pos: %ld \n", val, pos);
						A->vpos[pos] = val;
						A->vval[pos] = stenc_v[i] + band_width * value;
						//A->vval[pos] = stenc_v[i] + (2 * band_width * value) / 20;
						pos++;
					}	else {
						//printf("No Diagonal val: %d, j: %ld \n", val, pos);
						A->vpos[pos] = val;
						A->vval[pos] = stenc_v[i] + value;
						pos++;
					}
//				} else if (val < (j-band_width)) {
				} else if (val < (jjj-band_width)) {
					// Choosing the correct value to val 
					A->vpos[pos] = val;
					A->vval[pos] = stenc_v[i];
					pos++;
					// prv = val + 1;
				} else {
					// Adding the elements in the band which are previous to val
//					int kk2 = ((j+band_width) >= dim)? dim-1: (j+band_width);
					int kk2 = ((jjj+band_width) >= dim)? dim-1: (jjj+band_width);
//					if (prv == 0) prv = j+1;
					if (prv == 0) prv = jjj+1;
					//printf("No entra en if prv: %d, kk2: %d, val: %d, j: %d, band_width: %d, j-band_width: %d, j+band_width: %d \n", prv, kk2, val, j, band_width, j-band_width, j+band_width);
					for (iii=prv; iii<=kk2; iii++) {
						A->vpos[pos] = iii;
						A->vval[pos] = value;
						pos++;
					}
					prv = kk2 + 1;
					// Choosing the correct value to val 
					A->vpos[pos] = val;
					A->vval[pos] = stenc_v[i];
					pos++;
				}
/*
//				A->vpos[pos] = j + stenc_c[i] + dspL;
				A->vpos[pos] = val;
				A->vval[pos] = stenc_v[i];
		//		printf("j: %d, stenc_c[%d]: %d, A->vpos[%d]: %d, A->vval[%d]: %f\n", j, i, stenc_c[i], pos, A->vpos[pos], pos, A->vval[pos]);
				pos++;
*/
			}
		}
//		if (prv <= (j+band_width)) {
		if (prv <= (jjj+band_width)) {
//			int kk2 = ((j+band_width) >= dim)? dim-1: (j+band_width);
			int kk2 = ((jjj+band_width) >= dim)? dim-1: (jjj+band_width);
//			if (prv == 0) prv = j+1;
			if (prv == 0) prv = jjj+1;
			//printf("No entra en if prv: %d, kk2: %d, j: %d, band_width: %d, j-band_width: %d, j+band_width: %d \n", prv, kk2, j, band_width, j-band_width, j+band_width);
			for (iii=prv; iii<=kk2; iii++) {
				A->vpos[pos] = iii;
				A->vval[pos] = value;
				pos++;
			}
		}
	}

	// point to just beyond last element
	//A->vptr[j] = pos;
	vptr[j] = pos;
	//A->elemc = pos;
//	A->vptr[j-start_row] = pos;
	printf("FIN Generate matrix ---- dim: %d, dimL: %d, dspL: %d, band_width: %d \n", dim, dimL, dspL, band_width);
}

void generate_Poisson3D_perm(ptr_SparseMatrix A, const int p, const int stencil_points, int init, int step, int dimL, int dim)
{
	int p2 = p * p, i, j=0; //, pos=0;
	int pos = 0;
	int *vptr = A->vptr;
	const int    *stenc_c;
	const double *stenc_v;

	const int    stenc_c7[]  = { -p2,  -p,  -1,   0,   1,   p,  p2};
	//const double stenc_v7[]  = { 1.0, 1.0, 1.0,-6.0, 1.0, 1.0, 1.0};
	const double stenc_v7[]  = { -1.0, -1.0, -1.0, 6.0, -1.0, -1.0, -1.0};

	const double r = 1.0;
	const int    stenc_c19[] =
	{
		       -p2-p,          -p2-1,  -p2+0, -p2+1,          -p2+p,
		 -p-1,    -p,    -p+1,    -1,      0,     1,     p-1,     p,     p+1,
		        p2-p,           p2-1,   p2+0,  p2+1,           p2+p
	};
	//const double stenc_v19[] =
//	{
//		     1+r,      1+r,    8*r-4,   1+r,      1+r,
//		2, 6-2*r, 2, 6-2*r, -32-16*r, 6-2*r, 2, 6-2*r, 2,
//		     1+r,      1+r,    8*r-4,   1+r,      1+r
//	};
	const double stenc_v19[] =
	{
		     -(1+r),      -(1+r),    -(8*r-4),   -(1+r),      -(1+r),
		-2, -(6-2*r), -2, -(6-2*r), -(-32-16*r), -(6-2*r), -2, -(6-2*r), -2,
		     -(1+r),      -(1+r),    -(8*r-4),   -(1+r),      -(1+r)
	};

	const int    stenc_c27[] =
	{
		-p2-p-1, -p2-p, -p2-p+1, -p2-1,  -p2+0, -p2+1, -p2+p-1, -p2+p, -p2+p+1,
		   -p-1,    -p,    -p+1,    -1,      0,     1,     p-1,     p,     p+1,
		 p2-p-1,  p2-p,  p2-p+1,  p2-1,   p2+0,  p2+1,  p2+p-1,  p2+p,  p2+p+1
	};
	//const double stenc_v27[] =
//	{
//		   2+r,  8-10*r,    2+r,  8-10*r,   100*r-40,  8-10*r,    2+r,  8-10*r,    2+r,
//		20-2*r, 80-20*r, 20-2*r, 80-20*r, -400-200*r, 80-20*r, 20-2*r, 80-20*r, 20-2*r,
//		   2+r,  8-10*r,    2+r,  8-10*r,   100*r-40,  8-10*r,    2+r,  8-10*r,    2+r
//	};
	const double stenc_v27[] =
	{
		  -(2+r),  -(8-10*r),    -(2+r),  -(8-10*r),   -(100*r-40),  -(8-10*r),    -(2+r),  -(8-10*r),    -(2+r),
		-(20-2*r), -(80-20*r), -(20-2*r), -(80-20*r), -(-400-200*r), -(80-20*r), -(20-2*r), -(80-20*r), -(20-2*r),
		   -(2+r),  -(8-10*r),    -(2+r),  -(8-10*r),   -(100*r-40),  -(8-10*r),    -(2+r),  -(8-10*r),    -(2+r)
	};

	if( stencil_points == 7 )
	{
		stenc_c = stenc_c7;
		stenc_v = stenc_v7;
	}
	else if( stencil_points == 19 )
	{
		stenc_c = stenc_c19;
		stenc_v = stenc_v19;
	}
	else if( stencil_points == 27 )
	{
		stenc_c = stenc_c27;
		stenc_v = stenc_v27;
	}
	else
		// this should be impossible, but silences compiler warnings
		return;

	// to compute the nnz, we just need to know that each stencil point at distance |d| from the diagonal
	// will be excluded from the matrix on d lines, otherwise each stencil point is on each line
	/*A->nnz = stencil_points * p3 ; 
	for(i=0; i<stencil_points; i++)
		A->nnz -= abs(stenc_c[i]);*/ //M
	//A->nnz = A->nnz / nProcs; //M

	// let's only do the part here.
	//for(j=start_row; j<end_row; j++)
	int dimB = (step == 1) ? 0: (dim / step);
	int resB = (step == 1) ? 0: (dim % step);
	int row=init;
//	printf ("init = %d , step = %d , dimB = %d , dim = %d\n", init, step, dimB, dim);
	for(j=0; j<dimL; j++) //M
//	for(j=init; j<dim; j+=step) //M
	{
		//A->vptr[j-start_row] = pos;
		//MA->vptr[j] = pos;
		vptr[j] = pos;
	//	printf("A->vptr[%d]: %d \n", j, A->vptr[j]);
//		printf ("(%d) row = %d , pos = %d\n", init, row-1, pos);
		for(i=0; i<stencil_points; i++){
			int val = row + stenc_c[i];
			int val1 = (val / step);
			int val2 = (val % step);
//			int k = (val2 * dimB + val1) ;
			int k = (val2 * dimB + val1 + ((val2 < resB)? val2: resB)) ;
			//long k = (val2 * dimB + val1 + ((val2 < resB)? val2: resB)) ;
		//	printf("j: %d dspL: %d, stenc_c[%d]: %d, dim: %d \n", j, dspL, i, stenc_c[i], dim);
//			if( k >= 0 && k < dim )
			if( val >= 0 && val < dim )
			{
				A->vpos[pos] = k;
				A->vval[pos] = stenc_v[i];
		//		printf("j: %d, stenc_c[%d]: %d, A->vpos[%d]: %d, A->vval[%d]: %f\n", j, i, stenc_c[i], pos, A->vpos[pos], pos, A->vval[pos]);
				pos++;
			}
		}
		row += step;
	}

	// point to just beyond last element
	//M A->vptr[j] = pos;
	vptr[j] = pos;
	//A->elemc = pos;
//	A->vptr[j-start_row] = pos;
}

void allocate_matrix(const int m, const int n, const int nnz, ptr_SparseMatrix A)
{
	A->dim1 = m;
	A->dim2 = n;
	//A->elemc = nnz;
	//long *vptr = A->vptr;

	//A->vptr = (int*)calloc((n+1), sizeof(int));
	A->vptr = (int*)calloc((n+1), sizeof(int));
	//vptr = (long*)calloc((n+1), sizeof(long));

	A->vpos = (int*)calloc(nnz, sizeof(int));
	//A->vpos = (long *)calloc(nnz, sizeof(long));
	A->vval = (double*)calloc(nnz, sizeof(double));

/*	if( ! A->vval ) 
	{
		fprintf(stderr, "Allocating vval failed !\n");
		exit(2);
	}

	if (! A->vpos ) 
	{
		fprintf(stderr, "Allocating vpos failed !\n");
		exit(2);
	}

	if (! A->vptr )
	{
		fprintf(stderr, "Allocating vptr failed !\n");
		exit(2);
	}*/
	if( ! A->vval || ! A->vpos || ! A->vptr )
	{
		fprintf(stderr, "Allocating sparse matrix of size %d rows and %d non-zeros failed !\n", n, nnz);
		exit(2);
	}
	fprintf(stderr, "Matrix allocated\n");
}

void ScaleFirstRowCol(SparseMatrix A, int despL, int dimL, int myId, int root, double factor){   
// To generate ill-conditioned matrices
  int i;

  if (myId == root) {
    for (i=A.vptr[0]; i< A.vptr[1]; i++)
       A.vval[i] *= factor;
  }
  if (despL == 0) {
    i = 0;
    while((i < dimL) && (A.vpos[A.vptr[i]] == 0)) {
      A.vval[A.vptr[i]] *= factor;
      i++;
    }
  }
}

/*
void deallocate_matrix(Matrix *A)
{
	free(A->r);
	free(A->c);

	if( A->v )
		free(A->v);
}*/


