#ifndef MATRIX_H_INCLUDED
#define MATRIX_H_INCLUDED

#include <stdio.h>
#include <SparseProduct.h>

typedef struct Matrix
{
	int n, m, *c, *r;
	long nnz;
	double *v;
} Matrix;

typedef enum
{
	FROM_FILE = 0,
	POISSON3D
} matrix_type;

void generate_Poisson3D(ptr_SparseMatrix A, const int p, const int stencil_points, int dspL, int dimL, int dim);

// memory utility functions
void allocate_matrix(const int n, const int m, const int nnz, ptr_SparseMatrix A);

void generate_Poisson3D_filled(ptr_SparseMatrix A, const int p, const int stencil_points, int band_width, int dspL, int dimL, int dim);

void generate_Poisson3D_perm(ptr_SparseMatrix A, const int p, const int stencil_points, int init, int step, int dimL, int dim);

void ScaleFirstRowCol(SparseMatrix A, int despL, int dimL, int myId, int root, double factor);
#endif // MATRIX_H_INCLUDED

