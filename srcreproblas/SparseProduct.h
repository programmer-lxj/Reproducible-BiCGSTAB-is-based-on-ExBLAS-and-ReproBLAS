#ifndef SparseProductTip

#define SparseProductTip 1

typedef struct
	{
		int dim1, dim2;
		int *vptr;
		int *vpos;
		double *vval;
	} SparseMatrix, *ptr_SparseMatrix;

/*********************************************************************************/

// This routine creates a sparseMatrix from the next parameters
// * numR defines the number of rows
// * numC defines the number of columns
// * numE defines the number of nonzero elements
// * msr indicates if the MSR is the format used to the sparse matrix
// If msr is actived, numE doesn't include the diagonal elements
// The parameter index indicates if 0-indexing or 1-indexing is used.
extern void CreateSparseMatrix (ptr_SparseMatrix p_spr, int index, int numR, int numC, int numE, 
																	int msr);

// This routine liberates the memory related to matrix spr
extern void RemoveSparseMatrix (ptr_SparseMatrix spr);

/*********************************************************************************/

// This routine creates de sparse matrix dst from the symmetric matrix spr.
// The parameters indexS and indexD indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices.
extern void DesymmetrizeSparseMatrices (SparseMatrix src, int indexS, ptr_SparseMatrix dst, 
																					int indexD);

/*********************************************************************************/

extern int ReadMatrixHB (char *filename, ptr_SparseMatrix p_spr);

/*********************************************************************************/

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVector2 (SparseMatrix spr, int index, double *vec, double *res);

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVectorByRows (SparseMatrix spr, int index, double *vec, double *res);

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVectorByRows_OMP (SparseMatrix spr, int index, double *vec, double *res);

/*********************************************************************************/

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVectorByCols (SparseMatrix spr, int index, double *vec, double *res);

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
extern void ProdSparseMatrixVectorByCols_OMP (SparseMatrix spr, int index, double *vec, double *res);

/*********************************************************************************/

extern void GetDiagonalSparseMatrix2 (SparseMatrix spr, int shft, double *diag, int *posd);

/*********************************************************************************/

#endif
