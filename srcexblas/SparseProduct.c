#include <stdio.h>
#include <stdlib.h>
#include <math.h>

/// #include "InputOutput.h"
#include "ScalarVectors.h"
#include "hb_io.h"
#include "SparseProduct.h"

/*********************************************************************************/

// This routine creates a sparseMatrix from the next parameters
// * numR defines the number of rows
// * numC defines the number of columns
// * numE defines the number of nonzero elements
// * msr indicates if the MSR is the format used to the sparse matrix
// If msr is actived, numE doesn't include the diagonal elements
// The parameter index indicates if 0-indexing or 1-indexing is used.
void CreateSparseMatrix (ptr_SparseMatrix p_spr, int index, int numR, int numC, int numE, int msr) {
//	printf (" index = %d , numR = %d , numC = %d , numE = %d\n", index, numR, numC, numE);
	// The scalar components of the structure are initiated
	p_spr->dim1 = numR; p_spr->dim2 = numC; 
	// Only one malloc is made for the vectors of indices
	CreateInts (&(p_spr->vptr), numE+numR+1);
	// The first component of the vectors depends on the used format
	*(p_spr->vptr) = ((msr)? (numR+1): 0) + index;
	p_spr->vpos = p_spr->vptr + ((msr)? 0: (numR+1));
	// The number of nonzero elements depends on the format used
	CreateDoubles (&(p_spr->vval), numE+(numR+1)*msr);
}

// This routine liberates the memory related to matrix spr
void RemoveSparseMatrix (ptr_SparseMatrix spr) {
	// First the scalar are initiated
	spr->dim1 = -1; spr->dim2 = -1; 
	// The vectors are liberated
	RemoveInts (&(spr->vptr)); RemoveDoubles (&(spr->vval)); 
}

/*********************************************************************************/

// This routine creates de sparse matrix dst from the symmetric matrix spr.
// The parameters indexS and indexD indicate, respectivaly, if 0-indexing or 1-indexing is used
// to store the sparse matrices.
void DesymmetrizeSparseMatrices (SparseMatrix src, int indexS, ptr_SparseMatrix dst, int indexD) {
	int n = src.dim1, nnz = 0;
	int *sizes = NULL;
	int *pp1 = NULL, *pp2 = NULL, *pp3 = NULL, *pp4 = NULL, *pp5 = NULL;
	int i, j, dim, indexDS = indexD - indexS;
	double *pd3 = NULL, *pd4 = NULL;

	// The vector sizes is created and initiated
	CreateInts (&sizes, n); InitInts (sizes, n, 0, 0);
	// This loop counts the number of elements in each row
	pp1 = src.vptr; pp3 = src.vpos + *pp1 - indexS;
	pp2 = pp1 + 1 ; pp4 = sizes - indexS;
	for (i=indexS; i<(n+indexS); i++) {
		// The size of the corresponding row is accumulated
		dim = (*pp2 - *pp1); pp4[i] += dim;
		// Now each component of the row is analyzed
		for (j=0; j<dim; j++) {
			// The nondiagonals elements define another element in the graph
			if (*pp3 != i) pp4[*pp3]++;
			pp3++;
		}
		pp1 = pp2++; 
	}
	
	// Compute the number of nonzeros of the new sparse matrix
	nnz = AddInts (sizes, n);
	// Create the new sparse matrix
	CreateSparseMatrix (dst, indexD, n, n, nnz, 0);
	// Fill the vector of pointers
	CopyInts (sizes, (dst->vptr) + 1, n);
	dst->vptr[0] = indexD; TransformLengthtoHeader (dst->vptr, n);
	// The vector sizes is initiated with the beginning of each row
	CopyInts (dst->vptr, sizes, n);
	// This loop fills the contents of vector vpos
	pp1 = src.vptr; pp3 = src.vpos + *pp1 - indexS; 
	pp2 = pp1 + 1 ; pp4 = dst->vpos - indexD; pp5 = sizes - indexS;
	pd3 = src.vval  + *pp1 - indexS; pd4 = dst->vval - indexD;
	for (i=indexS; i<(n+indexS); i++) {
		dim = (*pp2 - *pp1);
		for (j=0; j<dim; j++) {
			// The elements in the i-th row
			pp4[pp5[i]  ] = *pp3+indexDS; 
			pd4[pp5[i]++] = *pd3; 
			if (*pp3 != i) {
				// The nondiagonals elements define another element in the graph
				pp4[pp5[*pp3]  ] = i+indexDS;
				pd4[pp5[*pp3]++] = *pd3;
			}
			pp3++; pd3++;
		}
		pp1 = pp2++;
	}
	// The memory related to the vector sizes is liberated
	RemoveInts (&sizes);
}

/*********************************************************************************/

int ReadMatrixHB (char *filename, ptr_SparseMatrix p_spr) {
  int *colptr = NULL;
  double *exact = NULL;
  double *guess = NULL;
  int indcrd;
  char *indfmt = NULL;
  FILE *input;
  char *key = NULL;
  char *mxtype = NULL;
  int ncol;
  int neltvl;
  int nnzero;
  int nrhs;
  int nrhsix;
  int nrow;
  int ptrcrd;
  char *ptrfmt = NULL;
  int rhscrd;
  char *rhsfmt = NULL;
  int *rhsind = NULL;
  int *rhsptr = NULL;
  char *rhstyp = NULL;
  double *rhsval = NULL;
  double *rhsvec = NULL;
  int *rowind = NULL;
  char *title = NULL;
  int totcrd;
  int valcrd;
  char *valfmt = NULL;
  double *values = NULL;

	printf ("\nTEST09\n");
	printf ("  HB_FILE_READ reads all the data in an HB file.\n");
	printf ("  HB_FILE_MODULE is the module that stores the data.\n");

	input = fopen (filename, "r");
	if ( !input ) {
		printf ("\n TEST09 - Warning!\n Error opening the file %s .\n", filename);
		return -1;
	}

	hb_file_read ( input, &title, &key, &totcrd, &ptrcrd, &indcrd,
									&valcrd, &rhscrd, &mxtype, &nrow, &ncol, &nnzero, &neltvl,
									&ptrfmt, &indfmt, &valfmt, &rhsfmt, &rhstyp, &nrhs, &nrhsix,
									&colptr, &rowind, &values, &rhsval, &rhsptr, &rhsind, &rhsvec,
									&guess, &exact );
	fclose (input);

	// Conversion Fortran to C
	CopyShiftInts (colptr, colptr, nrow+1, -1);
	CopyShiftInts (rowind, rowind, nnzero, -1);

	//  Data assignment
	p_spr->dim1 = nrow  ; p_spr->dim2 = ncol  ; 
	p_spr->vptr = colptr; p_spr->vpos = rowind; p_spr->vval = values; 

	//  Memory liberation
	free (exact ); free (guess ); free (indfmt);
	free (key   ); free (mxtype); free (ptrfmt);
	free (rhsfmt); free (rhsind); free (rhsptr);
	free (rhstyp); free (rhsval); free (rhsvec);
	free (title ); free (valfmt);
	
	return 0;
}

/*********************************************************************************/

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVector2 (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j;
	int *pp1 = spr.vptr, *pp2 = pp1+1, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pvec = vec - index, *pd2 = res;
	double *pd1 = spr.vval + *pp1 - index;

	// If the MSR format is used, first the diagonal has to be processed
	if (spr.vptr == spr.vpos)
		VvecDoubles (1.0, spr.vval, vec, 1.0, res, spr.dim1);

	for (i=0; i<spr.dim1; i++) {
		// The dot product between the row i and the vector vec is computed
		aux = 0.0;
		for (j=*pp1; j<*pp2; j++)
			aux += *(pd1++) * pvec[*(pi1++)];
//		for (j=spr.vptr[i]; j<spr.vptr[i+1]; j++)
//			aux += spr.vval[j] * pvec[spr.vpos[j]];
		// Accumulate the obtained value on the result
		*(pd2++) += aux; pp1 = pp2++;
	}
}

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVectorByRows (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j, dim = spr.dim1;
	int *pp1 = spr.vptr, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pvec = vec + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index;

	// Process all the rows of the matrix
	for (i=0; i<dim; i++) {
		// The dot product between the row i and the vector vec is computed
		aux = 0.0;
		for (j=pp1[i]; j<pp1[i+1]; j++)
			aux = fma(pd1[j], pvec[pi1[j]], aux);
		// Accumulate the obtained value on the result
		res[i] += aux; 
	}
}

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVectorByRows_OMP (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j, dim = spr.dim1;
	int *pp1 = spr.vptr, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pvec = vec + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index;

	// Process all the rows of the matrix
	#pragma omp parallel for private(j, aux)
	for (i=0; i<dim; i++) {
		// The dot product between the row i and the vector vec is computed
		aux = 0.0;
		for (j=pp1[i]; j<pp1[i+1]; j++)
			aux += pd1[j] * pvec[pi1[j]];
		// Accumulate the obtained value on the result
		res[i] += aux; 
	}
}

/*void ProdSparseMatrixVectorByRows_OMPTasks (SparseMatrix spr, int index, double *vec, double *res, int bm) {
	int i, dim = spr.dim1;

	// Process all the rows of the matrix
	//#pragma omp taskloop grainsize(bm) 
	for ( i=0; i<dim; i+=bm) {
		int cs = dim - i;
		int c = cs < bm ? cs : bm;
//	for (i=0; i<dim; i++) {
	  #pragma omp task depend(inout:res[i:i+c-1]) //shared(c)
		{
//	printf("Task SPMV ---- i: %d, c: %d \n", i, c);
		  int *pp1 = spr.vptr, *pi1 = spr.vpos + *pp1 - index;
	    double aux, *pvec = vec + *pp1 - index;
	    double *pd1 = spr.vval + *pp1 - index;
	  	// The dot product between the row i and the vector vec is computed
		  aux = 0.0;
			for(int idx=i; idx < i+c; idx++){
			//	printf("Task SPMV ---- idx: %d\n", idx);
		  	for (int j=pp1[idx]; j<pp1[idx+1]; j++)
			  	aux += pd1[j] * pvec[pi1[j]];
		  	// Accumulate the obtained value on the result
		  	res[idx] += aux; 
	  	}
		}
	}
}
*/

void ProdSparseMatrixVectorByRows_OMPTasks (SparseMatrix spr, int index, double *vec, double *res, int bm) {
	int i, j, dim = spr.dim1;
	int *pp1 = spr.vptr, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pvec = vec + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index;

	// Process all the rows of the matrix
	#pragma omp taskloop grainsize(bm) 
	for ( i=0; i<dim; i++ ) {
	  	// The dot product between the row i and the vector vec is computed
		  aux = 0.0;
		  for (j=pp1[i]; j<pp1[i+1]; j++)
			  aux += pd1[j] * pvec[pi1[j]];
		  // Accumulate the obtained value on the result
		  res[i] += aux; 
	}
}


/*********************************************************************************/

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVectorByCols (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j, dim = spr.dim1;
	int *pp1 = spr.vptr, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pres = res + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index;

	// Process all the columns of the matrix
	for (i=0; i<dim; i++) {
		// The result is scaled by the column i and the scalar vec[i]
		aux = vec[i];
		for (j=pp1[i]; j<pp1[i+1]; j++)
			pres[pi1[j]] += pd1[j] * aux;
	}
}

// This routine computes the product { res += spr * vec }.
// The parameter index indicates if 0-indexing or 1-indexing is used,
void ProdSparseMatrixVectorByCols_OMP (SparseMatrix spr, int index, double *vec, double *res) {
	int i, j, dim = spr.dim1;
	int *pp1 = spr.vptr, *pi1 = spr.vpos + *pp1 - index;
	double aux, *pres = res + *pp1 - index;
	double *pd1 = spr.vval + *pp1 - index;

	// Process all the columns of the matrix
	#pragma omp parallel for private(j, aux)
	for (i=0; i<dim; i++) {
		// The result is scaled by the column i and the scalar vec[i]
		for (j=pp1[i]; j<pp1[i+1]; j++) {
			aux = vec[i] * pd1[j];
			#pragma omp atomic
				pres[pi1[j]] += aux;
		}
	}
}

/*********************************************************************************/

void GetDiagonalSparseMatrix2 (SparseMatrix spr, int shft, double *diag, int *posd) {
	int i, j, dim = (spr.dim1 < spr.dim2)? spr.dim1 : spr.dim2;
	int *pp1 = NULL, *pp2 = NULL, *pi1 = NULL, *pi2 = posd; 
	double *pd1 = NULL, *pd2 = diag;

	if (spr.vptr == spr.vpos)
		CopyDoubles (spr.vval, diag, spr.dim1);
	else {
		pp1 = spr.vptr; pp2 = pp1+1; j = (*pp2-*pp1);
		pi1 = spr.vpos+*pp1; pd1 = spr.vval+*pp1; 
		for (i=0; i<dim; i++) {
//			while ((j > 0) && (*pi1 < i))
			while ((j > 0) && (*pi1 < (i+shft)))
				{ pi1++; pd1++; j--; }
//			*(pd2++) = ((j > 0) && (*pi1 == i))? *pd1: 0.0;
			*(pd2++) = ((j > 0) && (*pi1 == (i+shft)))? *pd1: 0.0;
//			*(pi2++) = ((j > 0) && (*pi1 == i))? *pp2-j: -1;
			*(pi2++) = ((j > 0) && (*pi1 == (i+shft)))? *pp2-j: -1;
			pi1 += j; pd1 += j; pp1 = (pp2++); j = (*pp2-*pp1);
		}
	}
}

/*********************************************************************************/
