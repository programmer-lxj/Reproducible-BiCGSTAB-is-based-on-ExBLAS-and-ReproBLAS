# Reproducible-BiCGSTAB-is-based-on-ExBLAS-and-ReproBLAS
# Our work is based on https://github.com/riakymch/ReproCG#building-reprocg (ReproCG)

Building ReproBiCGSTAB

1. clone the git-repository into <ReproBiCGSTAB>

2. inside the src directory, invoke make to create CG_MPI executable


Example

automatically generated matrix arising from the finite-difference method of a 3D Poisson’s equation with 27 stencil points. This matrix has the number of rows/columns of the matrix equal to 256^3.
  
mpirun -np P --bind-to core ./ReproBiCGSTAB/src/CG_MPI ../Matrices/$mat 0 2 159 27

matrix from the Suite Sparse Matrix Collection

mpirun -np P --bind-to core ./ReproBiCGSTAB/src/CG_MPI MAT.rb 1

paper
Comparison of Reproducible Parallel Preconditioned BiCGSTAB Algorithm Based on ExBLAS and ReproBLAS.
https://dl.acm.org/doi/abs/10.1145/3578178.3578234
