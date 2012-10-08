These files are Fortran and Julia implementations of a noise whitening
procedure based on a particular form of covariance.  The procedure consists
of a Cholesky factorization of a covariance matrix and matrix-vector
multiplication to a signal vector of the lower triangular factor or its
inverse.  The special form of the covariance enables these routines to
take O(n) operations, rather than O(n^3) for the factorization and O(n^2)
for the multiplications.

The aim is to compare Julia performance to Fortran performance for this
simple code.  The first Julia file covar1.jl contains a naive implementation,
while the second covar2.jl makes most loops explicit (to improve performance).
