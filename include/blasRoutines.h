#ifndef BLAS_ROUTINES_H
#define BLAS_ROUTINES_H

#define F77NAME(x) x##_

extern "C" {
    // copy x into y
    void F77NAME(dcopy) (
        const int& n,
        const double* x,
        const int& incx,
        double* y,
        const int& incy
    );

    // x = alpha * x
    void F77NAME(dscal) (
        const int& n,
        const double& alpha,
        double* x,
        const int& incx
    );

    // y = a*x + y
    void F77NAME(daxpy) (
        const int& n,
        const double& a,
        const double* x,
        const int& incx,
        double* y,
        const int& incy
    );

    // general matrix vector multiply
    void F77NAME(dgemv) (
        const char& trans,
        const int& m,
        const int& n,
        const double& alpha,
        const double* A,
        const int& lda,
        const double* x,
        const int& incx,
        const double& beta,
        double* y,
        const int& incy
    );

    void F77NAME(dgbmv)(
        const char& trans,
        const int& m,
        const int& n,
        const int& kl,
        const int& ku,
        const double& alpha,
        const double* A,
        const int& lda,
        const double* x,
        const int& incx,
        const int& beta,
        double* y,
        const int& incy
    );

    void F77NAME(dtpmv) (
        const char& uplo,
        const char& trans,
        const char& diag,
        const int& n,
        const double* A,
        double* x,
        const int& incx
    );

    // symmetric banded matrix vector product
    void F77NAME(dsbmv) (
        const char& UPLO, // 'U' for upper diagonal, 'L' for lower diagonal
        const int& N, // order of matrix A
        const int& K, // number of super-diagonals of the matrix A
        const double& ALPHA, // scalar multiplying A
        const double* A, // array of dimensions (LDA, A)
        const int& LDA,
        const double* X, // vector X
        const int& INCX, // increment of X
        const double& BETA, // scalar multiplying Y (careful to Y when non zero)
        double* Y, // vector Y
        const int& INCY // increment of Y
    );

}

#endif // BLAS_ROUTINES_H
