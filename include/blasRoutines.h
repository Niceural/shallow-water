#ifndef BLAS_ROUTINES_H
#define BLAS_ROUTINES_H

#define F77NAME(x) x##_

extern "C" {
    // x = alpha * x
    void F77NAME(dscal) (
        const int& n,
        const double& alpha,
        double* x,
        const int& incx
    );

    void F77NAME(dgbmv)(
        const char& trans,
        const int& m,
        const int& n,
        const int& kl,
        const int& ku,
        const double& alpha,
        double* A,
        const int& lda,
        double* x,
        const int& incx,
        const int& beta,
        double* y,
        const int& incy
    );
}

#endif // BLAS_ROUTINES_H
