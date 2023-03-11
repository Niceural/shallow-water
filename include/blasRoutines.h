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
}

#endif // BLAS_ROUTINES_H
