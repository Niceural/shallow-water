#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include "./matrices/GeneralMatrix.h"
#include "./matrices/SquareBandedMatrix.h"
#include <omp.h>

/// @brief Performs finite difference with a BLAS or a loop implementation.
class FiniteDifference {
    private:
        // variables for the loop implementation
        const double _invdx;
        const double _invdy;
        // loop
        const double _clx[6];
        const double _cly[6];

        // variables for the blas implementation
        // blas wrt x
        SquareBandedMatrix _dx_d;
        GeneralMatrix _dx_t1;
        GeneralMatrix _dx_t2;
        void _generateDx(const double dx);

        // blas wrt y
        SquareBandedMatrix _dy_d;
        GeneralMatrix _dy_t1;
        GeneralMatrix _dy_t2;
        void _generateDy(const double dy);
    
    public:
        FiniteDifference(const int m, const int n, const double dx, const double dy);
        ~FiniteDifference();

        void loop(
            const GeneralMatrix& U,
            const GeneralMatrix& V,
            const GeneralMatrix& H,
            GeneralMatrix& dUdx,
            GeneralMatrix& dUdy,
            GeneralMatrix& dVdx,
            GeneralMatrix& dVdy,
            GeneralMatrix& dHdx,
            GeneralMatrix& dHdy
        );

        void blas(
            const GeneralMatrix& U,
            const GeneralMatrix& V,
            const GeneralMatrix& H,
            GeneralMatrix& dUdx,
            GeneralMatrix& dUdy,
            GeneralMatrix& dVdx,
            GeneralMatrix& dVdy,
            GeneralMatrix& dHdx,
            GeneralMatrix& dHdy
        );
};

#endif // FINITE_DIFFERENCE_H