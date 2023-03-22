#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include "./matrices/GeneralMatrix.h"
#include "./matrices/SquareBandedMatrix.h"
#include <omp.h>

class FiniteDifference {
    private:
        // const int _m;
        // const int _n;
        // const double _dx;
        // const double _dy;
        const double _clx[6];
        const double _cly[6];

        // SquareBandedMatrix _dx_d;
        // GeneralMatrix _dx_t1;
        // GeneralMatrix _dx_t2;
        // void _generateDx();

        // SquareBandedMatrix _dy_d;
        // GeneralMatrix _dy_t1;
        // GeneralMatrix _dy_t2;
        // void _generateDy();
    
        // void _performWrtXLoop(const GeneralMatrix& A, GeneralMatrix& dAdx);
        // void _performWrtYLoop(const GeneralMatrix& A, GeneralMatrix& dAdy);

        // void _performWrtXBlas(const GeneralMatrix& A, GeneralMatrix& dAdx);
        // void _performWrtYBlas(const GeneralMatrix& A, GeneralMatrix& dAdy);

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