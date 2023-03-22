#ifndef FINITE_DIFFERENCE_H
#define FINITE_DIFFERENCE_H

#include "MultiQuantityMatrix.h"
#include "./matrices/GeneralMatrix.h"
#include "./matrices/SquareBandedMatrix.h"

class FiniteDifference {
    private:
        const int _m;
        const int _n;
        const double _dx;
        const double _dy;

        void _centralDifferenceLoop(MultiQuantityMatrix& grid);

        SquareBandedMatrix _cd_d;
        GeneralMatrix _cd_t1;
        GeneralMatrix _cd_t2;
        void _centralDifferenceBlas(MultiQuantityMatrix& grid);

    public:
        FiniteDifference(const int m, const int n, const double dx, const double dy);
        ~FiniteDifference();

        void centralDifference(const bool loopBlas, MultiQuantityMatrix& grid);
};

#endif // FINITE_DIFFERENCE_H