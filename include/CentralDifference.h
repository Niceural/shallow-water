#ifndef CENTRAL_DIFFERENCE_H
#define CENTRAL_DIFFERENCE_H

#include "./matrices/GeneralMatrix.h"
#include "./matrices/SquareBandedMatrix.h"
#include <omp.h>

class CentralDifference {
    private:
        const int _m;
        const int _n;
        const double _dx;
        const double _dy;

        SquareBandedMatrix _dx_d;
        GeneralMatrix _dx_t1;
        GeneralMatrix _dx_t2;
        void _generateDx();

        SquareBandedMatrix _dy_d;
        GeneralMatrix _dy_t1;
        GeneralMatrix _dy_t2;
        void _generateDy();
    
        void _performWrtXLoop(const GeneralMatrix& A, GeneralMatrix& dAdx);
        void _performWrtYLoop(const GeneralMatrix& A, GeneralMatrix& dAdy);

        void _performWrtXBlas(const GeneralMatrix& A, GeneralMatrix& dAdx);
        void _performWrtYBlas(const GeneralMatrix& A, GeneralMatrix& dAdy);

    public:
        CentralDifference(const int m, const int n, const double dx, const double dy);
        ~CentralDifference();

        void performWrtX(const bool loopBlas, const GeneralMatrix& A, GeneralMatrix& dAdx);
        void performWrtY(const bool loopBlas, const GeneralMatrix& A, GeneralMatrix& dAdy);
};

#endif // CENTRAL_DIFFERENCE_H