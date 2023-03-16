#ifndef CENTRAL_DIFFERENCE_H
#define CENTRAL_DIFFERENCE_H

#include "./matrices/GeneralMatrix.h"
#include "./matrices/SquareBandedMatrix.h"

class CentralDifference {
    private:
        const int _m;
        const int _n;
        const double _dx;
        const double _dy;

        SquareBandedMatrix _dx_d;
        GeneralMatrix _dx_t1;
        GeneralMatrix _dx_t2;
        void _generateDx(const double dx);

        SquareBandedMatrix _dy_d;
        GeneralMatrix _dy_t1;
        GeneralMatrix _dy_t2;
        void _generateDy(const double dy);
    
    public:
        CentralDifference(const int m, const int n, const double dx, const double dy);
        ~CentralDifference();

        void performWrtXLoop(const GeneralMatrix& A, GeneralMatrix& dAdx);
        void performWrtYLoop(const GeneralMatrix& A, GeneralMatrix& dAdy);

        void performWrtXBlas(const GeneralMatrix& A, GeneralMatrix& dAdx);
        void performWrtYBlas(const GeneralMatrix& A, GeneralMatrix& dAdy);
};

#endif // CENTRAL_DIFFERENCE_H