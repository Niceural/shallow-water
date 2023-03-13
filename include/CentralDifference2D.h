#ifndef CENTRAL_DIFFERENCE_2D_H
#define CENTRAL_DIFFERENCE_2D_H

#include "blasRoutines.h"
#include "GeneralMatrix.h"
#include "GeneralBandedMatrix.h"
#include "TriangularPackedMatrix.h"

class CentralDifference2D {
    private:
        int _m;
        int _n;

        GeneralBandedMatrix _cd_x_d;
        // TriangularPackedMatrix _cd_x_t1;
        GeneralMatrix _cd_x_t1;
        // TriangularPackedMatrix _cd_x_t2;
        GeneralMatrix _cd_x_t2;

        GeneralMatrix _dUdx;
        GeneralMatrix _dUdy;
        GeneralMatrix _dVdx;
        GeneralMatrix _dVdy;
        GeneralMatrix _dHdx;
        GeneralMatrix _dHdy;

        int _gbTo1d(int i, int j);
    
    public:
        CentralDifference2D(int m, int n);
        ~CentralDifference2D();

        // void integrate(const double* U, const double* V, const double* H);
        // void integrateU(const double* U);
        // void integrateV(const double* V);
        // void integrateH(const double* H);

        void integrateWrtX(GeneralMatrix& A, GeneralMatrix& dAdx);
};

#endif // CENTRAL_DIFFERENCE_2D_H