#ifndef CENTRAL_DIFFERENCE_2D_H
#define CENTRAL_DIFFERENCE_2D_H

#include "blasRoutines.h"
#include "GeneralMatrix.h"
#include "SquareBandedMatrix.h"
#include "TriangularPackedMatrix.h"
#include <cmath>

class CentralDifference2D {
    private:
        // parameters
        const double _dt; /// Time-step.
        const double _t; /// Total integration time.
        const int _nx; /// Number of grid points in x.
        const int _ny; /// Number of grid points in y.
        const int _n; /// Total number of grid points.
        const int _ic; /// Index of the initial condition to use (1-4).
        const double _dx; /// Constant point spacing along x.
        const double _dy; /// Constant point spacing along y.
        const double _g; /// Acceleration due to gravity.

        // grid
        GeneralMatrix _U; /// Matrix (nx * ny) of x-component of velocity.
        GeneralMatrix _V; /// Matrix (nx * ny) of y-component of velocity.
        GeneralMatrix _H; /// Matrix (nx * ny) of surface height.

        // central difference wrt x
        SquareBandedMatrix _cd_x_d;
        GeneralMatrix _cd_x_t1;
        GeneralMatrix _cd_x_t2;
        // central difference wrt y
        SquareBandedMatrix _cd_y_d;
        GeneralMatrix _cd_y_t1;
        GeneralMatrix _cd_y_t2;

        // central difference matrices
        GeneralMatrix _dUdx;
        GeneralMatrix _dUdy;
        GeneralMatrix _dVdx;
        GeneralMatrix _dVdy;
        GeneralMatrix _dHdx;
        GeneralMatrix _dHdy;

        int _gbTo1d(int i, int j);
    
    public:
        CentralDifference2D(const double dt, const double t, const int nx, const int ny, const int nc);
        ~CentralDifference2D();

        void setInitialConditions();
        void timeIntegrate();

        void integrateWrtX(GeneralMatrix& A, GeneralMatrix& dAdx);
        void integrateWrtY(GeneralMatrix& A, GeneralMatrix& dAdy);
};

#endif // CENTRAL_DIFFERENCE_2D_H
