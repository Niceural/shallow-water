#ifndef SHALLOW_WATER_H
#define SHALLOW_WATER_H

#include <cmath>
#include "blasRoutines.h"
#include <stdexcept>
#include "GeneralBandedMatrix.h"

class ShallowWater {
    private:
        // inputs
        double _dt; /// Time-step to use.
        double _t; /// Total integration time.
        int _nx; /// Number of grid points in x.
        int _ny; /// Number of grid points in y.
        int _ic; /// Index of the initial condition to use (1-4).

        int _n; /// Total number of grid points (nx * ny).
        double _dx; /// Constant point spacing in x.
        double _dy; /// Constant point spacing in y.
        double _g; /// Acceleration due to gravity.

        // 
        double* _u; /// Matrix of size (nx * ny) of x-component of velocity at each grid point.
        double* _v; /// Matrix of size (nx * ny) of y-component of velocity at each grid point.
        double* _h; /// Matrix of size (nx * ny) of surface height, relative to a mean height of zero, at each grid point.

        // central difference
        // GeneralBandedMatrix _cdx; /// Banded matrix to perform the 6th-order central difference scheme with respect to x (size: nx-2 x ny).
        // GeneralBandedMatrix _cdy; /// Banded matrix to perform the 6th-order central difference scheme with respect to x (size: nx-2 x ny).

        inline int _colMajToArrId(int i, int j);

    public:
        ShallowWater(double dt, double t, int nx, int ny, int nc);
        ~ShallowWater();

        void setInitialConditions();
        void timeIntegrate();

        // accessors (for testing)
        int get_nx();
        int get_ny();

        double get_dx();
        double get_dy();

        double get_u(int i, int j);
        double get_v(int i, int j);
        double get_h(int i, int j);
};

#endif // SHALLOW_WATER_H