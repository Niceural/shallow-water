#ifndef SHALLOW_WATER_H
#define SHALLOW_WATER_H

#include <cmath>
#include "blasRoutines.h"
#include <stdexcept>

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
        int _ncd;
        int _mcd;
        int _ldcd;
        double* _cd; /// Banded matrix to perform the 6th-order central difference scheme with respect to x (size: nx-2 x ny).

    public:
        ShallowWater(double dt, double t, int nx, int ny, int nc);
        ~ShallowWater();

        void setInitialConditions();
};

#endif // SHALLOW_WATER_H