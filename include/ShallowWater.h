#ifndef SHALLOW_WATER_H
#define SHALLOW_WATER_H

#include "blasRoutines.h"
#include "matrices/GeneralMatrix.h"
#include "matrices/SquareBandedMatrix.h"
#include "CentralDifference.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

class ShallowWater {
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

        CentralDifference _cd;
    
        void _timeIntegrateLoop();
        void _timeIntegrateBlas();

    public:
        ShallowWater(const double dt, const double t, const int nx, const int ny, const int nc);
        ~ShallowWater();

        void setInitialConditions();
        void timeIntegrate();

        void exportData(const std::string& fname);

        // getters
        void test();
        GeneralMatrix getH() const;
};

#endif // SHALLOW_WATER_H
