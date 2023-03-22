#ifndef SHALLOW_WATER_H
#define SHALLOW_WATER_H

#include "blasRoutines.h"
#include "matrices/GeneralMatrix.h"
#include "matrices/SquareBandedMatrix.h"
#include "FiniteDifference.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>

class ShallowWater {
    private:
        // parameters
        const double _dx; /// Constant point spacing along x.
        const double _dy; /// Constant point spacing along y.

        // grid
        GeneralMatrix _U; /// Matrix (nx * ny) of x-component of velocity.
        GeneralMatrix _V; /// Matrix (nx * ny) of y-component of velocity.
        GeneralMatrix _H; /// Matrix (nx * ny) of surface height.

        FiniteDifference _fd;
    
    public:
        ShallowWater(const int nx, const int ny, const double dx, const double dy);
        ~ShallowWater();

        void setInitialConditions(const int ic);
        void timeIntegrate(const bool loopBlas, const double dt, const double t);

        void exportData(const std::string& fname);
        void test();
};

#endif // SHALLOW_WATER_H
