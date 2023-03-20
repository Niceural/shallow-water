#ifndef SHALLOW_WATER_H
#define SHALLOW_WATER_H
#define CONST_G 9.81

#include "blasRoutines.h"
#include "matrices/GeneralMatrix.h"
#include "matrices/SquareBandedMatrix.h"
#include "FiniteDifference.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "MultiQuantityMatrix.h"

class ShallowWater {
    private:
        const double _dx;
        const double _dy;

        MultiQuantityMatrix _grid;
        FiniteDifference _fd;

    public:
        ShallowWater(const int nx, const int ny, const double dx, const double dy);

        void setInitialConditions(const int ic, const double meanH);
        void timeIntegrate(const bool loopBlas, const double dt, const double t);

        void exportGrid(const std::string& fname);
        void test();
};

#endif // SHALLOW_WATER_H
