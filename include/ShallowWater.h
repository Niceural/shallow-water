#ifndef SHALLOW_WATER_H
#define SHALLOW_WATER_H
#define CONST_G 9.81

#include "blasRoutines.h"
#include "matrices/GeneralMatrix.h"
#include "matrices/SquareBandedMatrix.h"
#include "CentralDifference.h"
#include <cmath>
#include <string>
#include <iostream>
#include <fstream>
#include "Domain.h"

class ShallowWater {
    private:
        Domain _domain;

    public:
        ShallowWater(const int nx, const int ny, const double dx, const double dy);

        void setInitialConditions(const int ic, const double meanH);
        void timeIntegrate();

        void exportData(const std::string& fname);

        // getters
        void test();
        GeneralMatrix getH() const;
};

#endif // SHALLOW_WATER_H
