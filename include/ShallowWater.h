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

/// @brief Solves the shallow water equations.
class ShallowWater {
    private:
        // parameters
        const double _dx; /// Constant point spacing along x.
        const double _dy; /// Constant point spacing along y.

        // grid
        GeneralMatrix _U; /// Matrix of x-component of velocity of size (N_x by N_y).
        GeneralMatrix _V; /// Matrix of y-component of velocity of size (N_x by N_y).
        GeneralMatrix _H; /// Matrix of surface height of size (N_x by N_y).

        // finite difference
        FiniteDifference _fd; /// Object responsible for performing finite difference.
    
    public:
        /// @brief Constructor.
        /// @param nx Number of points in the x direction N_x.
        /// @param ny Number of points in the y direction N_y.
        /// @param dx Constant point spacing in x.
        /// @param dy Constant point spacing in y.
        ShallowWater(const int nx, const int ny, const double dx, const double dy);

        /// @brief Sets U, V, and H to the initial conditions.
        /// @param ic Initial condition id (1-4).
        void setInitialConditions(const int ic);

        /// @brief Solves the equations.
        /// @param loopBlas 0 -> use the loop implementation of finite difference, 1 -> use the BLAS implementation of finite difference.
        /// @param dt Time step.
        /// @param t Total integration time.
        void timeIntegrate(const bool loopBlas, const double dt, const double t);

        /// @brief Exports the grid values U, V and H in the format prescribed by the brief.
        /// @param fname Name of the file the values are written to. If the file already exists it will be overwritten.
        void exportData(const std::string& fname);
};

#endif // SHALLOW_WATER_H
