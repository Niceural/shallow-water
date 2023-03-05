#ifndef INITIAL_CONDITION_H
#define INITIAL_CONDITION_H

#include <cmath>

class InitialCondition {
    public:
        double g(double x, double y);
};

class PlaneWavesPropagatingInXInitialCondition: public InitialCondition {
    public:
        double g(double x, double y);
};

#endif // INITIAL_CONDITION_H