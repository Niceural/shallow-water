#include "../include/InitialCondition.h"

double PlaneWavesPropagatingInXInitialCondition::g(double x, double y) {
    double sq = x - 50.0;
    return std::exp(- sq * sq / 25.0);
}