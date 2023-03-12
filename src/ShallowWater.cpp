#include "../include/ShallowWater.h"

ShallowWater::ShallowWater(
    double dt, int t,
    int nx, int ny, 
    int ic
):
    _dt(dt), _t(t),
    _nx(nx), _ny(ny),
    _ic(ic),
    _n(nx*ny), _dx(1.0), _dy(1.0), _g(9.81),
    _u(new double[_n]),
    _v(new double[_n]),
    _h(new double[_n])
{}

ShallowWater::~ShallowWater() {
    delete[] _u;
    delete[] _v;
    delete[] _h;
}

void ShallowWater::setInitialConditions() {
    // set u to zero
    F77NAME(dscal)(_n, 0.0, _u, 1);

    // set v to zero
    F77NAME(dscal)(_n, 0.0, _v, 1);

    // set h to the initial surface height for each test cases
    switch (_ic) {
        case 1:
            for (int i = 0; i < _nx; i++) {
                double x = i * _dx;
                double h = std::exp(- (x-50.0) * (x-50.0) / 25);
                for (int j = 0; j < _ny; j++)
                    _h[i + j*_nx] = h;
            }
            break;
        case 2:
            break;
        case 3:
            break;
        default:
            break;
    }
}