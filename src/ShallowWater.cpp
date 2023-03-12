#include "../include/ShallowWater.h"

ShallowWater::ShallowWater(
    double dt, double t,
    int nx, int ny, 
    int ic
):
    _dt(dt), _t(t),
    _nx(nx), _ny(ny),
    _ic(ic),
    _n(nx*ny), _dx(1.0), _dy(1.0), _g(9.81),
    _u(new double[_n]), _v(new double[_n]), _h(new double[_n]),
    _ncd(_nx-2), _mcd(_nx), _ldcd(7), _cd(new double[_ncd*_mcd])
{
    // check inputs
    if (
        _dt < 0.0 || _t < 0.0 || _t < _dt ||
        _nx < 2 || _ny < 2 || _nx != _ny ||
        _ic < 1 || 4 < _ic
    ) throw std::invalid_argument("Invalid argument.");

    // central difference matrix
    double val[] = {
        1.0/60.0/_dx, -3.0/20.0/_dx, 3.0/4.0/_dx, 
        0.0,
        -3.0/4.0/_dx, 3.0/20.0/_dx, -1.0/60.0/_dx
    };
    for (int i = 0; i < _ncd; i++) {
        for (int j = 0; j < _ldcd; j++) {
            _cd[i*_ldcd + j] = val[j];
        }
    }
}

ShallowWater::~ShallowWater() {
    delete[] _u;
    delete[] _v;
    delete[] _h;
    delete[] _cd;
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