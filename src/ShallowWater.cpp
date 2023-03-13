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
    _u(new double[_n]), _v(new double[_n]), _h(new double[_n])
    // _ncd(_nx-2), _mcd(_nx), _ldcd(7), _cd(new double[_ncd*_mcd])
{
    // check inputs
    if (
        _dt < 0.0 || _t < 0.0 || _t < _dt ||
        _nx < 2 || _ny < 2 || // _nx != _ny ||
        _ic < 1 || 4 < _ic
    ) throw std::invalid_argument("Invalid argument.");

    // central difference

    setInitialConditions();
}

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
            // plane waves propagating in x
            for (int i = 0; i < _nx; i++) {
                double x = i * _dx;
                x -= 50.0;
                x *= x;
                double h = 10.0 + std::exp(-x / 25.0);
                for (int j = 0; j < _ny; j++) {
                    int id = _colMajToArrId(i, j);
                    _h[id] = h;
                }
            }
            break;
        case 2:
            // plane waves propagating in y
            for (int j = 0; j < _ny; j++) {
                double y = j * _dy;
                y -= 50.0;
                y *= y;
                double h = 10.0 + std::exp(-y / 25.0);
                for (int i = 0; i < _nx; i++) {
                    int id = _colMajToArrId(i, j);
                    _h[id] = h;
                }
            }
            break;
        case 3:
            // single droplet
            for (int j = 0; j < _ny; j++) {
                for (int i = 0; i < _nx; i++) {
                    double x = _dx * i; double y = _dy * j;
                    x -= 50.0; y -= 50.0;
                    x *= x; y *= y;
                    double h = 10.0 + std::exp(-(x + y) / 25.0);
                    int id = _colMajToArrId(i, j);
                    _h[id] = h;
                }
            }
            break;
        default:
            // double droplet
            for (int j = 0; j < _ny; j++) {
                for (int i = 0; i < _nx; i++) {
                    double x = _dx * i; double y = _dy * j;
                    double h = 10.0 + 
                        std::exp(-((x-25.0)*(x-25.0) + (y-25.0)*(y-25.0)) / 25.0) +
                        std::exp(-((x-75.0)*(x-75.0) + (y-75.0)*(y-75.0)) / 25.0);
                    int id = _colMajToArrId(i, j);
                    _h[id] = h;
                }
            }
            break;
    }
}

void ShallowWater::timeIntegrate() {
    int nxr = _nx-2; int nyr = _ny-2;
    double* dudx = new double[nxr * nyr];
    double* dudy = new double[nxr * nyr];
    double* dvdx = new double[nxr * nyr];
    double* dvdy = new double[nxr * nyr];
    double* dhdx = new double[nxr * nyr];
    double* dhdy = new double[nxr * nyr];

    double t = 0.0;
    while (t < _t) {

        t += _dt;
    }

    delete[] dudx; delete[] dudy; delete[] dvdx; delete[] dvdy; delete[] dhdx; delete[] dhdy;
}

inline int ShallowWater::_colMajToArrId(int i, int j) {
    return i + j*_nx;
}

// accessors

int ShallowWater::get_nx() {
    return _nx;
}

int ShallowWater::get_ny() {
    return _ny;
}

double ShallowWater::get_dx() {
    return _dx;
}

double ShallowWater::get_dy() {
    return _dy;
}

double ShallowWater::get_u(int i, int j) {
    return _u[i + j*_nx];
}

double ShallowWater::get_v(int i, int j) {
    return _v[i + j*_nx];
}

double ShallowWater::get_h(int i, int j) {
    return _h[i + j*_nx];
}
