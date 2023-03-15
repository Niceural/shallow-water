#include "../include/ShallowWater.h"

ShallowWater::ShallowWater(
    const double dt, const double t,
    const int nx, const int ny,
    const int ic
):
    // parameters
    _dt(dt), _t(t),
    _nx(nx), _ny(ny), _n(nx * ny),
    _ic(ic),
    _dx(1.0), _dy(1.0),
    _g(9.81),

    // grid
    _U(GeneralMatrix(nx, ny)), _V(GeneralMatrix(nx, ny)), _H(GeneralMatrix(nx, ny)),
    _cd(CentralDifference(_nx, _ny, _dx, _dy))
{
    if (
        _dt < 0.0 || _t < 0.0 || _t < _dt ||
        _nx < 2 || _ny < 2 ||
        _ic < 1 || 4 < _ic
    ) throw std::invalid_argument("Invalid argument.");
}

ShallowWater::~ShallowWater() {}

void ShallowWater::setInitialConditions() {
    // set u to zero
    F77NAME(dscal)(_n, 0.0, _U.getPointer(0), 1);

    // set v to zero
    F77NAME(dscal)(_n, 0.0, _V.getPointer(0), 1);

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
                    _H.set(i, j, h);
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
                    _H.set(i, j, h);
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
                    _H.set(i, j, h);
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
                    _H.set(i, j, h);
                }
            }
            break;
    }
}

void ShallowWater::timeIntegrate() {
    setInitialConditions();

    // central difference
    GeneralMatrix dUdx(_nx, _ny);
    GeneralMatrix dUdy(_nx, _ny);
    GeneralMatrix dVdx(_nx, _ny);
    GeneralMatrix dVdy(_nx, _ny);
    GeneralMatrix dHdx(_nx, _ny);
    GeneralMatrix dHdy(_nx, _ny);

    // runge kutta U
    GeneralMatrix k1U(_nx, _ny);
    GeneralMatrix k2U(_nx, _ny);
    GeneralMatrix k3U(_nx, _ny);
    GeneralMatrix k4U(_nx, _ny);

    // runge kutta V
    GeneralMatrix k1V(_nx, _ny);
    GeneralMatrix k2V(_nx, _ny);
    GeneralMatrix k3V(_nx, _ny);
    GeneralMatrix k4V(_nx, _ny);

    // runge kutta H
    GeneralMatrix k1H(_nx, _ny);
    GeneralMatrix k2H(_nx, _ny);
    GeneralMatrix k3H(_nx, _ny);
    GeneralMatrix k4H(_nx, _ny);

    // temporary variables
    GeneralMatrix tempU(_nx, _ny);
    GeneralMatrix tempV(_nx, _ny);
    GeneralMatrix tempH(_nx, _ny);

    double t = 0.0;
    while (t < _t) {
        // perform central difference
        _cd.performWrtX(_U, dUdx);
        _cd.performWrtX(_V, dVdx);
        _cd.performWrtX(_H, dHdx);
        _cd.performWrtY(_U, dUdy);
        _cd.performWrtY(_V, dVdy);
        _cd.performWrtY(_H, dHdy);

        //--------- k1
        // copy U, V, H, to tempU, tempV, tempH
        F77NAME(dcopy)(_U.size(), _U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(dcopy)(_V.size(), _V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(dcopy)(_H.size(), _H.getPointer(0), 1, tempH.getPointer(0), 1);
        // k1U
        F77NAME(dcopy)(dHdx.size(), dHdx.getPointer(0), 1, k1U.getPointer(0), 1); // k1U = dHdx
        F77NAME(dscal)(k1U.size(), -_g, k1U.getPointer(0), 1); // k1U *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dUdx, 1.0, k1U); // k1U -= U * dUdx
        elementWiseMultiplicationLoop(-1.0, tempV, dUdy, 1.0, k1U); // k1U -= V * dUdy
        // k1V
        F77NAME(dcopy)(dHdy.size(), dHdy.getPointer(0), 1, k1V.getPointer(0), 1); // k1V = dHdy
        F77NAME(dscal)(k1V.size(), -_g, k1V.getPointer(0), 1); // k1V *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dVdx, 1.0, k1V); // k1V -= U * dVdx
        elementWiseMultiplicationLoop(-1.0, tempV, dVdy, 1.0, k1V); // k1V -= V * dVdy
        // k1H
        elementWiseMultiplicationLoop(-1.0, tempU, dHdx, 0.0, k1H);
        elementWiseMultiplicationLoop(-1.0, tempH, dUdx, 1.0, k1H);
        elementWiseMultiplicationLoop(-1.0, tempV, dHdy, 1.0, k1H);
        elementWiseMultiplicationLoop(-1.0, tempH, dVdy, 1.0, k1H);

        //--------- k2
        // copy U, V, H, to tempU, tempV, tempH
        F77NAME(dcopy)(_U.size(), _U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(dcopy)(_V.size(), _V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(dcopy)(_H.size(), _H.getPointer(0), 1, tempH.getPointer(0), 1);
        // a_n += dt * k1 / 2
        F77NAME(daxpy)(tempU.size(), _dt*0.5, k1U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(daxpy)(tempV.size(), _dt*0.5, k1V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(daxpy)(tempH.size(), _dt*0.5, k1H.getPointer(0), 1, tempH.getPointer(0), 1);
        // k2U
        F77NAME(dcopy)(dHdx.size(), dHdx.getPointer(0), 1, k2U.getPointer(0), 1); // k2U = dHdx
        F77NAME(dscal)(k2U.size(), -_g, k2U.getPointer(0), 1); // k2U *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dUdx, 1.0, k2U); // k2U -= U * dUdx
        elementWiseMultiplicationLoop(-1.0, tempV, dUdy, 1.0, k2U); // k2U -= V * dUdy
        // k2V
        F77NAME(dcopy)(dHdy.size(), dHdy.getPointer(0), 1, k2V.getPointer(0), 1); // k2V = dHdy
        F77NAME(dscal)(k2V.size(), -_g, k2V.getPointer(0), 1); // k2V *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dVdx, 1.0, k2V); // k2V -= U * dVdx
        elementWiseMultiplicationLoop(-1.0, tempV, dVdy, 1.0, k2V); // k2V -= V * dVdy
        // k2H
        elementWiseMultiplicationLoop(-1.0, tempU, dHdx, 0.0, k2H);
        elementWiseMultiplicationLoop(-1.0, tempH, dUdx, 1.0, k2H);
        elementWiseMultiplicationLoop(-1.0, tempV, dHdy, 1.0, k2H);
        elementWiseMultiplicationLoop(-1.0, tempH, dVdy, 1.0, k2H);

        //--------- k3
        // copy U, V, H, to tempU, tempV, tempH
        F77NAME(dcopy)(_U.size(), _U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(dcopy)(_V.size(), _V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(dcopy)(_H.size(), _H.getPointer(0), 1, tempH.getPointer(0), 1);
        // a_n += dt * k2 / 2
        F77NAME(daxpy)(tempU.size(), _dt*0.5, k2U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(daxpy)(tempV.size(), _dt*0.5, k2V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(daxpy)(tempH.size(), _dt*0.5, k2H.getPointer(0), 1, tempH.getPointer(0), 1);
        // k3U
        F77NAME(dcopy)(dHdx.size(), dHdx.getPointer(0), 1, k3U.getPointer(0), 1); // k3U = dHdx
        F77NAME(dscal)(k3U.size(), -_g, k3U.getPointer(0), 1); // k3U *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dUdx, 1.0, k3U); // k3U -= U * dUdx
        elementWiseMultiplicationLoop(-1.0, tempV, dUdy, 1.0, k3U); // k3U -= V * dUdy
        // k3V
        F77NAME(dcopy)(dHdy.size(), dHdy.getPointer(0), 1, k3V.getPointer(0), 1); // k3V = dHdy
        F77NAME(dscal)(k3V.size(), -_g, k3V.getPointer(0), 1); // k3V *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dVdx, 1.0, k3V); // k3V -= U * dVdx
        elementWiseMultiplicationLoop(-1.0, tempV, dVdy, 1.0, k3V); // k3V -= V * dVdy
        // k3H
        elementWiseMultiplicationLoop(-1.0, tempU, dHdx, 0.0, k3H);
        elementWiseMultiplicationLoop(-1.0, tempH, dUdx, 1.0, k3H);
        elementWiseMultiplicationLoop(-1.0, tempV, dHdy, 1.0, k3H);
        elementWiseMultiplicationLoop(-1.0, tempH, dVdy, 1.0, k3H);

        //--------- k4
        // copy U, V, H, to tempU, tempV, tempH
        F77NAME(dcopy)(_U.size(), _U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(dcopy)(_V.size(), _V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(dcopy)(_H.size(), _H.getPointer(0), 1, tempH.getPointer(0), 1);
        // a_n += dt * k3 / 2
        F77NAME(daxpy)(tempU.size(), _dt, k3U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(daxpy)(tempV.size(), _dt, k3V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(daxpy)(tempH.size(), _dt, k3H.getPointer(0), 1, tempH.getPointer(0), 1);
        // k4U
        F77NAME(dcopy)(dHdx.size(), dHdx.getPointer(0), 1, k4U.getPointer(0), 1); // k4U = dHdx
        F77NAME(dscal)(k4U.size(), -_g, k4U.getPointer(0), 1); // k4U *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dUdx, 1.0, k4U); // k4U -= U * dUdx
        elementWiseMultiplicationLoop(-1.0, tempV, dUdy, 1.0, k4U); // k4U -= V * dUdy
        // k4V
        F77NAME(dcopy)(dHdy.size(), dHdy.getPointer(0), 1, k4V.getPointer(0), 1); // k4V = dHdy
        F77NAME(dscal)(k4V.size(), -_g, k4V.getPointer(0), 1); // k4V *= -g
        elementWiseMultiplicationLoop(-1.0, tempU, dVdx, 1.0, k4V); // k4V -= U * dVdx
        elementWiseMultiplicationLoop(-1.0, tempV, dVdy, 1.0, k4V); // k4V -= V * dVdy
        // k4H
        elementWiseMultiplicationLoop(-1.0, tempU, dHdx, 0.0, k4H);
        elementWiseMultiplicationLoop(-1.0, tempH, dUdx, 1.0, k4H);
        elementWiseMultiplicationLoop(-1.0, tempV, dHdy, 1.0, k4H);
        elementWiseMultiplicationLoop(-1.0, tempH, dVdy, 1.0, k4H);

        // add k1
        F77NAME(dcopy)(k1U.size(), k1U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(dcopy)(k1V.size(), k1V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(dcopy)(k1H.size(), k1H.getPointer(0), 1, tempH.getPointer(0), 1);
        // add k2
        F77NAME(daxpy)(tempU.size(), 2.0, k2U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(daxpy)(tempV.size(), 2.0, k2V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(daxpy)(tempH.size(), 2.0, k2H.getPointer(0), 1, tempH.getPointer(0), 1);
        // add k3
        F77NAME(daxpy)(tempU.size(), 2.0, k3U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(daxpy)(tempV.size(), 2.0, k3V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(daxpy)(tempH.size(), 2.0, k3H.getPointer(0), 1, tempH.getPointer(0), 1);
        // add k4
        F77NAME(daxpy)(tempU.size(), 1.0, k4U.getPointer(0), 1, tempU.getPointer(0), 1);
        F77NAME(daxpy)(tempV.size(), 1.0, k4V.getPointer(0), 1, tempV.getPointer(0), 1);
        F77NAME(daxpy)(tempH.size(), 1.0, k4H.getPointer(0), 1, tempH.getPointer(0), 1);
        // add to U, V, and H
        F77NAME(daxpy)(_U.size(), 1.0/6.0*_dt, tempU.getPointer(0), 1, _U.getPointer(0), 1);
        F77NAME(daxpy)(_V.size(), 1.0/6.0*_dt, tempV.getPointer(0), 1, _V.getPointer(0), 1);
        F77NAME(daxpy)(_H.size(), 1.0/6.0*_dt, tempH.getPointer(0), 1, _H.getPointer(0), 1);

        t += _dt;
    }
}

void ShallowWater::exportData(const std::string& fname) {
    std::ofstream file;
    file.open(fname);

    for (int j = 0; j < _ny; j++) {
        double y = j * _dy;            
        for (int i = 0; i < _nx; i++) {
            double x = i * _dx;
            file << x << " " << y << " ";
            file << _U.get(i, j) << " ";
            file << _V.get(i, j) << " ";
            file << _H.get(i, j) << std::endl;
        }
        file << std::endl;
    }

    file.close();
}

GeneralMatrix ShallowWater::getH() const {
    return _H;
}
