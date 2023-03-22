#include "../include/ShallowWater.h"

ShallowWater::ShallowWater(
    const int nx, const int ny,
    const double dx, const double dy
):
    // parameters
    _dx(dx), _dy(dy),
    // grid
    _U(GeneralMatrix(nx, ny)), _V(GeneralMatrix(nx, ny)), _H(GeneralMatrix(nx, ny)),
    _fd(FiniteDifference(nx, ny, dx, dy))
{
    if (nx<2 || ny<2 ||
        dx<0.0 || dy<0.0
    ) throw std::invalid_argument("Invalid argument.");
}

ShallowWater::~ShallowWater() {}

void ShallowWater::setInitialConditions(const int ic) {
    for (int id = 0; id < _U.size(); id++) {
        _U[id] = 0.0;
        _V[id] = 0.0;
    }

    // set h to the initial surface height for each test cases
    switch (ic) {
        case 1:
            // plane waves propagating in x
            for (int i = 0; i<_H.m(); i++) {
                double x = i * _dx;
                x -= 50.0;
                x *= x;
                double h = 10.0 + std::exp(-x / 25.0);
                for (int j = 0; j < _H.n(); j++) {
                    _H.set(i, j, h);
                }
            }
            break;

        case 2:
            // plane waves propagating in y
            for (int j=0; j<_H.n(); j++) {
                double y = j * _dy;
                y -= 50.0;
                y *= y;
                double h = 10.0 + std::exp(-y / 25.0);
                for (int i=0; i<_H.m(); i++) {
                    _H.set(i, j, h);
                }
            }
            break;

        case 3:
            // single droplet
            for (int j = 0; j < _H.n(); j++) {
                for (int i = 0; i < _H.m(); i++) {
                    double x = _dx * i; double y = _dy * j;
                    x -= 50.0; y -= 50.0;
                    x *= x; y *= y;
                    double h = 10.0 + std::exp(-(x + y) / 25.0);
                    _H.set(i, j, h);
                }
            }
            break;

        case 4:
            // double droplet
            for (int j = 0; j < _H.n(); j++) {
                for (int i = 0; i < _H.m(); i++) {
                    double x = _dx * i; double y = _dy * j;
                    double h = 10.0 + 
                        std::exp(-((x-25.0)*(x-25.0) + (y-25.0)*(y-25.0)) / 25.0) +
                        std::exp(-((x-75.0)*(x-75.0) + (y-75.0)*(y-75.0)) / 25.0);
                    _H.set(i, j, h);
                }
            }
            break;

        default:
            throw std::invalid_argument("Invalid initial condition");
            break;
    }
}

void ShallowWater::timeIntegrate(const bool loopBlas, const double dt, const double t) {
    // finite difference
    GeneralMatrix dUdx(_U.m(), _U.n());
    GeneralMatrix dUdy(_U.m(), _U.n());
    GeneralMatrix dVdx(_U.m(), _U.n());
    GeneralMatrix dVdy(_U.m(), _U.n());
    GeneralMatrix dHdx(_U.m(), _U.n());
    GeneralMatrix dHdy(_U.m(), _U.n());

    // runge kutta U
    GeneralMatrix k1U(_U.m(), _U.n());
    GeneralMatrix k2U(_U.m(), _U.n());
    GeneralMatrix k3U(_U.m(), _U.n());
    GeneralMatrix k4U(_U.m(), _U.n());

    // runge kutta V
    GeneralMatrix k1V(_U.m(), _U.n());
    GeneralMatrix k2V(_U.m(), _U.n());
    GeneralMatrix k3V(_U.m(), _U.n());
    GeneralMatrix k4V(_U.m(), _U.n());

    // runge kutta H
    GeneralMatrix k1H(_U.m(), _U.n());
    GeneralMatrix k2H(_U.m(), _U.n());
    GeneralMatrix k3H(_U.m(), _U.n());
    GeneralMatrix k4H(_U.m(), _U.n());

    // temporary variables
    GeneralMatrix tU(_U.m(), _U.n());
    GeneralMatrix tV(_U.m(), _U.n());
    GeneralMatrix tH(_U.m(), _U.n());

    const int n = _U.size();
    const double _g = 9.81;
    double ct = 0;
    while (ct < t) {
        //--------- k1
        if (loopBlas)
            _fd.blas(_U, _V, _H, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        else
            _fd.loop(_U, _V, _H, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        // k1
        for (int i = 0; i < n; i++) {
            k1U[i] = - (_U[i]*dUdx[i] + _V[i]*dUdy[i] + _g*dHdx[i]);
            k1V[i] = - (_U[i]*dVdx[i] + _V[i]*dVdy[i] + _g*dHdy[i]);
            k1H[i] = - (_U[i]*dHdx[i] + _H[i]*dUdx[i] + _V[i]*dHdy[i] + _H[i]*dVdy[i]);
        }

        //--------- k2
        // update U, V, and H
        for (int i = 0; i < n; i++) {
            tU[i] = _U[i] + 0.5*dt*k1U[i];
            tV[i] = _V[i] + 0.5*dt*k1V[i];
            tH[i] = _H[i] + 0.5*dt*k1H[i];
        }
        // central difference
        if (loopBlas)
            _fd.blas(tU, tV, tH, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        else
            _fd.loop(tU, tV, tH, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        // k2
        for (int i = 0; i < n; i++) {
            k2U[i] = - (tU[i]*dUdx[i] + tV[i]*dUdy[i] + _g*dHdx[i]);
            k2V[i] = - (tU[i]*dVdx[i] + tV[i]*dVdy[i] + _g*dHdy[i]);
            k2H[i] = - (tU[i]*dHdx[i] + tH[i]*dUdx[i] + tV[i]*dHdy[i] + tH[i]*dVdy[i]);
        }

        //--------- k3
        // update U, V, and H
        for (int i = 0; i < n; i++) {
            tU[i] = _U[i] + 0.5*dt*k2U[i];
            tV[i] = _V[i] + 0.5*dt*k2V[i];
            tH[i] = _H[i] + 0.5*dt*k2H[i];
        }
        // central difference
        if (loopBlas)
            _fd.blas(tU, tV, tH, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        else
            _fd.loop(tU, tV, tH, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        // k3
        for (int i = 0; i < n; i++) {
            k3U[i] = - (tU[i]*dUdx[i] + tV[i]*dUdy[i] + _g*dHdx[i]);
            k3V[i] = - (tU[i]*dVdx[i] + tV[i]*dVdy[i] + _g*dHdy[i]);
            k3H[i] = - (tU[i]*dHdx[i] + tH[i]*dUdx[i] + tV[i]*dHdy[i] + tH[i]*dVdy[i]);
        }

        //--------- k4
        // update U, V, and H
        for (int i = 0; i < n; i++) {
            tU[i] = _U[i] + dt*k3U[i];
            tV[i] = _V[i] + dt*k3V[i];
            tH[i] = _H[i] + dt*k3H[i];
        }
        // central difference
        if (loopBlas)
            _fd.blas(tU, tV, tH, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        else
            _fd.loop(tU, tV, tH, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);
        // k4
        for (int i = 0; i < n; i++) {
            k4U[i] = - (tU[i]*dUdx[i] + tV[i]*dUdy[i] + _g*dHdx[i]);
            k4V[i] = - (tU[i]*dVdx[i] + tV[i]*dVdy[i] + _g*dHdy[i]);
            k4H[i] = - (tU[i]*dHdx[i] + tH[i]*dUdx[i] + tV[i]*dHdy[i] + tH[i]*dVdy[i]);
        }

        // update U, V, and H
        for (int i = 0; i < n; i++) {
            _U[i] += (k1U[i] + 2.0*k2U[i] + 2.0*k3U[i] + k4U[i]) *dt / 6.0;
            _V[i] += (k1V[i] + 2.0*k2V[i] + 2.0*k3V[i] + k4V[i]) *dt / 6.0;
            _H[i] += (k1H[i] + 2.0*k2H[i] + 2.0*k3H[i] + k4H[i]) *dt / 6.0;
        }

        ct += dt;
    }
}

void ShallowWater::exportData(const std::string& fname) {
    std::ofstream file;
    file.open(fname);

    for (int j = 0; j < _U.n(); j++) {
        double y = j * _dy;
        for (int i = 0; i < _U.m(); i++) {
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

void ShallowWater::test() {
    setInitialConditions(1);
    std::string fname("output.txt");
    std::ofstream file;
    file.open(fname);

    for (int j = 0; j < _H.n(); j++) {
        for (int i = 0; i < _H.m(); i++) {
            _V.set(i, j, 2.0*j);
        }
    }

    GeneralMatrix dUdx(_U.m(), _U.n());
    GeneralMatrix dUdy(_U.m(), _U.n());
    GeneralMatrix dVdx(_U.m(), _U.n());
    GeneralMatrix dVdy(_U.m(), _U.n());
    GeneralMatrix dHdx(_U.m(), _U.n());
    GeneralMatrix dHdy(_U.m(), _U.n());
    _fd.loop(_U, _V, _H, dUdx, dUdy, dVdx, dVdy, dHdx, dHdy);

    for (int j = 0; j < _H.n(); j++) {
        double y = j * _dy;
        for (int i = 0; i < _H.m(); i++) {
            double x = i * _dx;
            file << x << " " << y << " ";
            file << _U.get(i, j) << " ";
            file << _V.get(i, j) << " ";
            file << dVdy.get(i, j) << std::endl;
        }
        file << std::endl;
    }

    file.close();
}
