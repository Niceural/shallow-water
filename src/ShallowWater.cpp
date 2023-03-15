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

    // central difference with respect to x
    _cd_x_d(SquareBandedMatrix(nx, 3, 3, 7)),
    _cd_x_t1(GeneralMatrix(3, 3)),
    _cd_x_t2(GeneralMatrix(3, 3)),
    // central difference with respect to y
    _cd_y_d(SquareBandedMatrix(ny, 3, 3, 7)),
    _cd_y_t1(GeneralMatrix(3, 3)),
    _cd_y_t2(GeneralMatrix(3, 3)),

    // central difference matrices
    _dUdx(GeneralMatrix(nx, ny)),
    _dUdy(GeneralMatrix(nx, ny)),
    _dVdx(GeneralMatrix(nx, ny)),
    _dVdy(GeneralMatrix(nx, ny)),
    _dHdx(GeneralMatrix(nx, ny)),
    _dHdy(GeneralMatrix(nx, ny))
{
    if (
        _dt < 0.0 || _t < 0.0 || _t < _dt ||
        _nx < 2 || _ny < 2 ||
        _ic < 1 || 4 < _ic
    ) throw std::invalid_argument("Invalid argument.");

}

void ShallowWater::_generateCdWrtX() {
    double a = 3.0 / 4.0 / _dx;
    double b = - 3.0 / 20.0 / _dx;
    double c = 1.0 / 60.0 / _dx;

    // central difference with respect to x - banded matrix
    double val[] = {
        c, b, a,
        0.0,
        -a, -b, -c
    };
    for (int j = 0; j < _cd_x_d.n(); j++) {
        for (int i = 0; i < _cd_x_d.ld(); i++) {
            _cd_x_d.set(i, j, val[i]);
        }
    }

    // central difference with respect to x - top right triangular matrix
    double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
    for (int i = 0; i < _cd_x_t1.size(); i++)
        _cd_x_t1[i] = t1[i];

    // central difference with respect to x - bottom left triangular matrix
    double t2[] = { c, b, a, 0., c, b, 0., 0., c };
    for (int i = 0; i < _cd_x_t1.size(); i++)
        _cd_x_t2[i] = t2[i];
}

void ShallowWater::_generateCdWrtY() {
    double a = 3.0 / 4.0 / _dy;
    double b = - 3.0 / 20.0 / _dy;
    double c = 1.0 / 60.0 / _dy;

    // central difference with respect to x - banded matrix
    double val[] = {
        c, b, a,
        0.0,
        -a, -b, -c
    };
    for (int j = 0; j < _cd_y_d.n(); j++) {
        for (int i = 0; i < _cd_y_d.ld(); i++) {
            _cd_y_d.set(i, j, val[i]);
        }
    }

    // central difference with respect to x - top right triangular matrix
    double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
    for (int i = 0; i < _cd_y_t1.size(); i++)
        _cd_y_t1[i] = t1[i];

    // central difference with respect to x - bottom left triangular matrix
    double t2[] = { c, b, a, 0., c, b, 0., 0., c };
    for (int i = 0; i < _cd_y_t1.size(); i++)
        _cd_y_t2[i] = t2[i];
}

ShallowWater::~ShallowWater() {}

void ShallowWater::integrateWrtX(GeneralMatrix& A, GeneralMatrix& dAdx) {
    for (int i = 0; i < _ny; i++) {
        // std::cout << i << std::endl;
        F77NAME(dgbmv)('N', _nx, _nx, _cd_x_d.kl(), _cd_x_d.ku(), 1.0, _cd_x_d.getPointer(0), _cd_x_d.ld(), A.getPointer(i*_nx), 1, 0.0, dAdx.getPointer(i*_nx), 1);
    }
    // top right triangular matrix
    for (int i = 0; i < _ny; i++) {
        F77NAME(dgemv)('N', _cd_x_t1.m(), _cd_x_t1.n(), 1.0, _cd_x_t1.getPointer(0), 3, A.getPointer(_nx-3 + i*_nx), 1, 1.0, dAdx.getPointer(i*_nx), 1);
    }
    // bottom left triangular matrix
    for (int i = 0; i < _ny; i++) {
        F77NAME(dgemv)('N', _cd_x_t2.m(), _cd_x_t2.n(), 1.0, _cd_x_t2.getPointer(0), 3, A.getPointer(i*_nx), 1, 1.0, dAdx.getPointer(_nx-3 + i*_nx), 1);
    }
}

void ShallowWater::integrateWrtY(GeneralMatrix& A, GeneralMatrix& dAdy) {
    for (int i = 0; i < _nx; i++) {
        F77NAME(dgbmv)('N', _cd_y_d.n(), _cd_y_d.n(), _cd_y_d.kl(), _cd_y_d.ku(), 1.0, _cd_y_d.getPointer(0), _cd_y_d.ld(), A.getPointer(i), A.m(), 0.0, dAdy.getPointer(i), A.m());
    }
    // top right triangular matrix
    for (int i = 0; i < _nx; i++) {
        F77NAME(dgemv)('N', _cd_y_t1.m(), _cd_y_t1.n(), 1.0, _cd_y_t1.getPointer(0), 3, A.getPointer(A.n()-3 + i*A.n()), A.m(), 1.0, dAdy.getPointer(i*A.n()), A.m());
    }
    // bottom left triangular matrix
    for (int i = 0; i < _nx; i++) {
        F77NAME(dgemv)('N', _cd_x_t2.m(), _cd_x_t2.n(), 1.0, _cd_x_t2.getPointer(0), 3, A.getPointer(i*A.n()), A.m(), 1.0, dAdy.getPointer(A.n()-3 + i*A.n()), A.m());
    }
}

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

    GeneralMatrix k1 = GeneralMatrix(_nx, _ny);
    GeneralMatrix k2 = GeneralMatrix(_nx, _ny);
    GeneralMatrix k3 = GeneralMatrix(_nx, _ny);
    GeneralMatrix k4 = GeneralMatrix(_nx, _ny);

    double t = 0.0;
    while (t < _t) {
        // perform central difference
        integrateWrtX(_U, _dUdx);
        integrateWrtY(_U, _dUdy);

        // k1 = f(u_n, v_n, k_n)
        _U.elementwiseMultiplication(-1.0, _dUdx, 0.0, k1); // k1 = - u * dUdx
        _V.elementwiseMultiplication(-1.0, _dUdy, 1.0, k1); // k1 = k1 - v * dUdy
        F77NAME(daxpy)(k1.size(), -_g, _dHdx.getPointer(0), 1, k1.getPointer(0), 1); // k1 = k1 - g*dHdx

        // k2 = f(u_n + dt*k1/2, v_n, h_n)
        
        t += _dt;
    }
}

GeneralMatrix ShallowWater::getH() const {
    return _H;
}
