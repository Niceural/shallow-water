#include "../include/CentralDifference2D.h"

CentralDifference2D::CentralDifference2D(int m, int n):
    _m(m), _n(n),
    // central difference with respect to x
    _cd_x_d(GeneralBandedMatrix(m, n, 3, 3, 7)),
    // _cd_x_t1(TriangularPackedMatrix(3, 'u', false)),
    _cd_x_t1(GeneralMatrix(3, 3)),
    // _cd_x_t2(TriangularPackedMatrix(3, 'l', false)),
    _cd_x_t2(GeneralMatrix(3, 3)),
    // derivatives
    _dUdx(GeneralMatrix(m, n)),
    _dUdy(GeneralMatrix(m, n)),
    _dVdx(GeneralMatrix(m, n)),
    _dVdy(GeneralMatrix(m, n)),
    _dHdx(GeneralMatrix(m, n)),
    _dHdy(GeneralMatrix(m, n))
{
    double a = 3.0 / 4.0;
    double b = - 3.0 / 20.0;
    double c = 1.0 / 60.0;

    // central difference with respect to x - banded matrix
    double val[] = {
        c, b, a,
        0.0,
        -a, -b, -c
    };
    for (int j = 0; j < _cd_x_d.n(); j++) {
        for (int i = 0; i < _cd_x_d.ld(); i++) {
            _cd_x_d[i + j*_cd_x_d.ld()] = val[i];
        }
    }
    // _cd_x_d.print();

    // central difference with respect to x - top right triangular matrix
    double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
    for (int i = 0; i < _cd_x_t1.size(); i++)
        _cd_x_t1[i] = t1[i];
    // _cd_x_t1.print();

    // central difference with respect to x - bottom left triangular matrix
    double t2[] = { c, b, a, 0., c, b, 0., 0., c };
    for (int i = 0; i < _cd_x_t1.size(); i++)
        _cd_x_t2[i] = t2[i];
    _cd_x_t2.print();
}

CentralDifference2D::~CentralDifference2D() {}

void CentralDifference2D::integrateWrtX(GeneralMatrix& A, GeneralMatrix& dAdx) {
    for (int i = 0; i < _n; i++) {
        F77NAME(dgbmv)('N', _m, _n, _cd_x_d.kl(), _cd_x_d.ku(), 1.0, _cd_x_d.getPointer(), _cd_x_d.ld(), A.getPointer() + i*_m, 1, 0.0, dAdx.getPointer() + i*_m, 1);
    }
    // top right triangular matrix
    for (int i = 0; i < _n; i++) {
        F77NAME(dgemv)('N', _cd_x_t1.m(), _cd_x_t1.n(), 1.0, _cd_x_t1.getPointer(), 3, A.getPointer() + _m-3 + i*_m, 1, 1.0, dAdx.getPointer() + _m-3 + i*_m, 1);
    }
    // bottom left triangular matrix
    for (int i = 0; i < _n; i++) {
        F77NAME(dgemv)('N', _cd_x_t2.m(), _cd_x_t2.n(), 1.0, _cd_x_t2.getPointer(), 3, A.getPointer() + i*_m, 1, 1.0, dAdx.getPointer() + i*_m, 1);
    }
}

// void CentralDifference2D::integrate(const double* U, const double* V, const double* H) {
//     integrateU(U);
//     integrateV(V);
//     integrateH(H);
// }

// void CentralDifference2D::integrateU(const double* U) {
//     // with respect to x
//     // with respect to y
//     for (int i = 0; i < _m; i++) { // iterate over each row
//         F77NAME(dgbmv)('N', _m, _n, 3, 3, 1.0, _cdx_d, 7, U + i*_n, _m, 0.0, _dUdy + i*_n, _m);
//     }
// }

// void CentralDifference2D::integrateV(const double* V) {
//     // with respect to x
//     for (int i = 0; i < _n; i++) { // iterate over each column
//         F77NAME(dgbmv)('N', _m, _n, 3, 3, 1.0, _cdx_d, 7, V + i*_m, 1, 0.0, _dVdx + i*_m, 1);
//     }
//     // with respect to y
//     for (int i = 0; i < _m; i++) { // iterate over each row
//         F77NAME(dgbmv)('N', _m, _n, 3, 3, 1.0, _cdx_d, 7, V + i*_n, _m, 0.0, _dVdy + i*_n, _m);
//     }
// }

// void CentralDifference2D::integrateH(const double* H) {
//     // with respect to x
//     for (int i = 0; i < _n; i++) { // iterate over each column
//         F77NAME(dgbmv)('N', _m, _n, 3, 3, 1.0, _cdx_d, 7, H + i*_m, 1, 0.0, _dHdx + i*_m, 1);
//     }
//     // with respect to y
//     for (int i = 0; i < _m; i++) { // iterate over each row
//         F77NAME(dgbmv)('N', _m, _n, 3, 3, 1.0, _cdx_d, 7, H + i*_n, _m, 0.0, _dHdy + i*_n, _m);
//     }
// }