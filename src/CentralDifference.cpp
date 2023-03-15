#include "../include/CentralDifference.h"

CentralDifference::CentralDifference(
    const int m, const int n,
    const double dx, const double dy
):
    _m(n), _n(n),
    _dx_d(SquareBandedMatrix(m, 3, 3, 7)),
    _dx_t1(GeneralMatrix(3, 3)),
    _dx_t2(GeneralMatrix(3, 3)),
    _dy_d(SquareBandedMatrix(m, 3, 3, 7)),
    _dy_t1(GeneralMatrix(3, 3)),
    _dy_t2(GeneralMatrix(3, 3))
{
    _generateDx(dx);
    _generateDy(dy);
}

CentralDifference::~CentralDifference() {}

void CentralDifference::_generateDx(const double dx) {
    double a = 3.0 / 4.0 / dx;
    double b = - 3.0 / 20.0 / dx;
    double c = 1.0 / 60.0 / dx;

    // banded matrix
    double val[] = {
        c, b, a,
        0.0,
        -a, -b, -c
    };
    for (int j = 0; j < _dx_d.n(); j++)
        for (int i = 0; i < _dx_d.ld(); i++)
            _dx_d.set(i, j, val[i]);

    // top right triangular matrix
    double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
    for (int i = 0; i < _dx_t1.size(); i++)
        _dx_t1[i] = t1[i];

    // bottom left triangular matrix
    double t2[] = { c, b, a, 0., c, b, 0., 0., c };
    for (int i = 0; i < _dx_t2.size(); i++)
        _dx_t2[i] = t2[i];
}

void CentralDifference::_generateDy(const double dy) {
    double a = 3.0 / 4.0 / dy;
    double b = - 3.0 / 20.0 / dy;
    double c = 1.0 / 60.0 / dy;

    // banded matrix
    double val[] = {
        c, b, a,
        0.0,
        -a, -b, -c
    };
    for (int j = 0; j < _dy_d.n(); j++)
        for (int i = 0; i < _dy_d.ld(); i++)
            _dy_d.set(i, j, val[i]);

    // top right triangular matrix
    double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
    for (int i = 0; i < _dy_t1.size(); i++)
        _dy_t1[i] = t1[i];

    // bottom left triangular matrix
    double t2[] = { c, b, a, 0., c, b, 0., 0., c };
    for (int i = 0; i < _dy_t2.size(); i++)
        _dy_t2[i] = t2[i];

}

void CentralDifference::performWrtX(const GeneralMatrix& A, GeneralMatrix& dAdx) {
    // banded matrix
    for (int i = 0; i < A.n(); i++) {
        F77NAME(dgbmv)(
            'N', _dx_d.n(), _dx_d.n(), _dx_d.kl(), _dx_d.ku(), 1.0, _dx_d.getPointer(0), _dx_d.ld(),
            A.getPointer(i*A.m()), 1,
            0.0, dAdx.getPointer(i*dAdx.m()), 1
        );
    }

    // top right triangular matrix
    for (int i = 0; i < A.n(); i++) {
        F77NAME(dgemv)(
            'N', _dx_t1.m(), _dx_t1.n(), 1.0, _dx_t1.getPointer(0), 3,
            A.getPointer(A.m()-3 + i*A.m()), 1,
            1.0, dAdx.getPointer(i*dAdx.m()), 1
        );
    }

    // bottom left triangular matrix
    for (int i = 0; i < A.n(); i++) {
        F77NAME(dgemv)(
            'N', _dx_t2.m(), _dx_t2.n(), 1.0, _dx_t2.getPointer(0), 3,
            A.getPointer(i*A.m()), 1,
            1.0, dAdx.getPointer(dAdx.m()-3 + i*dAdx.m()), 1
        );
    }
}

void CentralDifference::performWrtY(const GeneralMatrix& A, GeneralMatrix&dAdy) {
    // banded matrix
    for (int i = 0; i < A.m(); i++) {
        F77NAME(dgbmv)(
            'N', _dy_d.n(), _dy_d.n(), _dy_d.kl(), _dy_d.ku(), 1.0, _dy_d.getPointer(0), _dy_d.ld(),
            A.getPointer(i), A.m(),
            0.0, dAdy.getPointer(i), dAdy.m()
        );
    }

    // top right triangular matrix
    for (int i = 0; i < A.m(); i++) {
        F77NAME(dgemv)(
            'N', _dy_t1.m(), _dy_t1.n(), 1.0, _dy_t1.getPointer(0), 3,
            A.getPointer(A.n()-3 + i*A.n()), A.m(),
            1.0, dAdy.getPointer(i*dAdy.n()), dAdy.m()
        );
    }

    // bottom left triangular matrix
    for (int i = 0; i < A.m(); i++) {
        F77NAME(dgemv)(
            'N', _dx_t2.m(), _dx_t2.n(), 1.0, _dx_t2.getPointer(0), 3,
            A.getPointer(i*A.n()), A.m(),
            1.0, dAdy.getPointer(dAdy.n()-3 + i*dAdy.n()), dAdy.m()
        );
    }
}