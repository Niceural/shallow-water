#include "../include/CentralDifference.h"

CentralDifference::CentralDifference(
    const int m, const int n,
    const double dx, const double dy
):
    _m(m), _n(n),
    _dx(dx), _dy(dy),

    _dx_d(SquareBandedMatrix(m, 3, 3, 7)),
    _dx_t1(GeneralMatrix(3, 3)),
    _dx_t2(GeneralMatrix(3, 3)),

    _dy_d(SquareBandedMatrix(n, 3, 3, 7)),
    _dy_t1(GeneralMatrix(3, 3)),
    _dy_t2(GeneralMatrix(3, 3))
{
    // #pragma omp parallel default(shared)
    // {
    // #pragma omp sections
    // {
    
    // #pragma omp section
    _generateDx();

    // #pragma omp section
    _generateDy();

    // }
    // }
}

CentralDifference::~CentralDifference() {}

void CentralDifference::_generateDx() {
    double a = 3.0 / 4.0 / _dx;
    double b = - 3.0 / 20.0 / _dx;
    double c = 1.0 / 60.0 / _dx;

    // arrays to initialize the central difference matrices
    // banded matrix
    double val[] = {
        c, b, a,
        0.0,
        -a, -b, -c
    };
    // top right triangular matrix
    double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
    // bottom left triangular matrix
    double t2[] = { c, b, a, 0., c, b, 0., 0., c };

    // #pragma omp parallel default(shared)
    // {
    // #pragma omp sections
    // {

    // banded matrix
    // #pragma omp parallel for
    for (int j = 0; j < _dx_d.n(); j++)
        for (int i = 0; i < _dx_d.ld(); i++)
            _dx_d.set(i, j, val[i]);

    // top right triangular matrix
    // #pragma omp section
    for (int i = 0; i < _dx_t1.size(); i++)
        _dx_t1[i] = t1[i];

    // bottom left triangular matrix
    // #pragma omp section
    for (int i = 0; i < _dx_t2.size(); i++)
        _dx_t2[i] = t2[i];

    // }
    // }
}

void CentralDifference::_generateDy() {
    double a = 3.0 / 4.0 / _dy;
    double b = - 3.0 / 20.0 / _dy;
    double c = 1.0 / 60.0 / _dy;

    // banded matrix
    double val[] = {
        c, b, a,
        0.0,
        -a, -b, -c
    };
    // top right triangular matrix
    double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
    // bottom left triangular matrix
    double t2[] = { c, b, a, 0., c, b, 0., 0., c };

    // #pragma omp parallel default(shared)
    // {
    // #pragma omp sections
    // {

    // banded matrix
    // #pragma omp parallel for
    for (int j = 0; j < _dy_d.n(); j++)
        for (int i = 0; i < _dy_d.ld(); i++)
            _dy_d.set(i, j, val[i]);

    // top right triangular matrix
    // #pragma omp section
    for (int i = 0; i < _dy_t1.size(); i++)
        _dy_t1[i] = t1[i];

    // bottom left triangular matrix
    // #pragma omp section
    for (int i = 0; i < _dy_t2.size(); i++)
        _dy_t2[i] = t2[i];

    // }
    // }
}

void CentralDifference::_performWrtXLoop(const GeneralMatrix& A, GeneralMatrix& dAdx) {
    double a = 3.0 / 4.0 / _dx;
    double b = - 3.0 / 20.0 / _dx;
    double c = 1.0 / 60.0 / _dx;

    // #pragma omp parallel for // num_threads(10)
    for (int j = 0; j < A.n(); j++) {
        // i = 0
        dAdx.set(0, j, -c*A.get(A.m()-3,j) -b*A.get(A.m()-2,j) -a*A.get(A.m()-1,j)
            +a*A.get(1,j) +b*A.get(2,j) + c*A.get(3,j));

        // i = 1
        dAdx.set(1, j, -c*A.get(A.m()-2,j) -b*A.get(A.m()-1,j) -a*A.get(0,j)
            +a*A.get(2,j) +b*A.get(3,j) + c*A.get(4,j));

        // i = 2
        dAdx.set(2, j, -c*A.get(A.m()-1,j) -b*A.get(0,j) -a*A.get(1,j)
            +a*A.get(3,j) +b*A.get(4,j) + c*A.get(5,j));

        for (int i = 3; i < A.m()-3; i++) {
            dAdx.set(i, j, -c*A.get(i-3,j) -b*A.get(i-2,j) -a*A.get(i-1,j)
                +a*A.get(i+1,j) +b*A.get(i+2,j) +c*A.get(i+3,j));
        }

        // i = nx-3
        dAdx.set(A.m()-3, j, -c*A.get(A.m()-6,j) -b*A.get(A.m()-5,j) -a*A.get(A.m()-4,j)
            +a*A.get(A.m()-2,j) +b*A.get(A.m()-1,j) +c*A.get(0,j));

        // i = nx-2
        dAdx.set(A.m()-2, j, -c*A.get(A.m()-5,j) -b*A.get(A.m()-4,j) -a*A.get(A.m()-3,j)
            +a*A.get(A.m()-1,j) +b*A.get(0,j) +c*A.get(1,j));

        // i = nx-1
        dAdx.set(A.m()-1, j, -c*A.get(A.m()-4,j) -b*A.get(A.m()-3,j) -a*A.get(A.m()-2,j)
            +a*A.get(0,j) +b*A.get(1,j) +c*A.get(2,j));
    }
}

void CentralDifference::_performWrtYLoop(const GeneralMatrix& A, GeneralMatrix& dAdy) {
    double a = 3.0 / 4.0 / _dy;
    double b = - 3.0 / 20.0 / _dy;
    double c = 1.0 / 60.0 / _dy;

    // #pragma omp parallel for // num_threads(10)
    for (int i = 0; i < A.m(); i++) {
        // j = 0
        dAdy.set(i, 0, -c*A.get(i,A.n()-3) -b*A.get(i,A.n()-2) -a*A.get(i,A.n()-1)
            +a*A.get(i,1) +b*A.get(i,2) +c*A.get(i,3));

        // j = 1
        dAdy.set(i, 1, -c*A.get(i,A.n()-2) -b*A.get(i,A.n()-1) -a*A.get(i,0)
            +a*A.get(i,2) +b*A.get(i,3) +c*A.get(i,4));

        // j = 2
        dAdy.set(i, 2, -c*A.get(i,A.n()-1) -b*A.get(i,0) -a*A.get(i,1)
            +a*A.get(i,3) +b*A.get(i,4) +c*A.get(i,5));

        // j not depending on bc
        for (int j = 3; j < A.n()-3; j++) {
            dAdy.set(i, j, -c*A.get(i,j-3) -b*A.get(i,j-2) -a*A.get(i,j-1)
                +a*A.get(i,j+1) +b*A.get(i,j+2) +c*A.get(i,j+3));
        }

        // j = ny-3
        dAdy.set(i, A.n()-3, -c*A.get(i,A.n()-6) -b*A.get(i,A.n()-5) -a*A.get(i,A.n()-4)
            +a*A.get(i,A.n()-2) +b*A.get(i,A.n()-1) + c*A.get(i,0));

        // j = ny-2
        dAdy.set(i, A.n()-2, -c*A.get(i,A.n()-5) -b*A.get(i,A.n()-4) -a*A.get(i,A.n()-3)
            +a*A.get(i,A.n()-1) +b*A.get(i,0) + c*A.get(i,1));

        // j = ny-1
        dAdy.set(i, A.n()-1, -c*A.get(i,A.n()-4) -b*A.get(i,A.n()-3) -a*A.get(i,A.n()-2)
            +a*A.get(i,0) +b*A.get(i,1) + c*A.get(i,2));
    }
}

void CentralDifference::_performWrtXBlas(const GeneralMatrix& A, GeneralMatrix& dAdx) {
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

void CentralDifference::_performWrtYBlas(const GeneralMatrix& A, GeneralMatrix&dAdy) {
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
            'N', _dy_t1.m(), _dy_t1.n(), 1.0, _dy_t1.getPointer(0), _dy_t1.m(),
            A.getPointer(A.m() * (A.n()-3) + i), A.m(),
            1.0, dAdy.getPointer(i), dAdy.m()
        );
    }

    // bottom left triangular matrix
    for (int i = 0; i < A.m(); i++) {
        F77NAME(dgemv)(
            'N', _dy_t2.m(), _dy_t2.n(), 1.0, _dy_t2.getPointer(0), _dy_t1.m(),
            A.getPointer(i), A.m(),
            1.0, dAdy.getPointer(A.m() * (A.n()-3) + i), dAdy.m()
        );
    }
}

void CentralDifference::performWrtX(const bool loopBlas, const GeneralMatrix& A, GeneralMatrix& dAdx) {
    if (loopBlas) {
        _performWrtXBlas(A, dAdx);
    } else {
        _performWrtXLoop(A, dAdx);
    }
}

void CentralDifference::performWrtY(const bool loopBlas, const GeneralMatrix& A, GeneralMatrix& dAdy) {
    if (loopBlas) {
        _performWrtYBlas(A, dAdy);
    } else {
        _performWrtYLoop(A, dAdy);
    }
}
