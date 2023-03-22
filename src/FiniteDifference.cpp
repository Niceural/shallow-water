#include "../include/FiniteDifference.h"

//------------------------------------- misc

FiniteDifference::FiniteDifference(
    const int m, const int n,
    const double dx, const double dy
): 
    _m(m), _n(n),
    _dx(dx), _dy(dy),
    _cd_d(SquareBandedMatrix(m < n ? n : m, 3, 3, 7)),
    _cd_t1(GeneralMatrix(3, 3)),
    _cd_t2(GeneralMatrix(3, 3))
{
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

    // banded matrix
    for (int j = 0; j < _cd_d.n(); j++)
        for (int i = 0; i < _cd_d.ld(); i++)
            _cd_d.set(i, j, val[i]);

    // top right triangular matrix
    for (int i = 0; i < _cd_t1.size(); i++)
        _cd_t1[i] = t1[i];

    // bottom left triangular matrix
    for (int i = 0; i < _cd_t2.size(); i++)
        _cd_t2[i] = t2[i];
}

FiniteDifference::~FiniteDifference() {
}

//------------------------------------- central difference

void FiniteDifference::centralDifference(const bool loopBlas, MultiQuantityMatrix& grid) {
    if (loopBlas) {
        _centralDifferenceBlasX(grid);
        _centralDifferenceBlasY(grid);
    } else {
        _centralDifferenceLoop(grid);
        // _centralDifferenceLoopY(grid);
    }
}

//------------------------------------- central difference loop

void FiniteDifference::_centralDifferenceLoop(MultiQuantityMatrix& grid) {
    const int m = grid.m();
    const int n = grid.n();

    const double ax = 3.0 / 4.0 / _dx;
    const double bx = - 3.0 / 20.0 / _dx;
    const double cx = 1.0 / 60.0 / _dx;

    const double ay = 3.0 / 4.0 / _dy;
    const double by = - 3.0 / 20.0 / _dy;
    const double cy = 1.0 / 60.0 / _dy;

    #pragma omp parallel default(shared)
    {

    // with respect to x, boundaries
    #pragma omp for schedule(static) nowait
    for (int j = 0; j < n; j++) {
        for (int q = 0; q < 9; q+=3) {
            // i = 0
            grid.set(0, j, q+1,
                -cx*grid.get(m-3,j,q) -bx*grid.get(m-2,j,q) -ax*grid.get(m-1,j,q)
                +ax*grid.get(1,j,q) +bx*grid.get(2,j,q) +cx*grid.get(3,j,q));
            // i = 1
            grid.set(1, j, q+1,
                -cx*grid.get(m-2,j,q) -bx*grid.get(m-1,j,q) -ax*grid.get(0,j,q)
                +ax*grid.get(2,j,q) +bx*grid.get(3,j,q) + cx*grid.get(4,j,q));
            // i = 2
            grid.set(2, j, q+1,
                -cx*grid.get(m-1,j,q) -bx*grid.get(0,j,q) -ax*grid.get(1,j,q)
                +ax*grid.get(3,j,q) +bx*grid.get(4,j,q) +cx*grid.get(5,j,q));
            // i = nx-3
            grid.set(m-3, j, q+1,
                -cx*grid.get(m-6,j,q) -bx*grid.get(m-5,j,q) -ax*grid.get(m-4,j,q)
                +ax*grid.get(m-2,j,q) +bx*grid.get(m-1,j,q) +cx*grid.get(0,j,q));
            // i = nx-2
            grid.set(m-2, j, q+1,
                -cx*grid.get(m-5,j,q) -bx*grid.get(m-4,j,q) -ax*grid.get(m-3,j,q)
                +ax*grid.get(m-1,j,q) +bx*grid.get(0,j,q) +cx*grid.get(1,j,q));
            // i = nx-1
            grid.set(m-1, j, q+1,
                -cx*grid.get(m-4,j,q) -bx*grid.get(m-3,j,q) -ax*grid.get(m-2,j,q)
                +ax*grid.get(0,j,q) +bx*grid.get(1,j,q) +cx*grid.get(2,j,q));
        }
    }

    // with respect to x, boundaries independent
    #pragma omp for schedule(static) nowait
    for (int j = 0; j < n; j++) {
        for (int i = 3; i < m-3; i++) {
            for (int q = 0; q < 9; q+=3) {
                grid.set(i, j, q+1,
                    -cx*grid.get(i-3,j,q) -bx*grid.get(i-2,j,q) -ax*grid.get(i-1,j,q)
                    +ax*grid.get(i+1,j,q) +bx*grid.get(i+2,j,q) +cx*grid.get(i+3,j,q));
            }
        }
    }

    // with respect to y, boundaries
    #pragma omp for schedule(static) nowait
    for (int i = 0; i < m; i++) {
        for (int q = 0; q < 9; q += 3) {
            // j = 0
            grid.set(i, 0, q+2,
                -cy*grid.get(i,n-3,q) -by*grid.get(i,n-2,q) -ay*grid.get(i,n-1,q)
                +ay*grid.get(i,1,q) +by*grid.get(i,2,q) +cy*grid.get(i,3,q));
            // j = 1
            grid.set(i, 1, q+2,
                -cy*grid.get(i,n-2,q) -by*grid.get(i,n-1,q) -ay*grid.get(i,0,q)
                +ay*grid.get(i,2,q) +by*grid.get(i,3,q) +cy*grid.get(i,4,q));
            // j = 2
            grid.set(i, 2, q+2,
                -cy*grid.get(i,n-1,q) -by*grid.get(i,0,q) -ay*grid.get(i,1,q)
                +ay*grid.get(i,3,q) +by*grid.get(i,4,q) +cy*grid.get(i,5,q));
            // j = ny-3
            grid.set(i, n-3, q+2,
                -cy*grid.get(i,n-6,q) -by*grid.get(i,n-5,q) -ay*grid.get(i,n-4,q)
                +ay*grid.get(i,n-2,q) +by*grid.get(i,n-1,q) +cy*grid.get(i,0,q));
            // j = ny-2
            grid.set(i, n-2, q+2,
                -cy*grid.get(i,n-5,q) -by*grid.get(i,n-4,q) -ay*grid.get(i,n-3,q)
                +ay*grid.get(i,n-1,q) +by*grid.get(i,0,q) + cy*grid.get(i,1,q));
            // j = ny-1
            grid.set(i, n-1, q+2,
                -cy*grid.get(i,n-4,q) -by*grid.get(i,n-3,q) -ay*grid.get(i,n-2,q)
                +ay*grid.get(i,0,q) +by*grid.get(i,1,q) + cy*grid.get(i,2,q));
        }
    }

    // with respect to y, boundaries independant
    #pragma omp for schedule(static) nowait
    for (int i = 0; i < m; i++) {
        for (int j = 3; j < n-3; j++) {
            for (int q = 0; q < 9; q += 3) {
                grid.set(i, j, q+2,
                    -cy*grid.get(i,j-3,q) -by*grid.get(i,j-2,q) -ay*grid.get(i,j-1,q)
                    +ay*grid.get(i,j+1,q) +by*grid.get(i,j+2,q) +cy*grid.get(i,j+3,q));
            }
        }
    }

    }
}

/*
void FiniteDifference::_centralDifferenceLoopY(MultiQuantityMatrix& grid) {
    const double ay = 3.0 / 4.0 / _dy;
    const double by = - 3.0 / 20.0 / _dy;
    const double cy = 1.0 / 60.0 / _dy;

    const int m = grid.m();
    const int n = grid.n();

    #pragma omp parallel for default(shared)
    for (int i = 0; i < m; i++) {
        for (int q = 0; q < 9; q += 3) {
            // j = 0
            grid.set(i, 0, q+2,
                -cy*grid.get(i,n-3,q) -by*grid.get(i,n-2,q) -ay*grid.get(i,n-1,q)
                +ay*grid.get(i,1,q) +by*grid.get(i,2,q) +cy*grid.get(i,3,q));
            // j = 1
            grid.set(i, 1, q+2,
                -cy*grid.get(i,n-2,q) -by*grid.get(i,n-1,q) -ay*grid.get(i,0,q)
                +ay*grid.get(i,2,q) +by*grid.get(i,3,q) +cy*grid.get(i,4,q));
            // j = 2
            grid.set(i, 2, q+2,
                -cy*grid.get(i,n-1,q) -by*grid.get(i,0,q) -ay*grid.get(i,1,q)
                +ay*grid.get(i,3,q) +by*grid.get(i,4,q) +cy*grid.get(i,5,q));
        }

        for (int j = 3; j < n-3; j++) {
            for (int q = 0; q < 9; q += 3) {
                grid.set(i, j, q+2,
                    -cy*grid.get(i,j-3,q) -by*grid.get(i,j-2,q) -ay*grid.get(i,j-1,q)
                    +ay*grid.get(i,j+1,q) +by*grid.get(i,j+2,q) +cy*grid.get(i,j+3,q));
            }
        }

        for (int q = 0; q < 9; q += 3) {
            // j = ny-3
            grid.set(i, n-3, q+2,
                -cy*grid.get(i,n-6,q) -by*grid.get(i,n-5,q) -ay*grid.get(i,n-4,q)
                +ay*grid.get(i,n-2,q) +by*grid.get(i,n-1,q) +cy*grid.get(i,0,q));
            // j = ny-2
            grid.set(i, n-2, q+2,
                -cy*grid.get(i,n-5,q) -by*grid.get(i,n-4,q) -ay*grid.get(i,n-3,q)
                +ay*grid.get(i,n-1,q) +by*grid.get(i,0,q) + cy*grid.get(i,1,q));
            // j = ny-1
            grid.set(i, n-1, q+2,
                -cy*grid.get(i,n-4,q) -by*grid.get(i,n-3,q) -ay*grid.get(i,n-2,q)
                +ay*grid.get(i,0,q) +by*grid.get(i,1,q) + cy*grid.get(i,2,q));
        }
    }
}
*/

//------------------------------------- central difference blas

void FiniteDifference::_centralDifferenceBlasX(MultiQuantityMatrix& grid) {
    // banded matrix
    // for (int i = 0; i < grid.n(); i++) {
    //     F77NAME(dgbmv)(
    //         'N', _cd_d.n(), _cd_d.n(), _cd_d.kl(), _cd_d.ku(), 1.0, _cd_d.getPointer(0), _cd_d.ld(),
    //         grid.getPointer(i*A.m()*)
    //     );
    // }
}

void FiniteDifference::_centralDifferenceBlasY(MultiQuantityMatrix& grid) {
}
