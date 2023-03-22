#include "../include/FiniteDifference.h"

FiniteDifference::FiniteDifference(
    const int m, const int n,
    const double dx, const double dy
):
    _clx{ -1.0/60.0/dx, 3.0/20.0/dx, -3.0/4.0/dx, 3.0/4.0/dx, -3.0/20.0/dx, 1.0/60.0/dx },
    _cly{ -1.0/60.0/dy, 3.0/20.0/dy, -3.0/4.0/dy, 3.0/4.0/dy, -3.0/20.0/dy, 1.0/60.0/dy }
    // _m(m), _n(n),
    // _dx(dx), _dy(dy),

    // _dx_d(SquareBandedMatrix(m, 3, 3, 7)),
    // _dx_t1(GeneralMatrix(3, 3)),
    // _dx_t2(GeneralMatrix(3, 3)),

    // _dy_d(SquareBandedMatrix(n, 3, 3, 7)),
    // _dy_t1(GeneralMatrix(3, 3)),
    // _dy_t2(GeneralMatrix(3, 3))
{
    // #pragma omp parallel default(shared)
    // {
    // #pragma omp sections
    // {
    
    // #pragma omp section
    // _generateDx();

    // #pragma omp section
    // _generateDy();

    // }
    // }
}

FiniteDifference::~FiniteDifference() {}

//------------------------------------- loop

void FiniteDifference::loop(
    const GeneralMatrix& U,
    const GeneralMatrix& V,
    const GeneralMatrix& H,
    GeneralMatrix& dUdx,
    GeneralMatrix& dUdy,
    GeneralMatrix& dVdx,
    GeneralMatrix& dVdy,
    GeneralMatrix& dHdx,
    GeneralMatrix& dHdy
) {
    double temp = 0.0;

    //------------- independent of boundary

    for (int j = 0; j < U.n(); j++) {
        for (int i = 0; i < U.m(); i++) {
            if (2<i && i<(U.m()-3)) {
                // dUdx
                temp=_clx[0]*U.get(i-3,j); temp+=_clx[1]*U.get(i-2,j); temp+=_clx[2]*U.get(i-1,j);
                temp+=_clx[3]*U.get(i+1,j); temp+=_clx[4]*U.get(i+2,j); temp+=_clx[5]*U.get(i+3,j);
                dUdx.set(i, j, temp);
                // dVdx
                temp=_clx[0]*V.get(i-3,j); temp+=_clx[1]*V.get(i-2,j); temp+=_clx[2]*V.get(i-1,j);
                temp+=_clx[3]*V.get(i+1,j); temp+=_clx[4]*V.get(i+2,j); temp+=_clx[5]*V.get(i+3,j);
                dVdx.set(i, j, temp);
                // dHdx
                temp=_clx[0]*H.get(i-3,j); temp+=_clx[1]*H.get(i-2,j); temp+=_clx[2]*H.get(i-1,j);
                temp+=_clx[3]*H.get(i+1,j); temp+=_clx[4]*H.get(i+2,j); temp+=_clx[5]*H.get(i+3,j);
                dHdx.set(i, j, temp);
            }

            if (2<j && j<(U.n()-3)) {
                // dUdy
                temp=_cly[0]*U.get(i,j-3); temp+=_cly[1]*U.get(i,j-2); temp+=_cly[2]*U.get(i,j-1);
                temp+=_cly[3]*U.get(i,j+1); temp+=_cly[4]*U.get(i,j+2); temp+=_cly[5]*U.get(i,j+3);
                dUdy.set(i, j, temp);
                // dVdy
                temp=_cly[0]*V.get(i,j-3); temp+=_cly[1]*V.get(i,j-2); temp+=_cly[2]*V.get(i,j-1);
                temp+=_cly[3]*V.get(i,j+1); temp+=_cly[4]*V.get(i,j+2); temp+=_cly[5]*V.get(i,j+3);
                dVdy.set(i, j, temp);
                // dHdy
                temp=_cly[0]*H.get(i,j-3); temp+=_cly[1]*H.get(i,j-2); temp+=_cly[2]*H.get(i,j-1);
                temp+=_cly[3]*H.get(i,j+1); temp+=_cly[4]*H.get(i,j+2); temp+=_cly[5]*H.get(i,j+3);
                dHdy.set(i, j, temp);
            }
        }
    }

    //------------- boundaries x

    for (int j = 0; j < U.n(); j++) {
        // U
        // i = 0
        temp=_clx[0]*U.get(U.m()-3,j); temp+=_clx[1]*U.get(U.m()-2,j); temp+=_clx[2]*U.get(U.m()-1,j);
        temp+=_clx[3]*U.get(1,j); temp+=_clx[4]*U.get(2,j); temp+=_clx[5]*U.get(3,j);
        dUdx.set(0, j, temp);
        // i = 1
        temp=_clx[0]*U.get(U.m()-2,j); temp+=_clx[1]*U.get(U.m()-1,j); temp+=_clx[2]*U.get(0,j);
        temp+=_clx[3]*U.get(2,j); temp+=_clx[4]*U.get(3,j); temp+=_clx[5]*U.get(4,j);
        dUdx.set(1, j, temp);
        // i = 2
        temp=_clx[0]*U.get(U.m()-1,j); temp+=_clx[1]*U.get(0,j); temp+=_clx[2]*U.get(1,j);
        temp+=_clx[3]*U.get(3,j); temp+=_clx[4]*U.get(4,j); temp+=_clx[5]*U.get(5,j);
        dUdx.set(2, j, temp);
        // i = nx-3
        temp=_clx[0]*U.get(U.m()-6,j); temp+=_clx[1]*U.get(U.m()-5,j); temp+=_clx[2]*U.get(U.m()-4,j);
        temp+=_clx[3]*U.get(U.m()-2,j); temp+=_clx[4]*U.get(U.m()-1,j); temp+=_clx[5]*U.get(0,j);
        dUdx.set(U.m()-3, j, temp);
        // i = nx-2
        temp=_clx[0]*U.get(U.m()-5,j); temp+=_clx[1]*U.get(U.m()-4,j); temp+=_clx[2]*U.get(U.m()-3,j);
        temp+=_clx[3]*U.get(U.m()-1,j); temp+=_clx[4]*U.get(0,j); temp+=_clx[5]*U.get(1,j);
        dUdx.set(U.m()-2, j, temp);
        // i = nx-1
        temp=_clx[0]*U.get(U.m()-4,j); temp+=_clx[1]*U.get(U.m()-3,j); temp+=_clx[2]*U.get(U.m()-2,j);
        temp+=_clx[3]*U.get(0,j); temp+=_clx[4]*U.get(1,j); temp+=_clx[5]*U.get(2,j);
        dUdx.set(U.m()-1, j, temp);

        // V
        // i = 0
        temp=_clx[0]*V.get(V.m()-3,j); temp+=_clx[1]*V.get(V.m()-2,j); temp+=_clx[2]*V.get(V.m()-1,j);
        temp+=_clx[3]*V.get(1,j); temp+=_clx[4]*V.get(2,j); temp+=_clx[5]*V.get(3,j);
        dVdx.set(0, j, temp);
        // i = 1
        temp=_clx[0]*V.get(V.m()-2,j); temp+=_clx[1]*V.get(V.m()-1,j); temp+=_clx[2]*V.get(0,j);
        temp+=_clx[3]*V.get(2,j); temp+=_clx[4]*V.get(3,j); temp+=_clx[5]*V.get(4,j);
        dVdx.set(1, j, temp);
        // i = 2
        temp=_clx[0]*V.get(V.m()-1,j); temp+=_clx[1]*V.get(0,j); temp+=_clx[2]*V.get(1,j);
        temp+=_clx[3]*V.get(3,j); temp+=_clx[4]*V.get(4,j); temp+=_clx[5]*V.get(5,j);
        dVdx.set(2, j, temp);
        // i = nx-3
        temp=_clx[0]*V.get(V.m()-6,j); temp+=_clx[1]*V.get(V.m()-5,j); temp+=_clx[2]*V.get(V.m()-4,j);
        temp+=_clx[3]*V.get(V.m()-2,j); temp+=_clx[4]*V.get(V.m()-1,j); temp+=_clx[5]*V.get(0,j);
        dVdx.set(V.m()-3, j, temp);
        // i = nx-2
        temp=_clx[0]*V.get(V.m()-5,j); temp+=_clx[1]*V.get(V.m()-4,j); temp+=_clx[2]*V.get(V.m()-3,j);
        temp+=_clx[3]*V.get(V.m()-1,j); temp+=_clx[4]*V.get(0,j); temp+=_clx[5]*V.get(1,j);
        dVdx.set(V.m()-2, j, temp);
        // i = nx-1
        temp=_clx[0]*V.get(V.m()-4,j); temp+=_clx[1]*V.get(V.m()-3,j); temp+=_clx[2]*V.get(V.m()-2,j);
        temp+=_clx[3]*V.get(0,j); temp+=_clx[4]*V.get(1,j); temp+=_clx[5]*V.get(2,j);
        dVdx.set(V.m()-1, j, temp);

        // H
        // i = 0
        temp=_clx[0]*H.get(H.m()-3,j); temp+=_clx[1]*H.get(H.m()-2,j); temp+=_clx[2]*H.get(H.m()-1,j);
        temp+=_clx[3]*H.get(1,j); temp+=_clx[4]*H.get(2,j); temp+=_clx[5]*H.get(3,j);
        dHdx.set(0, j, temp);
        // i = 1
        temp=_clx[0]*H.get(H.m()-2,j); temp+=_clx[1]*H.get(H.m()-1,j); temp+=_clx[2]*H.get(0,j);
        temp+=_clx[3]*H.get(2,j); temp+=_clx[4]*H.get(3,j); temp+=_clx[5]*H.get(4,j);
        dHdx.set(1, j, temp);
        // i = 2
        temp=_clx[0]*H.get(H.m()-1,j); temp+=_clx[1]*H.get(0,j); temp+=_clx[2]*H.get(1,j);
        temp+=_clx[3]*H.get(3,j); temp+=_clx[4]*H.get(4,j); temp+=_clx[5]*H.get(5,j);
        dHdx.set(2, j, temp);
        // i = nx-3
        temp=_clx[0]*H.get(H.m()-6,j); temp+=_clx[1]*H.get(H.m()-5,j); temp+=_clx[2]*H.get(H.m()-4,j);
        temp+=_clx[3]*H.get(H.m()-2,j); temp+=_clx[4]*H.get(H.m()-1,j); temp+=_clx[5]*H.get(0,j);
        dHdx.set(H.m()-3, j, temp);
        // i = nx-2
        temp=_clx[0]*H.get(H.m()-5,j); temp+=_clx[1]*H.get(H.m()-4,j); temp+=_clx[2]*H.get(H.m()-3,j);
        temp+=_clx[3]*H.get(H.m()-1,j); temp+=_clx[4]*H.get(0,j); temp+=_clx[5]*H.get(1,j);
        dHdx.set(H.m()-2, j, temp);
        // i = nx-1
        temp=_clx[0]*H.get(H.m()-4,j); temp+=_clx[1]*H.get(H.m()-3,j); temp+=_clx[2]*H.get(H.m()-2,j);
        temp+=_clx[3]*H.get(0,j); temp+=_clx[4]*H.get(1,j); temp+=_clx[5]*H.get(2,j);
        dHdx.set(H.m()-1, j, temp);
    }

    //------------- boundaries y

    for (int i = 0; i < U.m(); i++) {
        // U
        // j = 0
        temp=_cly[0]*U.get(i,U.n()-3); temp+=_cly[1]*U.get(i,U.n()-2); temp+=_cly[2]*U.get(i,U.n()-1);
        temp+=_cly[3]*U.get(i,1); temp+=_cly[4]*U.get(i,2); temp+=_cly[5]*U.get(i,3);
        dUdy.set(i, 0, temp);
        // j = 1
        temp=_cly[0]*U.get(i,U.n()-2); temp+=_cly[1]*U.get(i,U.n()-1); temp+=_cly[2]*U.get(i,0);
        temp+=_cly[3]*U.get(i,2); temp+=_cly[4]*U.get(i,3); temp+=_cly[5]*U.get(i,4);
        dUdy.set(i, 1, temp);
        // j = 2
        temp=_cly[0]*U.get(i,U.n()-1); temp+=_cly[1]*U.get(i,0); temp+=_cly[2]*U.get(i,1);
        temp+=_cly[3]*U.get(i,3); temp+=_cly[4]*U.get(i,4); temp+=_cly[5]*U.get(i,5);
        dUdy.set(i, 2, temp);
        // j = ny-3
        temp=_cly[0]*U.get(i,U.n()-6); temp+=_cly[1]*U.get(i,U.n()-5); temp+=_cly[2]*U.get(i,U.n()-4);
        temp+=_cly[3]*U.get(i,U.n()-2); temp+=_cly[4]*U.get(i,U.n()-1); temp+=_cly[5]*U.get(i,0);
        dUdy.set(i, U.n()-3, temp);
        // j = ny-2
        temp=_cly[0]*U.get(i,U.n()-5); temp+=_cly[1]*U.get(i,U.n()-4); temp+=_cly[2]*U.get(i,U.n()-3);
        temp+=_cly[3]*U.get(i,U.n()-1); temp+=_cly[4]*U.get(i,0); temp+=_cly[5]*U.get(i,1);
        dUdy.set(i, U.n()-2, temp);
        // j = ny-1
        temp=_cly[0]*U.get(i,U.n()-4); temp+=_cly[1]*U.get(i,U.n()-3); temp+=_cly[2]*U.get(i,U.n()-2);
        temp+=_cly[3]*U.get(i,0); temp+=_cly[4]*U.get(i,1); temp+=_cly[5]*U.get(i,2);
        dUdy.set(i, U.n()-1, temp);

        // V
        // j = 0
        temp=_cly[0]*V.get(i,V.n()-3); temp+=_cly[1]*V.get(i,V.n()-2); temp+=_cly[2]*V.get(i,V.n()-1);
        temp+=_cly[3]*V.get(i,1); temp+=_cly[4]*V.get(i,2); temp+=_cly[5]*V.get(i,3);
        dVdy.set(i, 0, temp);
        // j = 1
        temp=_cly[0]*V.get(i,V.n()-2); temp+=_cly[1]*V.get(i,V.n()-1); temp+=_cly[2]*V.get(i,0);
        temp+=_cly[3]*V.get(i,2); temp+=_cly[4]*V.get(i,3); temp+=_cly[5]*V.get(i,4);
        dVdy.set(i, 1, temp);
        // j = 2
        temp=_cly[0]*V.get(i,V.n()-1); temp+=_cly[1]*V.get(i,0); temp+=_cly[2]*V.get(i,1);
        temp+=_cly[3]*V.get(i,3); temp+=_cly[4]*V.get(i,4); temp+=_cly[5]*V.get(i,5);
        dVdy.set(i, 2, temp);
        // j = ny-3
        temp=_cly[0]*V.get(i,V.n()-6); temp+=_cly[1]*V.get(i,V.n()-5); temp+=_cly[2]*V.get(i,V.n()-4);
        temp+=_cly[3]*V.get(i,V.n()-2); temp+=_cly[4]*V.get(i,V.n()-1); temp+=_cly[5]*V.get(i,0);
        dVdy.set(i, V.n()-3, temp);
        // j = ny-2
        temp=_cly[0]*V.get(i,V.n()-5); temp+=_cly[1]*V.get(i,V.n()-4); temp+=_cly[2]*V.get(i,V.n()-3);
        temp+=_cly[3]*V.get(i,V.n()-1); temp+=_cly[4]*V.get(i,0); temp+=_cly[5]*V.get(i,1);
        dVdy.set(i, V.n()-2, temp);
        // j = ny-1
        temp=_cly[0]*V.get(i,V.n()-4); temp+=_cly[1]*V.get(i,V.n()-3); temp+=_cly[2]*V.get(i,V.n()-2);
        temp+=_cly[3]*V.get(i,0); temp+=_cly[4]*V.get(i,1); temp+=_cly[5]*V.get(i,2);
        dVdy.set(i, V.n()-1, temp);

        // H
        // j = 0
        temp=_cly[0]*H.get(i,H.n()-3); temp+=_cly[1]*H.get(i,H.n()-2); temp+=_cly[2]*H.get(i,H.n()-1);
        temp+=_cly[3]*H.get(i,1); temp+=_cly[4]*H.get(i,2); temp+=_cly[5]*H.get(i,3);
        dHdy.set(i, 0, temp);
        // j = 1
        temp=_cly[0]*H.get(i,H.n()-2); temp+=_cly[1]*H.get(i,H.n()-1); temp+=_cly[2]*H.get(i,0);
        temp+=_cly[3]*H.get(i,2); temp+=_cly[4]*H.get(i,3); temp+=_cly[5]*H.get(i,4);
        dHdy.set(i, 1, temp);
        // j = 2
        temp=_cly[0]*H.get(i,H.n()-1); temp+=_cly[1]*H.get(i,0); temp+=_cly[2]*H.get(i,1);
        temp+=_cly[3]*H.get(i,3); temp+=_cly[4]*H.get(i,4); temp+=_cly[5]*H.get(i,5);
        dHdy.set(i, 2, temp);
        // j = ny-3
        temp=_cly[0]*H.get(i,H.n()-6); temp+=_cly[1]*H.get(i,H.n()-5); temp+=_cly[2]*H.get(i,H.n()-4);
        temp+=_cly[3]*H.get(i,H.n()-2); temp+=_cly[4]*H.get(i,H.n()-1); temp+=_cly[5]*H.get(i,0);
        dHdy.set(i, H.n()-3, temp);
        // j = ny-2
        temp=_cly[0]*H.get(i,H.n()-5); temp+=_cly[1]*H.get(i,H.n()-4); temp+=_cly[2]*H.get(i,H.n()-3);
        temp+=_cly[3]*H.get(i,H.n()-1); temp+=_cly[4]*H.get(i,0); temp+=_cly[5]*H.get(i,1);
        dHdy.set(i, H.n()-2, temp);
        // j = ny-1
        temp=_cly[0]*H.get(i,H.n()-4); temp+=_cly[1]*H.get(i,H.n()-3); temp+=_cly[2]*H.get(i,H.n()-2);
        temp+=_cly[3]*H.get(i,0); temp+=_cly[4]*H.get(i,1); temp+=_cly[5]*H.get(i,2);
        dHdy.set(i, H.n()-1, temp);
    }
}

//------------------------------------- loop

void FiniteDifference::blas(
    const GeneralMatrix& U,
    const GeneralMatrix& V,
    const GeneralMatrix& H,
    GeneralMatrix& dUdx,
    GeneralMatrix& dUdy,
    GeneralMatrix& dVdx,
    GeneralMatrix& dVdy,
    GeneralMatrix& dHdx,
    GeneralMatrix& dHdy
) {

}

// void FiniteDifference::_generateDx() {
//     double a = 3.0 / 4.0 / _dx;
//     double b = - 3.0 / 20.0 / _dx;
//     double c = 1.0 / 60.0 / _dx;

//     // arrays to initialize the central difference matrices
//     // banded matrix
//     double val[] = {
//         c, b, a,
//         0.0,
//         -a, -b, -c
//     };
//     // top right triangular matrix
//     double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
//     // bottom left triangular matrix
//     double t2[] = { c, b, a, 0., c, b, 0., 0., c };

//     // #pragma omp parallel default(shared)
//     // {
//     // #pragma omp sections
//     // {

//     // banded matrix
//     // #pragma omp parallel for
//     for (int j = 0; j < _dx_d.n(); j++)
//         for (int i = 0; i < _dx_d.ld(); i++)
//             _dx_d.set(i, j, val[i]);

//     // top right triangular matrix
//     // #pragma omp section
//     for (int i = 0; i < _dx_t1.size(); i++)
//         _dx_t1[i] = t1[i];

//     // bottom left triangular matrix
//     // #pragma omp section
//     for (int i = 0; i < _dx_t2.size(); i++)
//         _dx_t2[i] = t2[i];

//     // }
//     // }
// }

// void FiniteDifference::_generateDy() {
//     double a = 3.0 / 4.0 / _dy;
//     double b = - 3.0 / 20.0 / _dy;
//     double c = 1.0 / 60.0 / _dy;

//     // banded matrix
//     double val[] = {
//         c, b, a,
//         0.0,
//         -a, -b, -c
//     };
//     // top right triangular matrix
//     double t1[] = { -c, 0., 0. , -b, -c, 0., -a, -b, -c };
//     // bottom left triangular matrix
//     double t2[] = { c, b, a, 0., c, b, 0., 0., c };

//     // #pragma omp parallel default(shared)
//     // {
//     // #pragma omp sections
//     // {

//     // banded matrix
//     // #pragma omp parallel for
//     for (int j = 0; j < _dy_d.n(); j++)
//         for (int i = 0; i < _dy_d.ld(); i++)
//             _dy_d.set(i, j, val[i]);

//     // top right triangular matrix
//     // #pragma omp section
//     for (int i = 0; i < _dy_t1.size(); i++)
//         _dy_t1[i] = t1[i];

//     // bottom left triangular matrix
//     // #pragma omp section
//     for (int i = 0; i < _dy_t2.size(); i++)
//         _dy_t2[i] = t2[i];

//     // }
//     // }
// }

// void FiniteDifference::_performWrtXLoop(const GeneralMatrix& A, GeneralMatrix& dAdx) {
//     double a = 3.0 / 4.0 / _dx;
//     double b = - 3.0 / 20.0 / _dx;
//     double c = 1.0 / 60.0 / _dx;

//     // #pragma omp parallel for // num_threads(10)
//     for (int j = 0; j < A.n(); j++) {
//         // i = 0
//         dAdx.set(0, j, -c*A.get(A.m()-3,j) -b*A.get(A.m()-2,j) -a*A.get(A.m()-1,j)
//             +a*A.get(1,j) +b*A.get(2,j) + c*A.get(3,j));

//         // i = 1
//         dAdx.set(1, j, -c*A.get(A.m()-2,j) -b*A.get(A.m()-1,j) -a*A.get(0,j)
//             +a*A.get(2,j) +b*A.get(3,j) + c*A.get(4,j));

//         // i = 2
//         dAdx.set(2, j, -c*A.get(A.m()-1,j) -b*A.get(0,j) -a*A.get(1,j)
//             +a*A.get(3,j) +b*A.get(4,j) + c*A.get(5,j));

//         for (int i = 3; i < A.m()-3; i++) {
//             dAdx.set(i, j, -c*A.get(i-3,j) -b*A.get(i-2,j) -a*A.get(i-1,j)
//                 +a*A.get(i+1,j) +b*A.get(i+2,j) +c*A.get(i+3,j));
//         }

//         // i = nx-3
//         dAdx.set(A.m()-3, j, -c*A.get(A.m()-6,j) -b*A.get(A.m()-5,j) -a*A.get(A.m()-4,j)
//             +a*A.get(A.m()-2,j) +b*A.get(A.m()-1,j) +c*A.get(0,j));

//         // i = nx-2
//         dAdx.set(A.m()-2, j, -c*A.get(A.m()-5,j) -b*A.get(A.m()-4,j) -a*A.get(A.m()-3,j)
//             +a*A.get(A.m()-1,j) +b*A.get(0,j) +c*A.get(1,j));

//         // i = nx-1
//         dAdx.set(A.m()-1, j, -c*A.get(A.m()-4,j) -b*A.get(A.m()-3,j) -a*A.get(A.m()-2,j)
//             +a*A.get(0,j) +b*A.get(1,j) +c*A.get(2,j));
//     }
// }

// void FiniteDifference::_performWrtYLoop(const GeneralMatrix& A, GeneralMatrix& dAdy) {
//     double a = 3.0 / 4.0 / _dy;
//     double b = - 3.0 / 20.0 / _dy;
//     double c = 1.0 / 60.0 / _dy;

//     // #pragma omp parallel for // num_threads(10)
//     for (int i = 0; i < A.m(); i++) {
//         // j = 0
//         dAdy.set(i, 0, -c*A.get(i,A.n()-3) -b*A.get(i,A.n()-2) -a*A.get(i,A.n()-1)
//             +a*A.get(i,1) +b*A.get(i,2) +c*A.get(i,3));

//         // j = 1
//         dAdy.set(i, 1, -c*A.get(i,A.n()-2) -b*A.get(i,A.n()-1) -a*A.get(i,0)
//             +a*A.get(i,2) +b*A.get(i,3) +c*A.get(i,4));

//         // j = 2
//         dAdy.set(i, 2, -c*A.get(i,A.n()-1) -b*A.get(i,0) -a*A.get(i,1)
//             +a*A.get(i,3) +b*A.get(i,4) +c*A.get(i,5));

//         // j not depending on bc
//         for (int j = 3; j < A.n()-3; j++) {
//             dAdy.set(i, j, -c*A.get(i,j-3) -b*A.get(i,j-2) -a*A.get(i,j-1)
//                 +a*A.get(i,j+1) +b*A.get(i,j+2) +c*A.get(i,j+3));
//         }

//         // j = ny-3
//         dAdy.set(i, A.n()-3, -c*A.get(i,A.n()-6) -b*A.get(i,A.n()-5) -a*A.get(i,A.n()-4)
//             +a*A.get(i,A.n()-2) +b*A.get(i,A.n()-1) + c*A.get(i,0));

//         // j = ny-2
//         dAdy.set(i, A.n()-2, -c*A.get(i,A.n()-5) -b*A.get(i,A.n()-4) -a*A.get(i,A.n()-3)
//             +a*A.get(i,A.n()-1) +b*A.get(i,0) + c*A.get(i,1));

//         // j = ny-1
//         dAdy.set(i, A.n()-1, -c*A.get(i,A.n()-4) -b*A.get(i,A.n()-3) -a*A.get(i,A.n()-2)
//             +a*A.get(i,0) +b*A.get(i,1) + c*A.get(i,2));
//     }
// }

// void FiniteDifference::_performWrtXBlas(const GeneralMatrix& A, GeneralMatrix& dAdx) {
//     // banded matrix
//     for (int i = 0; i < A.n(); i++) {
//         F77NAME(dgbmv)(
//             'N', _dx_d.n(), _dx_d.n(), _dx_d.kl(), _dx_d.ku(), 1.0, _dx_d.getPointer(0), _dx_d.ld(),
//             A.getPointer(i*A.m()), 1,
//             0.0, dAdx.getPointer(i*dAdx.m()), 1
//         );
//     }

//     // top right triangular matrix
//     for (int i = 0; i < A.n(); i++) {
//         F77NAME(dgemv)(
//             'N', _dx_t1.m(), _dx_t1.n(), 1.0, _dx_t1.getPointer(0), 3,
//             A.getPointer(A.m()-3 + i*A.m()), 1,
//             1.0, dAdx.getPointer(i*dAdx.m()), 1
//         );
//     }

//     // bottom left triangular matrix
//     for (int i = 0; i < A.n(); i++) {
//         F77NAME(dgemv)(
//             'N', _dx_t2.m(), _dx_t2.n(), 1.0, _dx_t2.getPointer(0), 3,
//             A.getPointer(i*A.m()), 1,
//             1.0, dAdx.getPointer(dAdx.m()-3 + i*dAdx.m()), 1
//         );
//     }
// }

// void FiniteDifference::_performWrtYBlas(const GeneralMatrix& A, GeneralMatrix&dAdy) {
//     // banded matrix
//     for (int i = 0; i < A.m(); i++) {
//         F77NAME(dgbmv)(
//             'N', _dy_d.n(), _dy_d.n(), _dy_d.kl(), _dy_d.ku(), 1.0, _dy_d.getPointer(0), _dy_d.ld(),
//             A.getPointer(i), A.m(),
//             0.0, dAdy.getPointer(i), dAdy.m()
//         );
//     }

//     // top right triangular matrix
//     for (int i = 0; i < A.m(); i++) {
//         F77NAME(dgemv)(
//             'N', _dy_t1.m(), _dy_t1.n(), 1.0, _dy_t1.getPointer(0), _dy_t1.m(),
//             A.getPointer(A.m() * (A.n()-3) + i), A.m(),
//             1.0, dAdy.getPointer(i), dAdy.m()
//         );
//     }

//     // bottom left triangular matrix
//     for (int i = 0; i < A.m(); i++) {
//         F77NAME(dgemv)(
//             'N', _dy_t2.m(), _dy_t2.n(), 1.0, _dy_t2.getPointer(0), _dy_t1.m(),
//             A.getPointer(i), A.m(),
//             1.0, dAdy.getPointer(A.m() * (A.n()-3) + i), dAdy.m()
//         );
//     }
// }

// void FiniteDifference::performWrtX(const bool loopBlas, const GeneralMatrix& A, GeneralMatrix& dAdx) {
//     if (loopBlas) {
//         _performWrtXBlas(A, dAdx);
//     } else {
//         _performWrtXLoop(A, dAdx);
//     }
// }

// void FiniteDifference::performWrtY(const bool loopBlas, const GeneralMatrix& A, GeneralMatrix& dAdy) {
//     if (loopBlas) {
//         _performWrtYBlas(A, dAdy);
//     } else {
//         _performWrtYLoop(A, dAdy);
//     }
// }
