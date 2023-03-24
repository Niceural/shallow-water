#include "../include/FiniteDifference.h"
#define CHUNK 8

FiniteDifference::FiniteDifference(
    const int m, const int n,
    const double dx, const double dy
):
    _cx1(-1.0/60.0/dx), _cx2(3.0/20.0/dx), _cx3(-3.0/4.0/dx), _cx4(3.0/4.0/dx), _cx5(-3.0/20.0/dx), _cx6(1.0/60.0/dx),
    _cy1(-1.0/60.0/dy), _cy2(3.0/20.0/dy), _cy3(-3.0/4.0/dy), _cy4(3.0/4.0/dy), _cy5(-3.0/20.0/dy), _cy6(1.0/60.0/dy),

    _dx_d(SquareBandedMatrix(m, 3, 3, 7)),
    _dx_t1(GeneralMatrix(3, 3)),
    _dx_t2(GeneralMatrix(3, 3)),

    _dy_d(SquareBandedMatrix(n, 3, 3, 7)),
    _dy_t1(GeneralMatrix(3, 3)),
    _dy_t2(GeneralMatrix(3, 3))
{
    _generateDx(dx);
    _generateDy(dy);
}

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
    #pragma omp parallel default(shared)
    {

    double temp = 0.0;
    //------------- independent of boundary

    // dUdx
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        for (int i = 3; i < U.m()-3; i++) {
            temp=_cx3*U.get(i-1,j);
            temp+=_cx1*U.get(i-3,j);  
            temp+=_cx2*U.get(i-2,j);
            temp+=_cx4*U.get(i+1,j);
            temp+=_cx5*U.get(i+2,j);
            temp+=_cx6*U.get(i+3,j);
            dUdx.set(i, j, temp);
        }
    }

    // dVdx
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        for (int i = 3; i < U.m()-3; i++) {
            temp=_cx3*V.get(i-1,j);
            temp+=_cx1*V.get(i-3,j);
            temp+=_cx2*V.get(i-2,j);
            temp+=_cx4*V.get(i+1,j);
            temp+=_cx5*V.get(i+2,j);
            temp+=_cx6*V.get(i+3,j);
            dVdx.set(i, j, temp);
        }
    }

    // dHdx
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        for (int i = 3; i < U.m()-3; i++) {
            temp=_cx3*H.get(i-1,j);
            temp+=_cx1*H.get(i-3,j);
            temp+=_cx2*H.get(i-2,j);
            temp+=_cx4*H.get(i+1,j);
            temp+=_cx5*H.get(i+2,j);
            temp+=_cx6*H.get(i+3,j);
            dHdx.set(i, j, temp);
        }
    }

    // dUdy
    #pragma omp for nowait schedule(static)
    for (int j = 3; j < U.n()-3; j++) {
        for (int i = 0; i < U.m(); i++) {
            temp=_cy3*U.get(i,j-1);
            temp+=_cy1*U.get(i,j-3);
            temp+=_cy2*U.get(i,j-2);
            temp+=_cy4*U.get(i,j+1);
            temp+=_cy5*U.get(i,j+2);
            temp+=_cy6*U.get(i,j+3);
            dUdy.set(i, j, temp);
        }
    }

    // dVdy
    #pragma omp for nowait schedule(static)
    for (int j = 3; j < U.n()-3; j++) {
        for (int i = 0; i < U.m(); i++) {
            temp=_cy3*V.get(i,j-1);
            temp+=_cy1*V.get(i,j-3);
            temp+=_cy2*V.get(i,j-2);
            temp+=_cy4*V.get(i,j+1);
            temp+=_cy5*V.get(i,j+2);
            temp+=_cy6*V.get(i,j+3);
            dVdy.set(i, j, temp);
        }
    }

    // dHdy
    #pragma omp for nowait schedule(static)
    for (int j = 3; j < U.n()-3; j++) {
        for (int i = 0; i < U.m(); i++) {
            temp=_cy3*H.get(i,j-1);
            temp+=_cy1*H.get(i,j-3);
            temp+=_cy2*H.get(i,j-2);
            temp+=_cy4*H.get(i,j+1);
            temp+=_cy5*H.get(i,j+2);
            temp+=_cy6*H.get(i,j+3);
            dHdy.set(i, j, temp);
        }
    }

    //------------- boundaries x

    // U
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        // i = 0
        temp=_cx1*U.get(U.m()-3,j); temp+=_cx2*U.get(U.m()-2,j); temp+=_cx3*U.get(U.m()-1,j);
        temp+=_cx4*U.get(1,j); temp+=_cx5*U.get(2,j); temp+=_cx6*U.get(3,j);
        dUdx.set(0, j, temp);
        // i = 1
        temp=_cx1*U.get(U.m()-2,j); temp+=_cx2*U.get(U.m()-1,j); temp+=_cx3*U.get(0,j);
        temp+=_cx4*U.get(2,j); temp+=_cx5*U.get(3,j); temp+=_cx6*U.get(4,j);
        dUdx.set(1, j, temp);
        // i = 2
        temp=_cx1*U.get(U.m()-1,j); temp+=_cx2*U.get(0,j); temp+=_cx3*U.get(1,j);
        temp+=_cx4*U.get(3,j); temp+=_cx5*U.get(4,j); temp+=_cx6*U.get(5,j);
        dUdx.set(2, j, temp);
    }

    // U
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        // i = nx-3
        temp=_cx1*U.get(U.m()-6,j); temp+=_cx2*U.get(U.m()-5,j); temp+=_cx3*U.get(U.m()-4,j);
        temp+=_cx4*U.get(U.m()-2,j); temp+=_cx5*U.get(U.m()-1,j); temp+=_cx6*U.get(0,j);
        dUdx.set(U.m()-3, j, temp);
        // i = nx-2
        temp=_cx1*U.get(U.m()-5,j); temp+=_cx2*U.get(U.m()-4,j); temp+=_cx3*U.get(U.m()-3,j);
        temp+=_cx4*U.get(U.m()-1,j); temp+=_cx5*U.get(0,j); temp+=_cx6*U.get(1,j);
        dUdx.set(U.m()-2, j, temp);
        // i = nx-1
        temp=_cx1*U.get(U.m()-4,j); temp+=_cx2*U.get(U.m()-3,j); temp+=_cx3*U.get(U.m()-2,j);
        temp+=_cx4*U.get(0,j); temp+=_cx5*U.get(1,j); temp+=_cx6*U.get(2,j);
        dUdx.set(U.m()-1, j, temp);
    }

    // V
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        // i = 0
        temp=_cx1*V.get(V.m()-3,j); temp+=_cx2*V.get(V.m()-2,j); temp+=_cx3*V.get(V.m()-1,j);
        temp+=_cx4*V.get(1,j); temp+=_cx5*V.get(2,j); temp+=_cx6*V.get(3,j);
        dVdx.set(0, j, temp);
        // i = 1
        temp=_cx1*V.get(V.m()-2,j); temp+=_cx2*V.get(V.m()-1,j); temp+=_cx3*V.get(0,j);
        temp+=_cx4*V.get(2,j); temp+=_cx5*V.get(3,j); temp+=_cx6*V.get(4,j);
        dVdx.set(1, j, temp);
        // i = 2
        temp=_cx1*V.get(V.m()-1,j); temp+=_cx2*V.get(0,j); temp+=_cx3*V.get(1,j);
        temp+=_cx4*V.get(3,j); temp+=_cx5*V.get(4,j); temp+=_cx6*V.get(5,j);
        dVdx.set(2, j, temp);
    }

    // V
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        // i = nx-3
        temp=_cx1*V.get(V.m()-6,j); temp+=_cx2*V.get(V.m()-5,j); temp+=_cx3*V.get(V.m()-4,j);
        temp+=_cx4*V.get(V.m()-2,j); temp+=_cx5*V.get(V.m()-1,j); temp+=_cx6*V.get(0,j);
        dVdx.set(V.m()-3, j, temp);
        // i = nx-2
        temp=_cx1*V.get(V.m()-5,j); temp+=_cx2*V.get(V.m()-4,j); temp+=_cx3*V.get(V.m()-3,j);
        temp+=_cx4*V.get(V.m()-1,j); temp+=_cx5*V.get(0,j); temp+=_cx6*V.get(1,j);
        dVdx.set(V.m()-2, j, temp);
        // i = nx-1
        temp=_cx1*V.get(V.m()-4,j); temp+=_cx2*V.get(V.m()-3,j); temp+=_cx3*V.get(V.m()-2,j);
        temp+=_cx4*V.get(0,j); temp+=_cx5*V.get(1,j); temp+=_cx6*V.get(2,j);
        dVdx.set(V.m()-1, j, temp);
    }

    // H
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        // i = 0
        temp=_cx1*H.get(H.m()-3,j); temp+=_cx2*H.get(H.m()-2,j); temp+=_cx3*H.get(H.m()-1,j);
        temp+=_cx4*H.get(1,j); temp+=_cx5*H.get(2,j); temp+=_cx6*H.get(3,j);
        dHdx.set(0, j, temp);
        // i = 1
        temp=_cx1*H.get(H.m()-2,j); temp+=_cx2*H.get(H.m()-1,j); temp+=_cx3*H.get(0,j);
        temp+=_cx4*H.get(2,j); temp+=_cx5*H.get(3,j); temp+=_cx6*H.get(4,j);
        dHdx.set(1, j, temp);
        // i = 2
        temp=_cx1*H.get(H.m()-1,j); temp+=_cx2*H.get(0,j); temp+=_cx3*H.get(1,j);
        temp+=_cx4*H.get(3,j); temp+=_cx5*H.get(4,j); temp+=_cx6*H.get(5,j);
        dHdx.set(2, j, temp);
    }

    // H
    #pragma omp for nowait schedule(static)
    for (int j = 0; j < U.n(); j++) {
        // i = nx-3
        temp=_cx1*H.get(H.m()-6,j); temp+=_cx2*H.get(H.m()-5,j); temp+=_cx3*H.get(H.m()-4,j);
        temp+=_cx4*H.get(H.m()-2,j); temp+=_cx5*H.get(H.m()-1,j); temp+=_cx6*H.get(0,j);
        dHdx.set(H.m()-3, j, temp);
        // i = nx-2
        temp=_cx1*H.get(H.m()-5,j); temp+=_cx2*H.get(H.m()-4,j); temp+=_cx3*H.get(H.m()-3,j);
        temp+=_cx4*H.get(H.m()-1,j); temp+=_cx5*H.get(0,j); temp+=_cx6*H.get(1,j);
        dHdx.set(H.m()-2, j, temp);
        // i = nx-1
        temp=_cx1*H.get(H.m()-4,j); temp+=_cx2*H.get(H.m()-3,j); temp+=_cx3*H.get(H.m()-2,j);
        temp+=_cx4*H.get(0,j); temp+=_cx5*H.get(1,j); temp+=_cx6*H.get(2,j);
        dHdx.set(H.m()-1, j, temp);
    }

    //------------- boundaries y

    // U
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.m(); i++) {
        // j = 0
        temp=_cy1*U.get(i,U.n()-3); temp+=_cy2*U.get(i,U.n()-2); temp+=_cy3*U.get(i,U.n()-1);
        temp+=_cy4*U.get(i,1); temp+=_cy5*U.get(i,2); temp+=_cy6*U.get(i,3);
        dUdy.set(i, 0, temp);
        // j = 1
        temp=_cy1*U.get(i,U.n()-2); temp+=_cy2*U.get(i,U.n()-1); temp+=_cy3*U.get(i,0);
        temp+=_cy4*U.get(i,2); temp+=_cy5*U.get(i,3); temp+=_cy6*U.get(i,4);
        dUdy.set(i, 1, temp);
        // j = 2
        temp=_cy1*U.get(i,U.n()-1); temp+=_cy2*U.get(i,0); temp+=_cy3*U.get(i,1);
        temp+=_cy4*U.get(i,3); temp+=_cy5*U.get(i,4); temp+=_cy6*U.get(i,5);
        dUdy.set(i, 2, temp);
    }

    // U
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.m(); i++) {
        // j = ny-3
        temp=_cy1*U.get(i,U.n()-6); temp+=_cy2*U.get(i,U.n()-5); temp+=_cy3*U.get(i,U.n()-4);
        temp+=_cy4*U.get(i,U.n()-2); temp+=_cy5*U.get(i,U.n()-1); temp+=_cy6*U.get(i,0);
        dUdy.set(i, U.n()-3, temp);
        // j = ny-2
        temp=_cy1*U.get(i,U.n()-5); temp+=_cy2*U.get(i,U.n()-4); temp+=_cy3*U.get(i,U.n()-3);
        temp+=_cy4*U.get(i,U.n()-1); temp+=_cy5*U.get(i,0); temp+=_cy6*U.get(i,1);
        dUdy.set(i, U.n()-2, temp);
        // j = ny-1
        temp=_cy1*U.get(i,U.n()-4); temp+=_cy2*U.get(i,U.n()-3); temp+=_cy3*U.get(i,U.n()-2);
        temp+=_cy4*U.get(i,0); temp+=_cy5*U.get(i,1); temp+=_cy6*U.get(i,2);
        dUdy.set(i, U.n()-1, temp);
    }

    // V
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.m(); i++) {
        // j = 0
        temp=_cy1*V.get(i,V.n()-3); temp+=_cy2*V.get(i,V.n()-2); temp+=_cy3*V.get(i,V.n()-1);
        temp+=_cy4*V.get(i,1); temp+=_cy5*V.get(i,2); temp+=_cy6*V.get(i,3);
        dVdy.set(i, 0, temp);
        // j = 1
        temp=_cy1*V.get(i,V.n()-2); temp+=_cy2*V.get(i,V.n()-1); temp+=_cy3*V.get(i,0);
        temp+=_cy4*V.get(i,2); temp+=_cy5*V.get(i,3); temp+=_cy6*V.get(i,4);
        dVdy.set(i, 1, temp);
        // j = 2
        temp=_cy1*V.get(i,V.n()-1); temp+=_cy2*V.get(i,0); temp+=_cy3*V.get(i,1);
        temp+=_cy4*V.get(i,3); temp+=_cy5*V.get(i,4); temp+=_cy6*V.get(i,5);
        dVdy.set(i, 2, temp);
    }

    // V
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.m(); i++) {
        // j = ny-3
        temp=_cy1*V.get(i,V.n()-6); temp+=_cy2*V.get(i,V.n()-5); temp+=_cy3*V.get(i,V.n()-4);
        temp+=_cy4*V.get(i,V.n()-2); temp+=_cy5*V.get(i,V.n()-1); temp+=_cy6*V.get(i,0);
        dVdy.set(i, V.n()-3, temp);
        // j = ny-2
        temp=_cy1*V.get(i,V.n()-5); temp+=_cy2*V.get(i,V.n()-4); temp+=_cy3*V.get(i,V.n()-3);
        temp+=_cy4*V.get(i,V.n()-1); temp+=_cy5*V.get(i,0); temp+=_cy6*V.get(i,1);
        dVdy.set(i, V.n()-2, temp);
        // j = ny-1
        temp=_cy1*V.get(i,V.n()-4); temp+=_cy2*V.get(i,V.n()-3); temp+=_cy3*V.get(i,V.n()-2);
        temp+=_cy4*V.get(i,0); temp+=_cy5*V.get(i,1); temp+=_cy6*V.get(i,2);
        dVdy.set(i, V.n()-1, temp);
    }

    // H
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.m(); i++) {
        // j = 0
        temp=_cy1*H.get(i,H.n()-3); temp+=_cy2*H.get(i,H.n()-2); temp+=_cy3*H.get(i,H.n()-1);
        temp+=_cy4*H.get(i,1); temp+=_cy5*H.get(i,2); temp+=_cy6*H.get(i,3);
        dHdy.set(i, 0, temp);
        // j = 1
        temp=_cy1*H.get(i,H.n()-2); temp+=_cy2*H.get(i,H.n()-1); temp+=_cy3*H.get(i,0);
        temp+=_cy4*H.get(i,2); temp+=_cy5*H.get(i,3); temp+=_cy6*H.get(i,4);
        dHdy.set(i, 1, temp);
        // j = 2
        temp=_cy1*H.get(i,H.n()-1); temp+=_cy2*H.get(i,0); temp+=_cy3*H.get(i,1);
        temp+=_cy4*H.get(i,3); temp+=_cy5*H.get(i,4); temp+=_cy6*H.get(i,5);
        dHdy.set(i, 2, temp);
    }

    // H
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.m(); i++) {
        // j = ny-3
        temp=_cy1*H.get(i,H.n()-6); temp+=_cy2*H.get(i,H.n()-5); temp+=_cy3*H.get(i,H.n()-4);
        temp+=_cy4*H.get(i,H.n()-2); temp+=_cy5*H.get(i,H.n()-1); temp+=_cy6*H.get(i,0);
        dHdy.set(i, H.n()-3, temp);
        // j = ny-2
        temp=_cy1*H.get(i,H.n()-5); temp+=_cy2*H.get(i,H.n()-4); temp+=_cy3*H.get(i,H.n()-3);
        temp+=_cy4*H.get(i,H.n()-1); temp+=_cy5*H.get(i,0); temp+=_cy6*H.get(i,1);
        dHdy.set(i, H.n()-2, temp);
        // j = ny-1
        temp=_cy1*H.get(i,H.n()-4); temp+=_cy2*H.get(i,H.n()-3); temp+=_cy3*H.get(i,H.n()-2);
        temp+=_cy4*H.get(i,0); temp+=_cy5*H.get(i,1); temp+=_cy6*H.get(i,2);
        dHdy.set(i, H.n()-1, temp);
    }

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
    #pragma omp parallel default(shared)
    {

    // dUdx
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.n(); i++) {
        // banded matrix
        F77NAME(dgbmv)(
            'N', _dx_d.n(), _dx_d.n(), _dx_d.kl(), _dx_d.ku(), 1.0, _dx_d.getPointer(0), _dx_d.ld(),
            U.getPointer(i*U.m()), 1,
            0.0, dUdx.getPointer(i*dUdx.m()), 1
        );
        // top right triangular matrix
        F77NAME(dgemv)(
            'N', _dx_t1.m(), _dx_t1.n(), 1.0, _dx_t1.getPointer(0), 3,
            U.getPointer(U.m()-3 + i*U.m()), 1,
            1.0, dUdx.getPointer(i*dUdx.m()), 1
        );
        // bottom left triangular matrix
        F77NAME(dgemv)(
            'N', _dx_t2.m(), _dx_t2.n(), 1.0, _dx_t2.getPointer(0), 3,
            U.getPointer(i*U.m()), 1,
            1.0, dUdx.getPointer(dUdx.m()-3 + i*dUdx.m()), 1
        );
    }

    // dVdx
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < V.n(); i++) {
        // banded matrix
        F77NAME(dgbmv)(
            'N', _dx_d.n(), _dx_d.n(), _dx_d.kl(), _dx_d.ku(), 1.0, _dx_d.getPointer(0), _dx_d.ld(),
            V.getPointer(i*V.m()), 1,
            0.0, dVdx.getPointer(i*dVdx.m()), 1
        );
        // top right triangular matrix
        F77NAME(dgemv)(
            'N', _dx_t1.m(), _dx_t1.n(), 1.0, _dx_t1.getPointer(0), 3,
            V.getPointer(V.m()-3 + i*V.m()), 1,
            1.0, dVdx.getPointer(i*dVdx.m()), 1
        );
        // bottom left triangular matrix
        F77NAME(dgemv)(
            'N', _dx_t2.m(), _dx_t2.n(), 1.0, _dx_t2.getPointer(0), 3,
            V.getPointer(i*V.m()), 1,
            1.0, dVdx.getPointer(dVdx.m()-3 + i*dVdx.m()), 1
        );
    }

    // dHdx
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < H.n(); i++) {
        // banded matrix
        F77NAME(dgbmv)(
            'N', _dx_d.n(), _dx_d.n(), _dx_d.kl(), _dx_d.ku(), 1.0, _dx_d.getPointer(0), _dx_d.ld(),
            H.getPointer(i*H.m()), 1,
            0.0, dHdx.getPointer(i*dHdx.m()), 1
        );
        // top right triangular matrix
        F77NAME(dgemv)(
            'N', _dx_t1.m(), _dx_t1.n(), 1.0, _dx_t1.getPointer(0), 3,
            H.getPointer(H.m()-3 + i*H.m()), 1,
            1.0, dHdx.getPointer(i*dHdx.m()), 1
        );
        // bottom left triangular matrix
        F77NAME(dgemv)(
            'N', _dx_t2.m(), _dx_t2.n(), 1.0, _dx_t2.getPointer(0), 3,
            H.getPointer(i*H.m()), 1,
            1.0, dHdx.getPointer(dHdx.m()-3 + i*dHdx.m()), 1
        );
    }

    // dUdy
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < U.m(); i++) {
        // banded matrix
        F77NAME(dgbmv)(
            'N', _dy_d.n(), _dy_d.n(), _dy_d.kl(), _dy_d.ku(), 1.0, _dy_d.getPointer(0), _dy_d.ld(),
            U.getPointer(i), U.m(),
            0.0, dUdy.getPointer(i), dUdy.m()
        );
        // top right triangular matrix
        F77NAME(dgemv)(
            'N', _dy_t1.m(), _dy_t1.n(), 1.0, _dy_t1.getPointer(0), _dy_t1.m(),
            U.getPointer(U.m() * (U.n()-3) + i), U.m(),
            1.0, dUdy.getPointer(i), dUdy.m()
        );
        // bottom left triangular matrix
        F77NAME(dgemv)(
            'N', _dy_t2.m(), _dy_t2.n(), 1.0, _dy_t2.getPointer(0), _dy_t1.m(),
            U.getPointer(i), U.m(),
            1.0, dUdy.getPointer(U.m() * (U.n()-3) + i), dUdy.m()
        );
    }

    // dVdy
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < V.m(); i++) {
        // banded matrix
        F77NAME(dgbmv)(
            'N', _dy_d.n(), _dy_d.n(), _dy_d.kl(), _dy_d.ku(), 1.0, _dy_d.getPointer(0), _dy_d.ld(),
            V.getPointer(i), V.m(),
            0.0, dVdy.getPointer(i), dVdy.m()
        );
        // top right triangular matrix
        F77NAME(dgemv)(
            'N', _dy_t1.m(), _dy_t1.n(), 1.0, _dy_t1.getPointer(0), _dy_t1.m(),
            V.getPointer(V.m() * (V.n()-3) + i), V.m(),
            1.0, dVdy.getPointer(i), dVdy.m()
        );
        // bottom left triangular matrix
        F77NAME(dgemv)(
            'N', _dy_t2.m(), _dy_t2.n(), 1.0, _dy_t2.getPointer(0), _dy_t1.m(),
            V.getPointer(i), V.m(),
            1.0, dVdy.getPointer(V.m() * (V.n()-3) + i), dVdy.m()
        );
    }

    // dHdy
    #pragma omp for nowait schedule(static)
    for (int i = 0; i < H.m(); i++) {
        // banded matrix
        F77NAME(dgbmv)(
            'N', _dy_d.n(), _dy_d.n(), _dy_d.kl(), _dy_d.ku(), 1.0, _dy_d.getPointer(0), _dy_d.ld(),
            H.getPointer(i), H.m(),
            0.0, dHdy.getPointer(i), dHdy.m()
        );
        // top right triangular matrix
        F77NAME(dgemv)(
            'N', _dy_t1.m(), _dy_t1.n(), 1.0, _dy_t1.getPointer(0), _dy_t1.m(),
            H.getPointer(H.m() * (H.n()-3) + i), H.m(),
            1.0, dHdy.getPointer(i), dHdy.m()
        );
        // bottom left triangular matrix
        F77NAME(dgemv)(
            'N', _dy_t2.m(), _dy_t2.n(), 1.0, _dy_t2.getPointer(0), _dy_t1.m(),
            H.getPointer(i), H.m(),
            1.0, dHdy.getPointer(H.m() * (H.n()-3) + i), dHdy.m()
        );
    }

    } // omp parallel
}

void FiniteDifference::_generateDx(const double dx) {
    double a = 3.0 / 4.0 / dx;
    double b = - 3.0 / 20.0 / dx;
    double c = 1.0 / 60.0 / dx;
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
    for (int j = 0; j < _dx_d.n(); j++)
        for (int i = 0; i < _dx_d.ld(); i++)
            _dx_d.set(i, j, val[i]);
    // top right triangular matrix
    for (int i = 0; i < _dx_t1.size(); i++)
        _dx_t1[i] = t1[i];
    // bottom left triangular matrix
    for (int i = 0; i < _dx_t2.size(); i++)
        _dx_t2[i] = t2[i];
}

void FiniteDifference::_generateDy(const double dy) {
    double a = 3.0 / 4.0 / dy;
    double b = - 3.0 / 20.0 / dy;
    double c = 1.0 / 60.0 / dy;
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
    for (int j = 0; j < _dy_d.n(); j++)
        for (int i = 0; i < _dy_d.ld(); i++)
            _dy_d.set(i, j, val[i]);
    // top right triangular matrix
    for (int i = 0; i < _dy_t1.size(); i++)
        _dy_t1[i] = t1[i];
    // bottom left triangular matrix
    for (int i = 0; i < _dy_t2.size(); i++)
        _dy_t2[i] = t2[i];
}
