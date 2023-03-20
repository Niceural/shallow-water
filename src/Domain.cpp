#include "../include/Domain.h"

Domain::Domain(const int nx, const int ny, const double dx, const double dy): 
    _nx(nx), _ny(ny), _size(nx * ny * 9),
    _dx(dx), _dy(dy),
    _arr(new double[nx * ny * 9])
{}

Domain::~Domain() {
    delete[] _arr;
}

//------------------------------------- central difference

void Domain::centralDifference() {
    _cdLoopX();
    _cdLoopY();
}

void Domain::_cdLoopX() {
    const double a = 3.0 / 4.0 / get_dx();
    const double b = - 3.0 / 20.0 / get_dx();
    const double c = 1.0 / 60.0 / get_dx();

    const int m = get_nx();
    const int n = get_ny();

    for (int j = 0; j < n; j++) {
        // U, V, and H
        for (int q = 0; q < 9; q+=3) {
            // i = 0
            set(0, j, q+1,
                -c*get(m-3,j,q) -b*get(m-2,j,q) -a*get(m-1,j,q)
                +a*get(1,j,q) +b*get(2,j,q) +c*get(3,j,q));
            // i = 1
            set(1, j, q+1,
                -c*get(m-2,j,q) -b*get(m-1,j,q) -a*get(0,j,q)
                +a*get(2,j,q) +b*get(3,j,q) + c*get(4,j,q));
            // i = 2
            set(2, j, q+1,
                -c*get(m-1,j,q) -b*get(0,j,q) -a*get(1,j,q)
                +a*get(3,j,q) +b*get(4,j,q) +c*get(5,j,q));
        }

        for (int i = 3; i < m-3; i++) {
            for (int q = 0; q < 3; q++) {
                set(i, j, q+1,
                    -c*get(i-3,j,q) -b*get(i-2,j,q) -a*get(i-1,j,q)
                    +a*get(i+1,j,q) +b*get(i+2,j,q) +c*get(i+3,j,q));
            }
        }

        for (int q = 0; q < 9; q+=3) {
            // i = nx-3
            set(m-3, j, q+1,
                -c*get(m-6,j,q) -b*get(m-5,j,q) -a*get(m-4,j,q)
                +a*get(m-2,j,q) +b*get(m-1,j,q) +c*get(0,j,q));
            // i = nx-2
            set(m-2, j, q+1,
                -c*get(m-5,j,q) -b*get(m-4,j,q) -a*get(m-3,j,q)
                +a*get(m-1,j,q) +b*get(0,j,q) +c*get(1,j,q));
            // i = nx-1
            set(m-1, j, q+1,
                -c*get(m-4,j,q) -b*get(m-3,j,q) -a*get(m-2,j,q)
                +a*get(0,j,q) +b*get(1,j,q) +c*get(2,j,q));
        }
    }
}

void Domain::_cdLoopY() {
    const double a = 3.0 / 4.0 / get_dy();
    const double b = - 3.0 / 20.0 / get_dy();
    const double c = 1.0 / 60.0 / get_dy();

    const int m = get_nx();
    const int n = get_ny();

    for (int i = 0; i < m; i++) {
        for (int q = 0; q < 9; q += 3) {
            // j = 0
            set(i, 0, q+2,
                -c*get(i,n-3,q) -b*get(i,n-2,q) -a*get(i,n-1,q)
                +a*get(i,1,q) +b*get(i,2,q) +c*get(i,3,q));
            // j = 1
            set(i, 1, q+2,
                -c*get(i,n-2,q) -b*get(i,n-1,q) -a*get(i,0,q)
                +a*get(i,2,q) +b*get(i,3,q) +c*get(i,4,q));
            // j = 2
            set(i, 2, q+2,
                -c*get(i,n-1,q) -b*get(i,0,q) -a*get(i,1,q)
                +a*get(i,3,q) +b*get(i,4,q) +c*get(i,5,q));
        }

        for (int j = 3; i < n-3; j++) {
            for (int q = 0; q < 9; q += 3) {
                set(i, j, q+2,
                    -c*get(i,j-3,q) -b*get(i,j-2,q) -a*get(i,j-1,q)
                    +a*get(i,j+1,q) +b*get(i,j+2,q) +c*get(i,j+3,q));
            }
        }

        for (int q = 0; q < 9; q += 3) {
            // j = ny-3
            set(i, n-3, q+2,
                -c*get(i,n-6,q) -b*get(i,n-5,q) -a*get(i,n-4,q)
                +a*get(i,n-2,q) +b*get(i,n-1,q) +c*get(i,0,q));
            // j = ny-2
            set(i, n-2, q+2,
                -c*get(i,n-5,q) -b*get(i,n-4,q) -a*get(i,n-3,q)
                +a*get(i,n-1,q) +b*get(i,0,q) + c*get(i,1,q));
            // j = ny-1
            set(i, n-1, q+2,
                -c*get(i,n-4,q) -b*get(i,n-3,q) -a*get(i,n-2,q)
                +a*get(i,0,q) +b*get(i,1,q) + c*get(i,2,q));
        }
    }
}

//------------------------------------- setters

void Domain::set(const int i, const int j, const int quantity, const double val) {
    _arr[(i + j*_nx) * 9 + quantity] = val;
}

void Domain::set_U(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9] = val;
}

void Domain::set_dUdx(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 1] = val;
}

void Domain::set_dUdy(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 2] = val;
}

void Domain::set_V(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 3] = val;
}

void Domain::set_dVdx(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 4] = val;
}

void Domain::set_dVdy(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 5] = val;
}

void Domain::set_H(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 6] = val;
}

void Domain::set_dHdx(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 7] = val;
}

void Domain::set_dHdy(const int i, const int j, const double val) {
    _arr[(i + j*_nx) * 9 + 8] = val;
}

//------------------------------------- getters

int Domain::get_nx() const {
    return _nx;
}

int Domain::get_ny() const {
    return _ny;
}

int Domain::get_size() const {
    return _size;
}

double Domain::get_dx() const {
    return _dx;
}

double Domain::get_dy() const {
    return _dy;
}

double Domain::get(const int i, const int j, const int quantity) const {
    return _arr[(i + j*_nx) * 9 + quantity];
}

double Domain::get_U(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9];
}

double Domain::get_dUdx(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 1];
}

double Domain::get_dUdy(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 2];
}

double Domain::get_V(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 3];
}

double Domain::get_dVdx(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 4];
}

double Domain::get_dVdy(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 5];
}

double Domain::get_H(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 6];
}

double Domain::get_dHdx(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 7];
}

double Domain::get_dHdy(const int i, const int j) const {
    return _arr[(i + j*_nx) * 9 + 8];
}
