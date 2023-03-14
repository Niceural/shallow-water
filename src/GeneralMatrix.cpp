#include "../include/GeneralMatrix.h"

GeneralMatrix::GeneralMatrix(const int m, const int n):
    _m(m), _n(n),
    _mat(new double[_m * _n]) 
{
    if (_m < 0 || _n < 0)
        throw std::invalid_argument("Invalid argument of general matrix.");
}

GeneralMatrix::~GeneralMatrix() {
    delete[] _mat;
}

void GeneralMatrix::elementwiseMultiplication(const double alpha, GeneralMatrix& b, const double beta, GeneralMatrix& c) {
    for (int i = 0; i < size(); i++) {
        c[i] = beta*c[i]  + alpha*_mat[i]*b[i];
    }
}

inline int GeneralMatrix::_colMajToArrId(const int i, const int j) const {
    return i + j*_m;
}

void GeneralMatrix::set(int i, int j, double val) {
    _mat[_colMajToArrId(i, j)] = val;
}

double & GeneralMatrix::operator[](int id) {
    return _mat[id];
}

int GeneralMatrix::m() const {
    return _m;
}

int GeneralMatrix::n() const {
    return _n;
}

int GeneralMatrix::size() const {
    return _m * _n;
}

void GeneralMatrix::print() const {
    std::cout << std::endl;
    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < _n; j++) {
            std::cout << _mat[i + j*_m] << ", ";
        }
        std::cout << std::endl;
    }
}

double GeneralMatrix::operator[](int id) const {
    return _mat[id];
}

double GeneralMatrix::get(int i, int j) const {
    return _mat[_colMajToArrId(i, j)];
}

double* GeneralMatrix::getPointer(const int inc) {
    return _mat+inc;
}
