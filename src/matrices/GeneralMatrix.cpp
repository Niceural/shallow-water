#include "../../include/matrices/GeneralMatrix.h"

//------------------------------------- constructor & destructor

GeneralMatrix::GeneralMatrix(const int m, const int n):
    _m(m), _n(n), _size(_m * _n),
    _mat(new double[_size]) {}

GeneralMatrix::~GeneralMatrix() {
    delete[] _mat;
}

//------------------------------------- private functions

inline int GeneralMatrix::_colMajToArrId(const int i, const int j) const {
    return i + j*_m;
}

//------------------------------------- element wise multiplication

void elementWiseMultiplicationLoop(const double alpha, const GeneralMatrix& A, const GeneralMatrix& B, const double beta, GeneralMatrix& C) {
    for (int i = 0; i < A.size(); i++)
        C._mat[i] = beta*C._mat[i] + alpha*A._mat[i]*B._mat[i];
}

void elementWiseMultiplicationBlas(const double alpha, const GeneralMatrix& A, const GeneralMatrix& B, const double beta, GeneralMatrix& C) {
    F77NAME(dsbmv)('u', A.size(), 1, alpha, A.getPointer(0), 1, B.getPointer(0), 1, beta, C.getPointer(0), 1);
}

//------------------------------------- setters

void GeneralMatrix::set(const int i, const int j, const double val) {
    _mat[_colMajToArrId(i, j)] = val;
}

double & GeneralMatrix::operator[](const int id) {
    return _mat[id];
}

double* GeneralMatrix::getPointer(const int inc) {
    return _mat + inc;
}

//------------------------------------- getters

int GeneralMatrix::m() const {
    return _m;
}

int GeneralMatrix::n() const {
    return _n;
}

int GeneralMatrix::size() const {
    return _size;
}

double GeneralMatrix::operator[](const int id) const {
    return _mat[id];
}

double GeneralMatrix::get(const int i, const int j) const {
    return _mat[_colMajToArrId(i, j)];
}

const double* GeneralMatrix::getPointer(const int inc) const {
    return _mat + inc;
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
