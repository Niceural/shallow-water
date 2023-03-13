#include "../include/GeneralMatrix.h"

GeneralMatrix::GeneralMatrix(int m, int n, int lda):
    _m(m), _n(n), _lda(lda), 
    _mat(new double[m*n]) 
{
    if (_m < 0 || _n < 0 || _lda < 0)
        throw std::invalid_argument("Invalid argument of general matrix.");
}

GeneralMatrix::~GeneralMatrix() {
    delete[] _mat;
}

inline int GeneralMatrix::_colMajToArrId(int i, int j) {
    return i + j*_m;
}

double GeneralMatrix::getEl(int i, int j) {
    int id = _colMajToArrId(i, j);
    return _mat[id];
}

void GeneralMatrix::setEl(int i, int j, double val) {
    int id = _colMajToArrId(i, j);
    _mat[id] = val;
}

double* GeneralMatrix::getPointer() {
    return _mat;
}
