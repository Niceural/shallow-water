#include "../include/MultiQuantityMatrix.h"

//------------------------------------- misc

MultiQuantityMatrix::MultiQuantityMatrix(const int m, const int n, const int nq):
    _m(m), _n(n), _nq(nq),
    _arr(new double[m * n * nq])
{}

MultiQuantityMatrix::~MultiQuantityMatrix() {
    delete[] _arr;
}

int MultiQuantityMatrix::_get1DId(const int i, const int j, const int q) const {
    return (i + j*_m) * _nq + q;
}

//------------------------------------- setters

void MultiQuantityMatrix::set(const int i, const int j, const int q, const double val) {
    _arr[_get1DId(i, j, q)] = val;
}

void MultiQuantityMatrix::set(const int id, const int q, const double val) {
    _arr[id * _nq + q] = val;
}

void MultiQuantityMatrix::add(const int id, const int q, const double val) {
    _arr[id * _nq + q] += val;
}

//------------------------------------- getters

int MultiQuantityMatrix::m() const {
    return _m;
}

int MultiQuantityMatrix::n() const {
    return _n;
}

int MultiQuantityMatrix::nq() const {
    return _nq;
}

int MultiQuantityMatrix::size() const {
    return _m * _n;
}

double MultiQuantityMatrix::get(const int i, const int j, const int q) const {
    return _arr[_get1DId(i, j, q)];
}

double MultiQuantityMatrix::get(const int id, const int q) const {
    return _arr[id * _nq + q];
}
