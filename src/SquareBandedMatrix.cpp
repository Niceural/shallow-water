#include "../include/SquareBandedMatrix.h"

SquareBandedMatrix::SquareBandedMatrix(
    const int n, const int kl, const int ku, const int ld
):
    _n(n), _kl(kl), _ku(ku), _ld(ld),
    _mat(new double[_ld*_n])
{}

SquareBandedMatrix::~SquareBandedMatrix() {
    delete[] _mat;
}

// setters

double & SquareBandedMatrix::operator[](int id) {
    return _mat[id];
}

void SquareBandedMatrix::set(const int diag, const int i, const double val) {
    _mat[diag + i*_ld] = val;
}

// getters

int SquareBandedMatrix::n() const {
    return _n;
}

int SquareBandedMatrix::ld() const {
    return _ld;
}

int SquareBandedMatrix::kl() const {
    return _kl;
}

int SquareBandedMatrix::ku() const {
    return _ku;
}

void SquareBandedMatrix::print() const {
    std::cout << std::endl;
    for (int i = 0; i < _ld; i++) {
        for (int j = 0; j < _n; j++) {
            std::cout << _mat[i + j*_ld] << ", ";
        }
        std::cout << std::endl;
    }
}

double SquareBandedMatrix::get(const int diag, const int i) const {
    return _mat[diag + i*_ld];
}

double* SquareBandedMatrix::getPointer(const int inc) {
    return _mat + inc;
}

double SquareBandedMatrix::operator[](int id) const {
    return _mat[id];
}
