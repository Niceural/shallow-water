#include "../include/GeneralBandedMatrix.h"

GeneralBandedMatrix::GeneralBandedMatrix(int m, int n, int kl, int ku, int ld):
    _m(m), _n(n), _kl(kl), _ku(ku), _ld(ld),
    _mat(new double[_ld*_n])
{}

GeneralBandedMatrix::~GeneralBandedMatrix() {
    delete[] _mat;
}

double & GeneralBandedMatrix::operator[](int id) {
    return _mat[id];
}

int GeneralBandedMatrix::n() const {
    return _n;
}

int GeneralBandedMatrix::ld() const {
    return _ld;
}

int GeneralBandedMatrix::kl() const {
    return _kl;
}

int GeneralBandedMatrix::ku() const {
    return _ku;
}

void GeneralBandedMatrix::print() const {
    std::cout << std::endl;
    for (int i = 0; i < _ld; i++) {
        for (int j = 0; j < _n; j++) {
            std::cout << _mat[i + j*_ld] << ", ";
        }
        std::cout << std::endl;
    }
}

double GeneralBandedMatrix::getEl(int i, int j) {
    return _mat[i + j*_m];
}

double* GeneralBandedMatrix::getPointer() {
    return _mat;
}

double GeneralBandedMatrix::operator[](int id) const {
    return _mat[id];
}