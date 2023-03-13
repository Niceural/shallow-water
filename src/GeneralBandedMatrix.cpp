#include "../include/GeneralBandedMatrix.h"

GeneralBandedMatrix::GeneralBandedMatrix(int m, int n, int kl, int ku, int ld):
    _m(m), _n(n), _kl(kl), _ku(ku), _ld(ld),
    _mat(new double[_ld*_n])
{}

GeneralBandedMatrix::~GeneralBandedMatrix() {
    delete[] _mat;
}

double GeneralBandedMatrix::getEl(int i, int j) {
    return 0.0;
}