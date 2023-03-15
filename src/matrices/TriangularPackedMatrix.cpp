#include "../../include/matrices/TriangularPackedMatrix.h"

TriangularPackedMatrix::TriangularPackedMatrix(int n, char uplo, bool unitTri):
    _n(n), _size(_n * (_n+1) / 2), _uplo(uplo), _unitTri(unitTri),
    _mat(new double[_size])
{
    if (
        _n < 0 ||
        (_uplo != 'u' && _uplo != 'l')
    ) throw std::invalid_argument("Invalid argument to triangular packed matrix constructor.");
}

TriangularPackedMatrix::~TriangularPackedMatrix() {
    delete[] _mat;
}

double & TriangularPackedMatrix::operator[](int id) {
    return _mat[id];
}

int TriangularPackedMatrix::n() const {
    return _n;
}

int TriangularPackedMatrix::size() const {
    return _size;
}

double TriangularPackedMatrix::operator[](int id) const {
    return _mat[id];
}