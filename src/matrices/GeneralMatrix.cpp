#include "../../include/matrices/GeneralMatrix.h"

//------------------------------------- constructor & destructor

GeneralMatrix::GeneralMatrix(const int m, const int n):
    _m(m), _n(n),
    _mat(new double[m*n])
{}

GeneralMatrix::~GeneralMatrix() {
    delete[] _mat;
}

//------------------------------------- element wise multiplication

// void elementWiseMultiplicationLoop(const double alpha, const GeneralMatrix& A, const GeneralMatrix& B, const double beta, GeneralMatrix& C) {
//     for (int i = 0; i < A.size(); i++)
//         C._mat[i] = beta*C._mat[i] + alpha*A._mat[i]*B._mat[i];
// }

// void elementWiseMultiplicationBlas(const double alpha, const GeneralMatrix& A, const GeneralMatrix& B, const double beta, GeneralMatrix& C) {
//     F77NAME(dsbmv)('u', A.size(), 1, alpha, A.getPointer(0), 1, B.getPointer(0), 1, beta, C.getPointer(0), 1);
// }

//------------------------------------- getters

void GeneralMatrix::print() const {
    std::cout << std::endl;
    for (int i = 0; i < _m; i++) {
        for (int j = 0; j < _n; j++) {
            std::cout << _mat[i + j*_m] << ", ";
        }
        std::cout << std::endl;
    }
}
