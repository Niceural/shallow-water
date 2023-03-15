#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include <stdexcept>
#include <iostream>
#include "../blasRoutines.h"

class GeneralMatrix {
    private:
        const int _m; /// Number of rows of the matrix.
        const int _n; /// Number of columns of the matrix.
        const int _size; /// Number of allocated elements (m x n).
        double* _mat; /// Column major storage of the matrix.

        inline int _colMajToArrId(const int i, const int j) const;

    public:
        GeneralMatrix(const int m, const int n);
        ~GeneralMatrix();

        // element wise multiplication
        friend void elementWiseMultiplicationLoop(const double alpha, const GeneralMatrix& A, const GeneralMatrix& B, const double beta, GeneralMatrix& C);
        friend void elementWiseMultiplicationBlas(const double alpha, GeneralMatrix& A, GeneralMatrix& B, const double beta, GeneralMatrix& C);

        // setters
        double & operator[](int id);
        void set(int i, int j, double val);
        double* getPointer(const int inc);

        // getters
        int m() const;
        int n() const;
        int size() const;
        double operator[](const int id) const;
        double get(const int i, const int j) const;
        const double* getPointer(const int inc) const;
        void print() const;
};

#endif // GENERAL_MATRIX_H
