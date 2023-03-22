#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include <stdexcept>
#include <iostream>
#include "../blasRoutines.h"

class GeneralMatrix {
    private:
        const int _m; /// Number of rows of the matrix.
        const int _n; /// Number of columns of the matrix.
        double* _mat; /// Column major storage of the matrix.

        inline int _2dTo1d(const int i, const int j) const {
            return i + j*_m;
        }

    public:
        GeneralMatrix(const int m, const int n);
        ~GeneralMatrix();

        // element wise multiplication
        // friend void elementWiseMultiplicationLoop(const double alpha, const GeneralMatrix& A, const GeneralMatrix& B, const double beta, GeneralMatrix& C);
        // friend void elementWiseMultiplicationBlas(const double alpha, GeneralMatrix& A, GeneralMatrix& B, const double beta, GeneralMatrix& C);

        //----------------------------- setters

        inline double& operator[](const int id) { return _mat[id]; }
        inline void set(const int i, const int j, const double val) { _mat[_2dTo1d(i, j)] = val; }
        inline double* getPointer(const int inc) { return _mat + inc; }

        //----------------------------- getters

        inline int m() const { return _m; }
        inline int n() const { return _n; }
        inline int size() const { return _n * _m; }
        inline double operator[](const int id) const { return _mat[id]; }
        inline double get(const int i, const int j) const { return _mat[_2dTo1d(i, j)]; }
        inline const double* getPointer(const int inc) const { return _mat + inc; }
        void print() const;
};

#endif // GENERAL_MATRIX_H
