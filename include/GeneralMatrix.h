#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include <stdexcept>
#include <iostream>

class GeneralMatrix {
    private:
        const int _m; /// Number of rows of the matrix.
        const int _n; /// Number of columns of the matrix.
        double* _mat; /// Elements.

        inline int _colMajToArrId(const int i, const int j) const;

    public:
        GeneralMatrix(const int m, const int n);
        ~GeneralMatrix();

        void elementwiseMultiplication(const double alpha, GeneralMatrix& b, const double beta, GeneralMatrix& c);

        // setters
        double & operator[](int id);
        void set(int i, int j, double val);
        double* getPointer(const int inc);

        // getters
        int m() const;
        int n() const;
        void print() const;
        int size() const;
        double operator[](int id) const;
        double get(int i, int j) const;
};

#endif // GENERAL_MATRIX_H
