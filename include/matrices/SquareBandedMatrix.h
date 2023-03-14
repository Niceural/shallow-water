#ifndef SQUARE_BANDED_MATRIX_H
#define SQUARE_BANDED_MATRIX_H

#include <stdexcept>
#include <iostream>

class SquareBandedMatrix {
    private:
        const int _n; /// Number of columns of the matrix.
        const int _kl;
        const int _ku;
        const int _ld;
        double* _mat; /// Elements.

    public:
        SquareBandedMatrix(const int n, const int kl, const int ku, const int ld);
        ~SquareBandedMatrix();

        // setters
        double & operator[](int id);
        void set(const int diag, const int i, const double val);
        double* getPointer(const int inc);

        // getters
        double operator[](int id) const;
        int n() const;
        int ld() const;
        int kl() const;
        int ku() const;
        void print() const;
        double get(const int diag, const int i) const;
};

#endif // SQUARE_BANDED_MATRIX_H
