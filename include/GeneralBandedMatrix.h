#ifndef GENERAL_BANDED_MATRIX_H
#define GENERAL_BANDED_MATRIX_H

#include <stdexcept>
#include <iostream>

class GeneralBandedMatrix {
    private:
        int _m; /// Number of rows of the matrix.
        int _n; /// Number of columns of the matrix.
        int _kl;
        int _ku;
        int _ld;
        double* _mat; /// Elements.

        inline int _colMajToArrId(int, int);

    public:
        GeneralBandedMatrix(int, int, int, int, int);
        ~GeneralBandedMatrix();

        double & operator[](int id);
        void setEl(int i, int j, double);

        double operator[](int id) const;
        int n() const;
        int ld() const;
        int kl() const;
        int ku() const;
        void print() const;
        double getEl(int i, int j);
        double* getPointer();
};

#endif // GENERAL_BANDED_MATRIX_H