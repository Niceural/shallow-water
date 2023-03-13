#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include <stdexcept>
#include <iostream>

class GeneralMatrix {
    private:
        int _m; /// Number of rows of the matrix.
        int _n; /// Number of columns of the matrix.
        double* _mat; /// Elements.

        inline int _colMajToArrId(int i, int j);

    public:
        GeneralMatrix(int m, int n);
        ~GeneralMatrix();

        double & operator[](int id);
        void setEl(int i, int j, double val);

        int m() const;
        int n() const;
        void print() const;
        int size() const;
        double operator[](int id) const;
        double getEl(int i, int j);
        double* getPointer();

};

#endif // GENERAL_MATRIX_H