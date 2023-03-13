#ifndef GENERAL_MATRIX_H
#define GENERAL_MATRIX_H

#include <stdexcept>

class GeneralMatrix {
    private:
        int _m; /// Number of rows of the matrix.
        int _n; /// Number of columns of the matrix.
        int _lda;
        double* _mat; /// Elements.

        inline int _colMajToArrId(int, int);

    public:
        GeneralMatrix(int, int, int);
        ~GeneralMatrix();

        double getEl(int i, int j);
        void setEl(int i, int j, double);
        double* getPointer();

};

#endif // GENERAL_MATRIX_H