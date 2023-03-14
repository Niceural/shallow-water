#ifndef TRIANGULAR_PACKED_MATRIX_H
#define TRIANGULAR_PACKED_MATRIX_H

#include <stdexcept>

class TriangularPackedMatrix {
    private:
        int _n;
        int _size;
        char _uplo;
        bool _unitTri;
        double* _mat;
    
    public:
        TriangularPackedMatrix(int n, char uplo, bool unitTri);
        ~TriangularPackedMatrix();

        double & operator[](int id);

        int n() const;
        int size() const;
        double operator[](int id) const;
};

#endif // TRIANGULAR_PACKED_MATRIX_H