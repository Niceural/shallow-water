#ifndef MULTI_QUANTITY_MATRIX_H
#define MULTI_QUANTITY_MATRIX_H

class MultiQuantityMatrix {
    private:
        const int _m;
        const int _n;
        const int _nq; // Number of quantities per point.
        double* _arr;

        inline int _get1DId(const int i, const int j, const int q) const {
            return (i + j*_m) * _nq + q;
        }

    public:
        MultiQuantityMatrix(const int m, const int n, const int nq);
        ~MultiQuantityMatrix();

        // setters
        inline void set(const int i, const int j, const int q, const double val) {
            _arr[_get1DId(i, j, q)] = val;
        }
        inline void set(const int id, const int q, const double val) {
            _arr[id * _nq + q] = val;
        }
        inline void add(const int id, const int q, const double val) {
            _arr[id * _nq + q] += val;
        }
        inline double* getPointer(const int i, const int j, const int q) {
            return _arr + _get1DId(i, j, q);
        }

        // getters
        int m() const;
        int n() const;
        int nq() const;
        int size() const;
        inline double get(const int i, const int j, const int q) const {
            return _arr[_get1DId(i, j, q)];
        }
        inline double get(const int id, const int q) const {
            return _arr[id * _nq + q];
        }
        inline const double* getPointer(const int i, const int j, const int q) const {
            return _arr + _get1DId(i, j, q);
        }
};

#endif // MULTI_QUANTITY_MATRIX_H