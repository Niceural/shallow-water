#ifndef MULTI_QUANTITY_MATRIX_H
#define MULTI_QUANTITY_MATRIX_H

class MultiQuantityMatrix {
    private:
        const int _m;
        const int _n;
        const int _nq; // Number of quantities per point.
        double* _arr;

        int _get1DId(const int i, const int j, const int q) const;

    public:
        MultiQuantityMatrix(const int m, const int n, const int nq);
        ~MultiQuantityMatrix();

        // setters
        void set(const int i, const int j, const int q, const double val);
        void set(const int id, const int q, const double val);
        void add(const int id, const int q, const double val);

        // getters
        int m() const;
        int n() const;
        int nq() const;
        int size() const;
        double get(const int i, const int j, const int q) const;
        double get(const int id, const int q) const;
};

#endif // MULTI_QUANTITY_MATRIX_H