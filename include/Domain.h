#ifndef DOMAIN_H
#define DOMAIN_H

class Domain {
    private:
        const int _nx; const int _ny;
        const int _size;
        const double _dx; const double _dy;
        double* _arr;

        void _cdLoopX();
        void _cdLoopY();

    public:
        Domain(const int nx, const int ny, const double dx, const double dy);
        ~Domain();

        void centralDifference();

        // setters
        void set(const int i, const int j, const int quantity, const double val);
        void set_U(const int i, const int j, const double val);
        void set_dUdx(const int i, const int j, const double val);
        void set_dUdy(const int i, const int j, const double val);
        void set_V(const int i, const int j, const double val);
        void set_dVdx(const int i, const int j, const double val);
        void set_dVdy(const int i, const int j, const double val);
        void set_H(const int i, const int j, const double val);
        void set_dHdx(const int i, const int j, const double val);
        void set_dHdy(const int i, const int j, const double val);

        // getters
        int get_nx() const;
        int get_ny() const;
        int get_size() const;
        double get_dx() const;
        double get_dy() const;
        double get(const int i, const int j, const int quantity) const;
        double get_U(const int i, const int j) const;
        double get_dUdx(const int i, const int j) const;
        double get_dUdy(const int i, const int j) const;
        double get_V(const int i, const int j) const;
        double get_dVdx(const int i, const int j) const;
        double get_dVdy(const int i, const int j) const;
        double get_H(const int i, const int j) const;
        double get_dHdx(const int i, const int j) const;
        double get_dHdy(const int i, const int j) const;
};

#endif // DOMAIN_H