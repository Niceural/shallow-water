#ifndef DOMAIN_H
#define DOMAIN_H

class Domain {
    private:
        const int _n;
        const int _m;
        double* _arr;

    public:
        Domain(const int n, const int m);
        ~Domain();

        // setters
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