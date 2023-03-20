#include "../include/Domain.h"

Domain::Domain(const int n, const int m): 
    _n(n), _m(m), 
    _arr(new double[_n * _m * 9])
{}

Domain::~Domain() {
    delete[] _arr;
}

//------------------------------------- setters

void Domain::set_U(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9] = val;
}

void Domain::set_dUdx(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 1] = val;
}

void Domain::set_dUdy(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 2] = val;
}

void Domain::set_V(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 3] = val;
}

void Domain::set_dVdx(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 4] = val;
}

void Domain::set_dVdy(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 5] = val;
}

void Domain::set_H(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 6] = val;
}

void Domain::set_dHdx(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 7] = val;
}

void Domain::set_dHdy(const int i, const int j, const double val) {
    _arr[(i + j*_n) * 9 + 8] = val;
}

//------------------------------------- getters

double Domain::get_U(const int i, const int j) const {
    return _arr[(i + j*_n) * 9];
}

double Domain::get_dUdx(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 1];
}

double Domain::get_dUdy(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 2];
}

double Domain::get_V(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 3];
}

double Domain::get_dVdx(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 4];
}

double Domain::get_dVdy(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 5];
}

double Domain::get_H(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 6];
}

double Domain::get_dHdx(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 7];
}

double Domain::get_dHdy(const int i, const int j) const {
    return _arr[(i + j*_n) * 9 + 8];
}
