#include "vector.h"
#include "matrix.h"
#include <algorithm>
#include <cmath>
#include <cstdlib>
#include <iostream>

using namespace std;

double inner_product(const Vector& u, const Vector& v) {
    double res = 0.0;
    for (unsigned i = 0; i < u.size(); i++) {
        res += u[i]*v[i];
    }
    return res;
}

Matrix outer_product(const Vector& u, const Vector& v) {
    Matrix res(u.size(), v.size());
    for (unsigned i = 0; i < u.size(); i++) {
        for (unsigned j = 0; j < v.size(); j++) {
            res(i,j) = u[i]*v[j];
        }
    }
    return res;
}

double squared_two_norm(const Vector& x) {
    return inner_product(x,x);
}

double one_norm(const Vector& x) {
    double res = 0.0;
    for (unsigned i = 0; i < x.size(); i++) {
        res += fabs(x[i]);
    }
    return res;
}

double two_norm(const Vector& x) {
    return sqrt(squared_two_norm(x));
}

double squared_distance(const Vector& x, const Vector& y) {
    double res = 0.0;
    for (unsigned i = 0; i < x.size(); i++) {
        res += (x[i] - y[i])*(x[i] - y[i]);
    }
    return res;
}

void zero(Vector& x) {
    for (unsigned i = 0; i < x.size(); i++) {
        x[i] = 0.0;
    }
}


void operator+=(Vector& a, const Vector& b) {
    for (unsigned i = 0; i < a.size(); i++) {
        a[i] += b[i];
    }
}

void operator-=(Vector& a, const Vector& b) {
    for (unsigned i = 0; i < a.size(); i++) {
        a[i] -= b[i];
    }
}

void operator/=(Vector& x, double d) {
    for (unsigned i = 0; i < x.size(); i++) {
        x[i] /= d;
    }
}

Vector operator-(const Vector& a, const Vector& b) {
    Vector res(a.size());
    for (unsigned i = 0; i < a.size(); i++) {
        res[i] = a[i] - b[i];
    }
    return res;
}

Vector operator*(double k, const Vector& x) {
    Vector res(x.size());
    for (unsigned i = 0; i < x.size(); i++) {
        res[i] = k*x[i];
    }
    return res;
}

Vector operator/(const Vector& x, double k) {
    Vector res(x.size());
    for (unsigned i = 0; i < x.size(); i++) {
        res[i] = x[i]/k;
    }
    return res;
}


void randomize(Vector& x) {
    for (unsigned i = 0; i < x.size(); i++) {
        x[i] = (double)rand();
    }
}

double inner_product(const SparseVector& u, const SparseVector& v) {
    unsigned i = 0, j = 0;
    double res = 0.0;
    while (i < u.size() and j < v.size()) {
        if (u[i].first == v[j].first) {
            res += u[i].second * v[j].second;
            i++;
            j++;
        }
        else if (u[i].first < v[j].first) {
            i++;
        }
        else {
            j++;
        }
    }
    
    return res;
}