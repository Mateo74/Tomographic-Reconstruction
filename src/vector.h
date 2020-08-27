#ifndef VECTOR_H
#define VECTOR_H

#include <vector>
#include <utility>

class Matrix;

typedef std::vector<double> Vector;
typedef std::vector<std::pair<unsigned,double> > SparseVector;

double inner_product(const Vector& u, const Vector& v);
Matrix outer_product(const Vector& u, const Vector& v);
double squared_two_norm(const Vector& x);
double one_norm(const Vector& x);
double two_norm(const Vector& x);
double squared_distance(const Vector& x, const Vector& y);
void   zero(Vector& x);

void operator+=(Vector& a, const Vector& b);
void operator-=(Vector& a, const Vector& b);
void operator/=(Vector& x, double d);
Vector operator-(const Vector& a, const Vector& b);
Vector operator*(double k, const Vector& x);
Vector operator/(const Vector& x, double k);

void randomize(Vector& x);

double inner_product(const SparseVector& u, const SparseVector& v);


#endif