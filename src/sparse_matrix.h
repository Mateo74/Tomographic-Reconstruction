#ifndef SPARSE_MATRIX_H
#define SPARSE_MATRIX_H

#include "vector.h"

class Matrix;

class SparseMatrix
{

public:

    SparseMatrix() {}
    
    SparseMatrix(unsigned num_rows, unsigned num_columns)
    : _columns(num_columns), _num_rows(num_rows) {}
    
    unsigned num_rows() const {
        return _num_rows;
    }
    
    unsigned num_columns() const {
        return _columns.size();
    }
    
    SparseVector& get_column(unsigned j) {
        return _columns[j];
    }

    const SparseVector& get_column(unsigned j) const {
        return _columns[j];
    }
    
    /** Realiza la multiplicacion A^t*A (siendo A == *this) y 
     *  devuelve el resultado por copia. */
    Matrix get_AtA_product() const;
    
private:

    std::vector<SparseVector> _columns;
    unsigned _num_rows;
    
};

Vector operator*(const SparseMatrix& mat, const Vector& v);

struct Metrics;
std::vector<Vector> least_squares(const SparseMatrix& A, const std::vector<Vector>& bs, Metrics& metrics);

#endif