#ifndef MATRIX_H
#define MATRIX_H

#include "vector.h"

class Matrix
{
  
friend Matrix operator-(const Matrix& A, const Matrix& B);
friend Matrix operator*(double k, const Matrix& A);
friend Vector operator*(const Matrix& A, const Vector& x);

public:

    /** Crea una matriz vacia. */
    Matrix() {}

    /** Crea una matriz sin inicializar sus valores. */
    Matrix(unsigned rows, unsigned columns)
    : _mat(rows, std::vector<double>(columns)) {}
    
    unsigned num_rows() const {
        return _mat.size();
    }
    
    unsigned num_columns() const {
        return _mat[0].size();
    }
    
    /** Obtiene el elemento (i,j) de la matriz. */
    double& operator()(unsigned i, unsigned j) {
        return _mat[i][j];
    }

    /** Version const. */
    const double& operator()(unsigned i, unsigned j) const {
        return _mat[i][j];
    }
    
    void set_row(unsigned i, const Vector& row);
    
    void set_column(unsigned j, const Vector& column);
    
    void operator/=(double d);
    
    /** Devuelve los autovalores y autovectores de la matriz. */
    void find_eigen(std::vector<double>& evalues, std::vector<Vector>& evectors) const;
    
private:

    /** Devuelve el autovalor de mayor magnitud y 
     *  su autovector asociado aplicando el metodo 
     *  de la potencia. */
    bool find_main_eigen(double& eigenval, Vector& eigenvec, const Matrix& original) const;
    
    bool find_main_eigen2(double& eigenval, Vector& eigenvec) const;

    std::vector<std::vector<double> > _mat;

};

#endif