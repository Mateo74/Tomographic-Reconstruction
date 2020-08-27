#ifndef DEBUG_H
#define DEBUG_H

#include "vector.h"
#include "matrix.h"
#include "sparse_matrix.h"
#include <fstream>
#include <iostream>

inline void print(const Vector& vec)
{
    for (unsigned i = 0; i < vec.size(); i++) {
        std::cout << vec[i] << std::endl;
    }
    std::cout << std::endl;
}

inline void print(const Matrix& mat)
{
    for (unsigned i = 0; i < mat.num_rows(); i++) {
        for (unsigned j = 0; j < mat.num_columns(); j++) {
            std::cout << mat(i,j) << '\t';
        }
        std::cout << std::endl;
    }
}

inline void print(const std::vector<Vector>& vs)
{
    for (unsigned i = 0; i < vs.size(); i++) {
        for (unsigned j = 0; j < vs[i].size(); j++) {
            std::cout << vs[i][j] << '\t';
        }
        std::cout << std::endl;
    }
}

inline void print(const Matrix& mat, const std::string& file)
{
    std::ofstream ofile(file);
    for (unsigned i = 0; i < mat.num_rows(); i++) {
        for (unsigned j = 0; j < mat.num_columns(); j++) {
            ofile << mat(i,j) << " ";
        }
        ofile << std::endl;
    }
}

inline void print(const Vector& vec, const std::string& file)
{
    std::ofstream ofile(file);
    for (unsigned i = 0; i < vec.size(); i++) {
        ofile << vec[i] << std::endl;
    }
}

inline void print(const SparseMatrix& D, unsigned discr_size, const std::string& file)
{
    unsigned rows = D.num_rows();
    unsigned cols = D.num_columns();

    Matrix M(rows, cols);
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            M(i,j) = 0.0;
        }
    }
    for (unsigned j = 0; j < cols; j++) {
        const SparseVector& col = D.get_column(j);
        for (unsigned elem = 0; elem < col.size(); elem++) {
            M(col[elem].first, j) = col[elem].second;
        }
    }
    
    std::ofstream ofile(file);
    for (unsigned i = 0; i < M.num_rows(); i++) {
        for (unsigned j = 0; j < discr_size; j++) {
            for (unsigned k = 0; k < discr_size; k++) {
                ofile << M(i, j*discr_size + k) << " ";
            }
            ofile << std::endl;
        }
        ofile << std::endl;
    }
}

inline void print(const SparseMatrix& D, const std::string& file)
{
    unsigned rows = D.num_rows();
    unsigned cols = D.num_columns();

    Matrix M(rows, cols);
    for (unsigned i = 0; i < rows; i++) {
        for (unsigned j = 0; j < cols; j++) {
            M(i,j) = 0.0;
        }
    }
    for (unsigned j = 0; j < cols; j++) {
        const SparseVector& col = D.get_column(j);
        for (unsigned elem = 0; elem < col.size(); elem++) {
            M(col[elem].first, j) = col[elem].second;
        }
    }
    
    std::ofstream ofile(file);
    for (unsigned i = 0; i < M.num_rows(); i++) {
        for (unsigned j = 0; j < M.num_columns(); j++) {
            ofile << M(i,j) << " ";
        }
        ofile << std::endl;
    }
}


#endif