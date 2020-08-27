#include "sparse_matrix.h"
#include "matrix.h"
#include "metrics.h"
#include <algorithm>
#include <cmath>
#include <ctime>
#include <limits>

#include <iostream>

using namespace std;

static const double epsilon = numeric_limits<double>::epsilon();


Matrix SparseMatrix::get_AtA_product() const
{
    unsigned n = _columns.size();
    Matrix AtA(n, n);
    
    for (unsigned i = 0; i < n; i++) {
        for (unsigned j = 0; j <= i; j++) {
            AtA(i,j) = inner_product(_columns[i], _columns[j]);
            AtA(j,i) = AtA(i,j);
        }
    }
    
    return AtA;
}

Vector operator*(const SparseMatrix& mat, const Vector& v)
{
    Vector res(mat.num_rows(), 0.0);
    
    for (unsigned j = 0; j < mat.num_columns(); j++) {
        const SparseVector& col = mat.get_column(j);
        for (unsigned elem = 0; elem < col.size(); elem++) {
            res[col[elem].first] += col[elem].second * v[j];
        }
    }
    
    return res;
}


vector<Vector> least_squares(const SparseMatrix& A, const vector<Vector>& bs, Metrics& metrics)
{
    clock_t start = clock();
    
    unsigned m = A.num_rows(); unsigned n = A.num_columns();
    
    Matrix AtA = A.get_AtA_product();
    
    vector<double> evalues;
    vector<Vector> evectors;
    AtA.find_eigen(evalues, evectors);
    
    Vector svalues(evalues.size());
    for (unsigned i = 0; i < svalues.size(); i++) {
        svalues[i] = sqrt(evalues[i]);
    }

    Matrix Ut(svalues.size(), m);
    for (unsigned i = 0; i < svalues.size(); i++) {
        Ut.set_row(i, A*evectors[i] / svalues[i]);
    }

    Matrix V(n, svalues.size());
    for (unsigned j = 0; j < svalues.size(); j++) {
        V.set_column(j, evectors[j]);
    }
    
    vector<Vector> results(bs.size());
    
    for (unsigned k = 0; k < bs.size(); k++) {
        Vector c = Ut*bs[k];

        Vector y(svalues.size());
        for (unsigned i = 0; i < y.size(); i++) {
            y[i] = c[i] / svalues[i];
        }
        
        results[k] = V*y;
    }

    metrics.reconstruction_time = (double)(clock() - start) / CLOCKS_PER_SEC;
    metrics.cond_number = evalues[0] / evalues.back();
    metrics.num_eigen_found = evalues.size();
    
    return results;
}