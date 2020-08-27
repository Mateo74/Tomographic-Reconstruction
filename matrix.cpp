#include "matrix.h"
#include <cmath>
#include <ctime>
#include <iostream>
#include <limits>

// Parametros asociados al criterio de convergencia del metodo de la potencia
#define PM_MAX_ERROR_1 0.01
#define PM_MAX_ERROR_2 0.1

using namespace std;

static const double epsilon = numeric_limits<double>::epsilon();

void Matrix::set_row(unsigned i, const Vector& row) {
    for (unsigned j = 0; j < num_columns(); j++) {
        _mat[i][j] = row[j];
    }
}

void Matrix::set_column(unsigned j, const Vector& column) {
    for (unsigned i = 0; i < num_rows(); i++) {
        _mat[i][j] = column[i];
    }
}

void Matrix::operator/=(double d) {
    for (unsigned i = 0; i < num_rows(); i++) {
        for (unsigned j = 0; j < num_columns(); j++) {
            _mat[i][j] /= d;
        }
    }
}

/* Calcula autovalores y autovectores, termina cuando se hayan 
 * hallado todos, cuando el calculo de uno de ellos falle, cuando 
 * un autovalor resulte menor a epsilon o negativo, o cuando un 
 * autovalor sea mucho mas grande al ultimo calculado. */
void Matrix::find_eigen(vector<double>& evalues, vector<Vector>& evectors) const {
    unsigned n = _mat.size();
    Matrix aux = *this;
    double eigenval;
    Vector eigenvec(n);
    for (unsigned k = 0; k < n; k++) {
        randomize(eigenvec);
        bool failed = !aux.find_main_eigen(eigenval, eigenvec, *this);
        if ( failed or eigenval < epsilon or (evalues.size() != 0 and evalues.back()/eigenval < 0.1 ) ) {
            return;
        }

        evalues.push_back(eigenval);
        evectors.push_back(eigenvec);
        
        // Deflacion
        for (unsigned i = 0; i < n; i++) {
            for (unsigned j = 0; j < n; j++) {
                aux(i,j) = aux(i,j) - eigenval*eigenvec[i]*eigenvec[j];
            }
        }

    }
}

/* La funcion toma como parametro adicional la verdadera matriz a la 
 * cual le estamos calculando el autovalor y autovector. La idea es que 
 * queremos parar de iterar una vez que el autovector actual 'v' y el 
 * autovalor actual 'k' cumplan que ||Av - kv|| sea muy chica comparada 
 * con ||kv||, siendo A la matriz original (que no es *this) y k = v^t Av
 * Esto va a evitar que al aplicar deflacion se vayan acumulando errores,
 * y por lo tanto los primeros autovalores/vectores van a ser tan precisos 
 * como los ultimos. Es mas, gracias a esto podemos ser menos estrictos 
 * con el máximo error tolerado (PM_MAX_ERROR_1), lo cual viene bien para 
 * el tiempo de ejecucion.
 *
 * Un inconveniente es que con esto (segun observaciones) puede pasar que 
 * no se logre un error menor a PM_MAX_ERROR_1 nunca, generando que la funcion 
 * se quede iterando infinitamente. Esto se soluciona viendo si el error 
 * es igual (con epsilon de C++) al ultimo obtenido, y en ese caso se 
 * devuelve lo que se tenga hasta ese momento. Al parecer cuando sucede 
 * esto no es muy grave y en general no vuelve a suceder con las proximas 
 * llamadas a esta funcion.
 *
 * Sin embargo, si detectamos que el error no cambio pero ademas es demasiado 
 * grande, es decir, mayor a PM_MAX_ERROR_2, entonces se concluye que el 
 * calculo del autovector/valor fallo y por lo tanto la funcion devuelve falso.
 
 * El calculo del error se hace luego de una cantidad fija de iteraciones, 
 * en vez de hacerlo en todas, para acelerar un poco el proceso. */
bool Matrix::find_main_eigen(double& eigenval, Vector& eigenvec, const Matrix& original) const {
    unsigned n = _mat.size();
    eigenvec /= two_norm(eigenvec);
    Vector temp;
    unsigned k = 0;
    double old_squared_error = 999999.0;
    while (true) {
        for (unsigned i = 0; i < 25; i++) {
            eigenvec = (*this)*eigenvec;
            eigenvec /= two_norm(eigenvec);
            k++;
        }
        
        temp = original*eigenvec;
        eigenval = inner_product(eigenvec, temp);
        
        double squared_error = 0.0;
        for (unsigned i = 0; i < n; i++) {
            squared_error += (temp[i] - eigenval * eigenvec[i])*(temp[i] - eigenval * eigenvec[i]);
        }
        squared_error /= eigenval*eigenval;
        
        if (squared_error <= PM_MAX_ERROR_1*PM_MAX_ERROR_1) {
            break;
        }
        else if (fabs(squared_error - old_squared_error) < epsilon) {
            if (squared_error <= PM_MAX_ERROR_2*PM_MAX_ERROR_2) {
                break;
            }
            else {
                return false;
            }
        }
        else {
            old_squared_error = squared_error;
        }
    }
    
    return true;
}

Matrix operator-(const Matrix& A, const Matrix& B) {
    Matrix res(A.num_rows(), A.num_columns());
    for (unsigned i = 0; i < res.num_rows(); i++) {
        for (unsigned j = 0; j < res.num_columns(); j++) {
            res(i,j) = A(i,j) - B(i,j);
        }
    }
    return res;
}

Matrix operator*(double k, const Matrix& A) {
    Matrix res(A.num_rows(), A.num_columns());
    for (unsigned i = 0; i < res.num_rows(); i++) {
        for (unsigned j = 0; j < res.num_columns(); j++) {
            res(i,j) = k * A(i,j);
        }
    }
    return res;
}

Vector operator*(const Matrix& A, const Vector& x) {
    Vector res(A.num_rows(), 0.0);
    for (unsigned i = 0; i < res.size(); i++) {
        double temp = 0.0;
        for (unsigned k = 0; k < A.num_columns(); k++) {
            temp += A(i,k) * x[k];
        }
        res[i] = temp;
    }
    return res;
}