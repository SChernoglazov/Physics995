#include "hw1.h"

double interpolation(const double x, const double* points, const int N, const double dx) {
    double term, xi, xj;
    double result = 0;
    int i, j;
    for (i = 0; i < N; i++) {
        term = 1.;
        xi = i * dx;
        for (j = 0; j < N; j++) {
            if (j != i) {
                xj = j * dx;
                term = term*(x - xj) / (xi - xj);
            }
        }
        term = term*points[i];
        result = result+term;
    }
    return result;
}

// this is procedure for problem 4
double interp_coefficient(const int i, const int N) {
    double term;
    int j, x, term1, term2;
    if( N % 2 == 0)
        x = N-1;
    else
        x = N;
    term1 = term2 = 1.;
    for (j = 0; j < N; j++) {
        if (j != i) {
            term1 = term1 * (x - 2 * j);
            term2 = term2 * (2 * (i - j));
        }
    }
    term= double(term1)/double(term2);
    return term;
}

double interpolation(const double* points, const int N){
    double result,coeff,term;
    result=0;
    for (int i=0; i<N; i++){
        coeff=interp_coefficient(i,N);
        term=coeff*points[i];
        result = result+term;
    }
    return result;
}
