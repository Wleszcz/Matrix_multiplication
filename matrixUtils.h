//
// Created by Wiktor on 21.04.2023.
//
#include "Matrix.h"
#ifndef NUMERYCZNE2_MATRIXUTILS_H
#define NUMERYCZNE2_MATRIXUTILS_H

Matrix *Jacobi(Matrix* A,Matrix* b,double tol);

Matrix *GaussSeidel(Matrix *A, Matrix *p, double tol);

Matrix *forwardSubstitution(Matrix *L, Matrix *b);
Matrix* backwardSubstitution(Matrix* A, Matrix* b);
#endif //NUMERYCZNE2_MATRIXUTILS_H
