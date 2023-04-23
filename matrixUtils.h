//
// Created by Wiktor on 21.04.2023.
//
#include "Matrix.h"
#ifndef NUMERYCZNE2_MATRIXUTILS_H
#define NUMERYCZNE2_MATRIXUTILS_H

const int MAX_iter = 1000;
const int index[] ={1,8,8,7,3,1};


Matrix *Jacobi(Matrix* A,Matrix* b,double tol,int * time);

Matrix *GaussSeidel(Matrix *A, Matrix *p, double tol, int * time);

Matrix ** LUFactorisation(Matrix *A);
Matrix * LUFactorisationSolving(Matrix * A, Matrix* b, int * time);

void compare(Matrix* A,Matrix* B,Matrix* C,double tol);
void compare(Matrix* A,Matrix* B,double tol);
void createBandMatrix(Matrix* A, int a1, int a2, int a3);

Matrix *forwardSubstitution(Matrix *L, Matrix *b);
Matrix* backwardSubstitution(Matrix* A, Matrix* b);
#endif //NUMERYCZNE2_MATRIXUTILS_H
