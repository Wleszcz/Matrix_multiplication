//
// Created by Wiktor on 21.04.2023.
//
#include "matrixUtils.h"
#include "Matrix.h"
#include "Matrix.h"
#include "stdio.h"
#include "cstdlib"

const int index[] ={1,8,8,7,3,1};

Matrix ** LUFactorisation(Matrix *A){

    int rows = A->Y;
    int cols = A->X;



    Matrix * U = A->copy();
    Matrix * L = new Matrix(rows,cols);


    // L - macierz jednostkowa
    for (int i = 0; i < rows; i++) {
        L->Mat[i][i] = 1.0;
    }

    for (int k = 0; k < rows-1; k++) {
        for (int i = k+1; i < rows; i++) {
            double factor = U->Mat[i][k] / U->Mat[k][k];

            // eliminacja Gaussa wiersza i z wykorzystaniem wiersza k
            for (int j = k; j < cols; j++) {
                U->Mat[i][j] -= factor * U->Mat[k][j];
            }

            L->Mat[i][k] = factor;
        }
    }


    Matrix** result = new Matrix*[2];
    result[0]=L;
    result[1]=U;

    return result;
}


Matrix *GaussSeidel(Matrix *A, Matrix *b, double tol) {
    int N=b->Y;
    Matrix* x = new Matrix(N,1);
    Matrix* x_pr = new Matrix(N,1);
    x->ones();

    Matrix* D = A->Diag();
    Matrix* L  = A->LDiag();
    Matrix* U  = A->UDiag();
    Matrix* DL = D->add(L);
    double err=1;

    int iter =0;
    while(err>tol){
        x_pr->eq(x);

        //solving  A = (Diag+LDiag)^-1 (b-UDiag*x_pr) using forward substitution

        Matrix *ux=U->mul(x_pr);
        Matrix * b_ux=b->sub(ux);

        x->del();
        x = (forwardSubstitution(DL,b_ux));

        err=((A->mul(x))->sub(b))->norm();
        iter++;

        delete(ux);
        delete(b_ux);

        if(iter>1000) break;
    }
    delete(L);
    delete(D);
    delete(U);
    delete(DL);
    delete(x_pr);
    return x;

}

Matrix* forwardSubstitution(Matrix* L, Matrix* b) {
    int N = b->Y;
    Matrix* x = new Matrix(N, 1);

    for (int i = 0; i < N; i++) {
        double sum = 0.0;
        for (int j = 0; j < i; j++) {
            sum += L->Mat[i][j] * x->Mat[j][0];
        }
        x->Mat[i][0] = (b->Mat[i][0]- sum) / L->Mat[i][i];
    }

    return x;
}

Matrix* backwardSubstitution(Matrix* A, Matrix* b) {
    int N = b->Y;
    Matrix* x = new Matrix(N, 1);

    for (int i = N - 1; i >= 0; i--) {
        double sum = 0.0;
        for (int j = i + 1; j < N; j++) {
            sum += A->Mat[i][j] * x->Mat[j][0];
        }
        x->Mat[i][0] = (b->Mat[i][0] - sum) / A->Mat[i][i];
    }


    return x;
}


Matrix* Jacobi(Matrix* A,Matrix* b,double tol){
    int N=b->Y;
    Matrix* x = new Matrix(N,1);
    x->ones();

    Matrix* D  = A->Diag();
    Matrix* DInv = D->inv();

    Matrix* L  = A->LDiag();
    Matrix* U  = A->UDiag();

    Matrix* LU = L->add(U);

    double err=1;

    int iter =0;

    while(err > tol) {
        Matrix* LUx=LU->mul(x);
        Matrix* DInvLux=DInv->mul(LUx);
        Matrix* DinvB=DInv->mul(b);

        x->del();
        x=DinvB->sub(DInvLux);

        err = (A->mul(x)->sub(b))->norm();

        delete(LUx);
        delete(DInvLux);
        delete(DinvB);

        iter++;
        if(iter>1000) break;

    }

    delete(L);
    delete(D);
    delete(U);
    delete(LU);

    return x;

}

Matrix *LUFactorisationSolving(Matrix *A, Matrix *b) {
    Matrix** result;
    result = LUFactorisation(A);
//    printf("wynik faktoryzacji LU\n");
//    result[0]->print();
//    result[1]->print();


    Matrix * y = forwardSubstitution(result[0], b);

    Matrix * x = backwardSubstitution(result[1], y);

    delete(y);
    delete(result[0]);
    delete(result[1]);
    delete(result);
    return x;
}

void compare(Matrix* A,Matrix*B,Matrix*C,float tol){

    if(A->Y==B->Y && A->Y==C->Y && A->X==B->X && A->X==C->X){
        for (int i = 0; i < A->Y; ++i) {
            for (int j = 0; j < A->X; ++j) {
                if(!(abs(A->Mat[i][j] - B->Mat[i][j]) < tol
                && abs(A->Mat[i][j]-C->Mat[i][j]) < tol)){
                    printf("Matrixes have different values");

                    return;
                }
            }
        }
        printf("Matrixes have the same values");
        return;
    }
    printf("Matrixes have different sizes");
}

void createBandMatrix(Matrix* A) {
    int a1 = 5 + index[3];
    int a2 = -1;
    int a3 = a2;
    A->isBand= true;
    for (int i = 0; i < A->Y; ++i) {
        for (int j = 0; j < A->X; ++j) {
            if (i == j) {
                A->Mat[i][j] = a1;
            } else if (i == j + 1 || i + 1 == j) {
                A->Mat[i][j] = a2;
            } else if (i == j + 2 || i + 2 == j) {
                A->Mat[i][j] = a3;
//            } else {
//                Mat[i][j] = 0;
//            }
            }
        }
    }
}
