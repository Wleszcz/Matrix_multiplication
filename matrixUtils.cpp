//
// Created by Wiktor on 21.04.2023.
//
#include "matrixUtils.h"
#include "Matrix.h"

Matrix *GaussSeidel(Matrix *A, Matrix *b, double tol) {
    int N=b->Y;
    Matrix* x = new Matrix(N,1);
    Matrix* x_pr = new Matrix(N,1);
    x->ones();

    Matrix* D  = A->D();
    Matrix* L  = A->L();
    Matrix* U  = A->U();
    Matrix* DL = D->add(L);
    double err=1;

    int iter =0;
    while(err>tol){
        x_pr->eq(x);

        //solving  A = (D+L)^-1 (b-U*x_pr) using forward substitution

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

    Matrix* D  = A->D();
    Matrix* DInv = D->inv();

    Matrix* L  = A->L();
    Matrix* U  = A->U();

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
