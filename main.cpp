#include <iostream>
#include "Matrix.h"
#include "math.h"


Matrix *Jacobi(Matrix* A,Matrix* b,double tol);

int main() {
    const int index[] ={1,8,8,7,3,1};
    const int len = sizeof(index)/ sizeof(int);
//    const int N = 100*9 +10*index[len-1]*index[len-2];
    const int N=5000;

    Matrix * A =new Matrix(N,N);

    A->createBandMatrix(index);

    Matrix * b = new Matrix(N,1);

    for (int i = 0; i < N; ++i) {
        b->Mat[i][0]=sin(i * (index[2] + 1));
    }
    Matrix * x = Jacobi(A,b,1e-14);

    x->print();

    A->del();
    b->del();

    return 0;
}

Matrix* Jacobi(Matrix* A,Matrix* b,double tol){
    int N=b->Y;
    Matrix* x = new Matrix(N,1);
    Matrix* x_pr = NULL;
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

        x=DinvB->sub(DInvLux);

        err = (A->mul(x)->sub(b))->norm();

        LUx->del();
        delete(LUx);

        DInvLux->del();
        delete(DInvLux);

        DinvB->del();
        delete(DinvB);


        iter++;
        if(iter>1000) break;

    }

    L->del();
    delete(L);

    D->del();
    delete(D);

    U->del();
    delete(U);

    return x;

}
