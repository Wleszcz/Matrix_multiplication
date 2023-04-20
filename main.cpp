#include <iostream>
#include "Matrix.h"
#include "math.h"


Matrix *Jacobi(Matrix* A,Matrix* b,double tol);

Matrix *GaussSeidel(Matrix *A, Matrix *p, double tol);

Matrix *forwardSubstitution(Matrix *L, Matrix *b);

int main() {
    const int index[] ={1,8,8,7,3,1};
    const int len = sizeof(index)/ sizeof(int);
//    const int N = 100*9 +10*index[len-1]*index[len-2];
    const int N=10;

    Matrix * A =new Matrix(N,N);

    A->createBandMatrix(index);
    A->print();

    Matrix * b = new Matrix(N,1);

    for (int i = 0; i < N; ++i) {
        b->Mat[i][0]=sin(i * (index[2] + 1));
    }
    Matrix * x = Jacobi(A,b,1e-14);
    x->print();


    //-0.017428
    A->print();
    Matrix * x2 = GaussSeidel(A,b,1e-14);


    x2->print();

    A->del();
    b->del();

    return 0;
}

Matrix *GaussSeidel(Matrix *A, Matrix *b, double tol) {
    int N=b->Y;
    Matrix* x = new Matrix(N,1);
    Matrix* x_pr;
    x->ones();

    Matrix* D  = A->D();
    Matrix* DInv = D->inv();

    Matrix* L  = A->L();
    Matrix* U  = A->U();
    Matrix* DL = D->add(L);
    double err=1;

    int iter =0;
    while(err>tol){

        x_pr=x->copy();

        // forward substitution
        Matrix *ux=U->mul(x_pr);

        Matrix * b_ux=b->sub(ux);

        Matrix* z = forwardSubstitution(DL,b_ux);

        //update solution
        x = forwardSubstitution(DL,z->neg());

        //x->print();

        err=x_pr->sub(x)->norm();
        iter++;
        if(iter>1000) break;
    }
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
        x->Mat[i][0]=  (b->Mat[i][0]- sum) / L->Mat[i][i];
    }

    return x;
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
