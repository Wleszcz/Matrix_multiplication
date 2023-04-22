#include <iostream>
#include "Matrix.h"
#include <cmath>
#include "matrixUtils.h"


int main() {
    const int index[] ={1,8,8,7,3,1};
    const int len = sizeof(index)/ sizeof(int);
//    const int N = 100*9 +10*index[len-1]*index[len-2];
    const int N=10;

    auto * A =new Matrix(N,N);
    A->createBandMatrix(index);

    auto * b = new Matrix(N,1);
    for (int i = 0; i < N; ++i) {
        b->Mat[i][0]=sin(i * (index[2] + 1));
    }


    Matrix * x = Jacobi(A,b,1e-14);
    //x->print();

    Matrix * x2 = GaussSeidel(A,b,1e-14);

    //x2->print();

    Matrix * he = new Matrix(3,3);
    he->Mat[0][0]=1;
    he->Mat[0][1]=-3;
    he->Mat[0][2]=0;

    he->Mat[1][0]=0;
    he->Mat[1][1]=1;
    he->Mat[1][2]=3;

    he->Mat[2][0]=2;
    he->Mat[2][1]=-10;
    he->Mat[2][2]=2;


    Matrix* x3 =LUFactorisationSolving(he,b);

   delete(A);
   delete(b);
   delete(x);
   delete(x2);

    return 0;
}


