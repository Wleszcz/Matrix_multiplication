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
    x->print();

    Matrix * x2 = GaussSeidel(A,b,1e-14);

    x2->print();


   delete(A);
   delete(b);
   delete(x);
   delete(x2);

    return 0;
}


