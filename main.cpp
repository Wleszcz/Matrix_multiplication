#include <iostream>
#include "Matrix.h"
#include <cmath>
#include "matrixUtils.h"


int main() {

    const int len = sizeof(index)/ sizeof(int);
//    const int N = 100*9 +10*index[len-1]*index[len-2];
    const int N=1000;

    auto * A =new Matrix(N,N);
    A->createBandMatrix(index);

    auto * b = new Matrix(N,1);
    for (int i = 0; i < N; ++i) {
        b->Mat[i][0]=sin(i * (index[2] + 1));
    }


    Matrix * x = Jacobi(A,b,1e-9);
    printf("Jacobi done\n");
    //x->print();

    Matrix * x2 = GaussSeidel(A,b,1e-9);
    printf("Gauss done\n");
    //x2->print();


    Matrix* x3 =LUFactorisationSolving(A,b);
    printf("Lu done\n");
    //x3->print();
    compare(x,x2,x3,1e-6);

   delete(A);
   delete(b);
   delete(x);
   delete(x2);
   delete(x3);

    return 0;
}


