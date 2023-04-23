#include <iostream>
#include "Matrix.h"
#include <cmath>
#include "matrixUtils.h"
#include <vector>
#include <fstream>

#include <array>
#include <memory>

std::string exec(const char* cmd);

int main() {
    int N_size = 931;
    double tol = 1e-9;
    int a1 = 5 + index[3];
    int a2 = -1;
    int a3 = a2;

    //ZADANIE A - stworzenie układu równań
    auto * A =new Matrix(N_size, N_size);
    createBandMatrix(A,a1,a2,a3);

    auto * b = new Matrix(N_size, 1);
    for (int i = 0; i < N_size; ++i) {
        b->Mat[i][0]=sin(i * (index[2] + 1));
    }

    //ZADANIE B - rozwiązanie układów metodami iteracyjnymi i porównanie ich

    int time;
    Matrix * x1 = Jacobi(A,b,tol,&time);
    //x1->print();

    Matrix * x2 = GaussSeidel(A,b,tol,&time);
    //x2->print();

    compare(x1,x2,1e-8);

    //ZADANIE C - stworzenie nowego układu i próba rozwiązania iteracyjnie
    Matrix * Anew = new Matrix(N_size,N_size);
    createBandMatrix(Anew,3,-1,-1);

    //Matrix * x3 = Jacobi(Anew,b,1e-9);
    //Matrix * x4 = GaussSeidel(Anew,b,1e-9);


    //ZADANIE  D - Faktoryzacja LU
    Matrix* x5 =LUFactorisationSolving(Anew,b,&time);

    double err=((A->mul(x5))->sub(b))->norm();
    printf("Error: %.*e\n",err);


    delete(A);
    delete(b);
    delete(x1);
    delete(x2);
    delete(x5);

    //ZADANIE E - czasy

    int sizes[]={100,250,500,1000,1500, 2000,2500,3000};

    int len = sizeof(sizes)/sizeof(sizes[0]);

    int JacobiTimes[len];
    int GaussTime[len];
    int LUTimes[len];

    for (int i = 0; i < len; ++i) {
        Matrix * A = new Matrix(sizes[i],sizes[i]);

        createBandMatrix(A,a1,a2,a3);

        Matrix * b = new Matrix(sizes[i],1);

        for (int j = 0; j < sizes[i]; ++j) {
            b->Mat[j][0]=sin(j * (index[2] + 1));
        }
        A->isBand = true;

        Matrix * x1 = Jacobi(A,b,tol,&time);
        JacobiTimes[i]= time;

        Matrix * x2 = GaussSeidel(A,b,tol,&time);
        GaussTime[i]= time;

        Matrix * x3 = LUFactorisationSolving(A,b,&time);
        LUTimes[i]= time;

        delete(A);
        delete(b);
        delete(x1);
        delete(x2);
        delete(x3);
    }



    std::ofstream file ("time.txt");
    if (file.is_open())
    {
        for (int i = 0; i < 4; ++i) {
            for (int j = 0; j <len; ++j) {
                if(i==0) {
                    file << sizes[j]<<";";
                }
                if(i==1){
                    file << JacobiTimes[j]<<";";
                }
                if(i==2){
                    file << GaussTime[j]<<";";
                }
                if(i==3){
                    file << LUTimes[j]<<";";
                }

            }
            file<<std::endl;

        }
    }

    exec("conda activate SI");
    exec("python -u graph.py");

    return 0;
}

std::string exec(const char* cmd) {
    std::array<char, 128> buffer;
    std::string result;
    std::unique_ptr<FILE, decltype(&pclose)> pipe(popen(cmd, "r"), pclose);
    if (!pipe) {
        throw std::runtime_error("popen() failed!");
    }
    while (fgets(buffer.data(), buffer.size(), pipe.get()) != nullptr) {
        result += buffer.data();
    }
    return result;
}

