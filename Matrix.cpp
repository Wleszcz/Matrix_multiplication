#include "Matrix.h"
#include "iostream"
#include <cmath>



Matrix::Matrix(int y,int x) {
    this->X = x;
    this->Y = y;
    this->N= -1;
    Mat = new double *[Y];

    for (int i = 0; i < Y; ++i) {
        Mat[i]= new double[x];
        for (int j = 0; j < x; ++j) {
            Mat[i][j]=0;
        }
    }
    isSquare=false;
    isBand= false;
    isDiagonal=false;

    if(X==Y){
        isSquare=true;
        N=Y;
    }
}

Matrix::~Matrix() {
    del();
}

void Matrix::print() const {
    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            printf("%f ", Mat[i][j]);
        }
        printf("\n");
    }
    printf("\n\n");
}

void Matrix::createBandMatrix(const int *index) {
    int a1 = 5 + index[3];
    int a2 = -1;
    int a3 = a2;
    isBand= true;
    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            if (i == j) {
                Mat[i][j] = a1;
            } else if (i == j + 1 || i + 1 == j) {
                Mat[i][j] = a2;
            } else if (i == j + 2 || i + 2 == j) {
                Mat[i][j] = a3;
//            } else {
//                Mat[i][j] = 0;
//            }
            }
        }
    }
}

Matrix* Matrix::Diag() {
    if(isSquare) {
        auto *D = new Matrix(N, N);
        for (int i = 0; i < N; ++i) {
            D->Mat[i][i] = this->Mat[i][i];
        }
        D->isDiagonal = true;
        return D;
    }
    throw std::invalid_argument("Diag func. for this type of matrix not implemented!");
}

Matrix* Matrix::LDiag() {
    if(isBand){
        Matrix *L = new Matrix(N,N);

        for (int i = 0; i < Y; ++i) {
            for (int j = 0; j < X; ++j) {
                if (i == j + 1 ) {
                    L->Mat[i][j] = Mat[i][j];
                }
                if (i == j + 2 ){
                    L->Mat[i][j] = Mat[i][j];
                }
            }
        }
        return L;
    }

    throw std::invalid_argument("LDiag func. for this type of matrix not implemented!");
}

Matrix* Matrix::UDiag() {
    if(isBand) {
        Matrix *U = new Matrix(N, N);

        for (int i = 0; i < Y; ++i) {
            for (int j = 0; j < X; ++j) {
                if (i + 1 == j) {
                    U->Mat[i][j] = Mat[i][j];
                }
                if (i + 2 == j) {
                    U->Mat[i][j] = Mat[i][j];
                }
            }
        }
        return U;
    }

    throw std::invalid_argument("UDiag func. for this type of matrix not implemented!");
}

Matrix* Matrix::inv() {
    if(isDiagonal){
        Matrix* res = new Matrix(Y,X);

        for (int i = 0; i < N; ++i) {
            res->Mat[i][i]=1/Mat[i][i];
        }
        return res;
    }
    throw std::invalid_argument("inverting this type of matrix not implemented");
}

Matrix* Matrix::add(Matrix *B) {

    if(this->X != B->X || this->Y != B->Y){
        throw std::invalid_argument("Wrong dimensions of matrix (add)");
    }

    Matrix* res = new Matrix(Y,X);

    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            res->Mat[i][j]= Mat[i][j]+B->Mat[i][j];
        }
    }
    return res;
}

Matrix * Matrix::sub(Matrix *B) {
    if(this->X != B->X || this->Y != B->Y){
        throw std::invalid_argument("Wrong dimensions of matrix (sub)");
    }
    Matrix* res = new Matrix(Y,X);

    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            res->Mat[i][j]= Mat[i][j]-B->Mat[i][j];
        }
    }
    return res;
}

Matrix *Matrix::mul(Matrix *B) {
    if(this->X!=B->Y){
        throw std::invalid_argument("Wrong dimensions of matrix (mull)");
    }
    Matrix* result = new Matrix(this->Y,B->X);
    for(int i = 0; i < this->Y; i++ ) {
        for (int j = 0; j < B->X; j++) {
            double s = 0;
            for (int k = 0; k < B->Y; k++){
                s += this->Mat[i][k] * B->Mat[k][j];
            }
            result->Mat[i][j] = s;
        }
    }
    return result;
}

void Matrix::ones() {
    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            Mat[i][j]=1;
        }
    }
}

void Matrix::del() {
    for (int i = 0; i < Y; ++i) {
        delete[] Mat[i];
    }
    delete[] Mat;
}

Matrix *Matrix::copy() {
    Matrix* nm = new Matrix(Y,X);
    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            nm->Mat[i][j]=Mat[i][j];
        }
    }
    return nm;
}

double Matrix::norm() {
    double sum = 0;
    for (int i = 0; i < Y; i++) {
        for (int j = 0; j < X; j++) {
            sum += Mat[i][j] * Mat[i][j];
        }
    }
    return sqrt(sum);
}

void Matrix::eq(Matrix *B) {
    if(X!=B->X || Y!=B->Y){
        throw std::invalid_argument("Wrong dimensions of matrix (equal)");
    }
    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            Mat[i][j]=B->Mat[i][j];
        }
    }
}

Matrix *Matrix::neg() {
    Matrix* nm = new Matrix(Y,X);
    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            nm->Mat[i][j]=-Mat[i][j];
        }
    }
    return nm;
}


