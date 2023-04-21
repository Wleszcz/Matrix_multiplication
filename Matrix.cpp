#include "Matrix.h"
#include "iostream"
#include "math.h"



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

void Matrix::print() {
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

Matrix* Matrix::D() {
    if(!isSquare){
        throw "Matrix is not square!";
    }
    Matrix *D = new Matrix(N,N);
    for (int i = 0; i < N; ++i) {
        D->Mat[i][i]=this->Mat[i][i];
    }
    D->isDiagonal=true;
    return D;
}

Matrix* Matrix::L() {
    if(!isBand){
        throw"not implemented!";
    }
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

Matrix* Matrix::U() {
    if(!isBand){
        throw"not implemented!";
    }
    Matrix *U = new Matrix(N,N);

    for (int i = 0; i < Y; ++i) {
        for (int j = 0; j < X; ++j) {
            if (i + 1 == j) {
                U->Mat[i][j] = Mat[i][j];
            }
            if (i + 2 == j){
                U->Mat[i][j] = Mat[i][j];
            }
        }
    }
    return U;
}

Matrix* Matrix::inv() {
    if(!isDiagonal){
        printf("not implemented");
        throw;
    }
    Matrix* res = new Matrix(Y,X);

    for (int i = 0; i < N; ++i) {
      res->Mat[i][i]=1/Mat[i][i];
    }
    return res;
}

Matrix* Matrix::add(Matrix *B) {

    if(this->X != B->X || this->Y != B->Y){
        printf("Wrong dimensions of matrix");
        return nullptr;
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
        printf("Wrong dimensions of matrix");
        return nullptr;
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
        printf("Wrong dimensions of matrix");
        return nullptr;
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





//dla macierzy wstęgowej bieżemy obecna iteracje, długość wstęgi

