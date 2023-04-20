
#ifndef NUMERYCZNE2_MATRIX_H
#define NUMERYCZNE2_MATRIX_H


class Matrix {
public:
    int Y;
    int X;
    int N;
    bool isSquare;
    bool isBand;
    bool isDiagonal;

    double **Mat;
    Matrix(int y,int x);
    void print();
    void createBandMatrix(const int index[]);
    Matrix* inv();
    Matrix* add(Matrix* B);
    Matrix* sub(Matrix* B);
    Matrix* mul(Matrix*B);
    Matrix* neg();
    void ones();
    void del();
    Matrix* copy();
    double norm();
    void eq(Matrix* B);



    Matrix* D();
    Matrix* L();
    Matrix* U();
};


#endif //NUMERYCZNE2_MATRIX_H
