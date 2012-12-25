// Some basic testing for matrix.c

#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "matrix.h"

void test1()
{
    double** A1 = calloc(2, sizeof(double*));
    A1[0] = calloc(2, sizeof(double));
    A1[1] = calloc(2, sizeof(double));
    A1[0][0] = 1;
    A1[0][1] = 2;
    A1[1][0] = 3;
    A1[1][1] = 4;
    Matrix M1 = matrix_new(A1, 2, 2);
    Matrix M2 = matrix_copy(M1);
    assert(matrix_equals(M1, M2));
    matrix_print(M1);
    printf("\n");
    matrix_print(M2);
    Matrix M3 = matrix_multiply(M1, M2);
    printf("\n");
    matrix_print(M3);
    Matrix M4 = matrix_multiply_scalar(M3, 10);
    printf("\n");
    matrix_print(M4);
    matrix_free(M1);
    matrix_free(M2);
    matrix_free(M3);
    matrix_free(M4);
}

void test2()
{
    double** A1 = calloc(2, sizeof(double*));
    A1[0] = calloc(2, sizeof(double));
    A1[1] = calloc(2, sizeof(double));
    A1[0][0] = 1;
    A1[0][1] = 2;
    A1[1][0] = 3;
    A1[1][1] = 4;
    Matrix M1 = matrix_new(A1, 2, 2);
    Matrix M2 = matrix_copy(M1);
    Matrix M4 = matrix_copy(M1);
    M2->A[1][0] = M2->A[0][1];
    M2->A[1][1] = M2->A[0][0];
    assert(!is_symmetric_matrix(M1));
    assert(is_symmetric_matrix(M2));

    Matrix M3 = matrix_copy(M2);
    M3->A[0][1] = 0;
    M3->A[1][0] = 0;
    matrix_print(M3);
    assert(!is_orthogonal_matrix(M2));
    assert(is_orthogonal_matrix(M3));

    printf("\n");
    matrix_print(M1);
    double* b = calloc(2, sizeof(double));
    b[0] = 1;
    b[1] = 1;
    double* x = solve_system(M1, b, 2);
    printf("%f\t%f\n", x[0], x[1]);
    printf("\n");
    matrix_print(M1);

    matrix_print(M4);
    printf("%f\n", matrix_determinant(M4));
}

void test3()
{
    double** A1 = calloc(3, sizeof(double*));
    A1[0] = calloc(3, sizeof(double));
    A1[1] = calloc(3, sizeof(double));
    A1[2] = calloc(3, sizeof(double));
    A1[0][0] = 1;
    A1[0][1] = 2;
    A1[0][2] = 0;
    A1[1][0] = -1;
    A1[1][1] = 1;
    A1[1][2] = 1;
    A1[2][0] = 1;
    A1[2][1] = 2;
    A1[2][2] = 3;
    Matrix M1 = matrix_new(A1, 3, 3);
    matrix_print(M1);
    printf("%f\n", matrix_determinant(M1));
    Matrix M2 = matrix_inverse(M1);
    matrix_print(M2);

    Matrix I = matrix_new_identity(3);
    matrix_print(I);
    assert(!(is_identity_matrix(M1)));
    assert(is_identity_matrix(I));
}

void test4()
{
    double** A1 = calloc(3, sizeof(double*));
    A1[0] = calloc(3, sizeof(double));
    A1[1] = calloc(3, sizeof(double));
    A1[2] = calloc(3, sizeof(double));
    A1[0][0] = 12;
    A1[0][1] = -51;
    A1[0][2] = 4;
    A1[1][0] = 6;
    A1[1][1] = 167;
    A1[1][2] = -68;
    A1[2][0] = -4;
    A1[2][1] = 24;
    A1[2][2] = -41;
    Matrix M1 = matrix_new(A1, 3, 3);
    matrix_print(M1);
    Matrix* QR = matrix_qr_decomposition(M1);
    matrix_print(QR[0]);
    matrix_print(QR[1]);
}

void test5()
{
    double** A1 = calloc(3, sizeof(double*));
    A1[0] = calloc(2, sizeof(double));
    A1[1] = calloc(2, sizeof(double));
    A1[2] = calloc(2, sizeof(double));
    A1[0][0] = 3;
    A1[0][1] = -6;
    A1[1][0] = 4;
    A1[1][1] = -8;
    A1[2][0] = 0;
    A1[2][1] = 1;
    Matrix A = matrix_new(A1, 3, 2);
    matrix_print(A);
    Matrix* QR = matrix_qr_decomposition(A);
    matrix_print(QR[0]);
    matrix_print(QR[1]);
    
    double* b = calloc(3, sizeof(double));
    b[0] = -1;
    b[1] = 7;
    b[2] = 2;
    double* x = least_squares_regression(A, b);
    for(int i = 0; i < 2; i++)
        printf("%g\n", x[i]);
}
int main()
{
    test1();
    test2();
    test3();
    test4();
    test5();
    return 0;
}