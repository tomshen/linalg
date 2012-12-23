#include <stdio.h>
#include <stdlib.h>
#include <assert.h>
#include "matrix.h"

int main()
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

    double** A3 = calloc(3, sizeof(double*));
    A3[0] = calloc(1, sizeof(double));
    A3[1] = calloc(1, sizeof(double));
    A3[2] = calloc(1, sizeof(double));
    A3[0][0] = 1;
    A3[1][0] = 2;
    A3[2][0] = 3;

    double** A4 = calloc(1, sizeof(double*));
    A4[0] = calloc(3, sizeof(double));
    A4[0][0] = 1;
    A4[0][1] = 2;
    A4[0][2] = 3;

    Matrix V1 = matrix_new(A3, 3, 1);
    Matrix V2 = matrix_new(A4, 1, 3);
    Matrix V3 = vector_normalize(V1);

    printf("\n");
    printf("%.4f\n", vector_dot_product(V1, V2));
    printf("\n");
    matrix_print(V3);
    Matrix V4 = matrix_transpose(V3);
    printf("\n");
    matrix_print(V4);
    matrix_free(V1);
    matrix_free(V2);
    matrix_free(V3);
    matrix_free(V4);
    return 0;
}