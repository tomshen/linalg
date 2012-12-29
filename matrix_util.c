// Utils
#include <stdlib.h>
#include <stdio.h>
#include <limits.h>
#include <assert.h>
#include "matrix.h"

void matrix_free(Matrix M)
{
    assert(is_matrix(M));
    for(int i = 0; i < M->r; i++)
        free(M->A[i]);
    free(M->A);
    free(M);
}
void matrix_print(Matrix M)
{
    matrix_print_format(M, 3);
}
void matrix_print_format(Matrix M, int n)
{
    assert(is_matrix(M));
    for(int k = 0; k < M->c; k++)
        printf("--------");
    printf("\n");

    // ensures you don't print something like -0.000
    int c = 1;
    for(int i = 0; i < n; i++) {
        assert(INT_MAX / 10 >= c);
        c *= 10;
    }

    for(int i = 0; i < M->r; i++) {
        for(int j = 0; j < M->c; j++) {
            if((int)(M->A[i][j] * c) == (int)(M->A[i][j]) * c)
                printf(" %-d\t", (int)M->A[i][j]);
            else
                printf(" %-.*f\t", n, M->A[i][j]);  
        }
        printf("\n");
    }
    for(int k = 0; k < M->c; k++)
        printf("--------");
    printf("\n");
}
void matrix_swap_row(Matrix M, double* b, int r1, int r2)
{
    assert(is_matrix(M));
    if (r1 == r2)
        return;
    int c = M->c;
    for (int i = 0; i < c; i++) {
        double temp = M->A[r1][i];
        M->A[r1][i] = M->A[r2][i];
        M->A[r2][i] = temp;
    }
    if(b != NULL) {
        double temp = b[r1];
        b[r1] = b[r2];
        b[r2] = temp;
    }
}
double* matrix_column_vector(Matrix M, int c)
{
    assert(is_matrix(M));
    double* v = calloc(M->r, sizeof(double));
    for(int i = 0; i < M->r; i++)
        v[i] = M->A[i][c];
    return v;
}