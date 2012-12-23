#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "matrix.h"

bool is_matrix(Matrix M)
{
    if(M == NULL)
        return false;
    else if(M->r <= 0)
        return false;
    else if(M->c <= 0)
        return false;
    else if(M->A == NULL)
        return false;
    else
        return true;
}
bool is_square_matrix(Matrix M)
{
    return is_matrix(M) && M->r == M->c;
}
Matrix matrix_new_empty(int row, int col)
{
    assert(row > 0 && col > 0);
    Matrix M = malloc(sizeof(matrix));
    M->r = row;
    M->c = col;
    double** A = calloc(row, sizeof(double*));
    for(int i = 0; i < row; i++) {
        A[i] = calloc(col, sizeof(double));
    }
    M->A = A;
    assert(is_matrix(M));
    return M;
}
Matrix matrix_new(double** array, int row, int col)
{
    assert(row > 0 && col > 0);
    assert(array != NULL);
    Matrix M = malloc(sizeof(matrix));
    M->r = row;
    M->c = col;
    M->A = array;
    assert(is_matrix(M));
    return M;
}
Matrix matrix_copy(Matrix M)
{
    assert(is_matrix(M));
    int r = M->r;
    int c = M->c;
    double** A = calloc(r, sizeof(double*));
    for(int i = 0; i < r; i++) {
        A[i] = calloc(c, sizeof(double));
        for(int j = 0; j < c; j++)
            A[i][j] = M->A[i][j];
    }
    Matrix C = matrix_new(A, r, c);
    assert(is_matrix(C));
    return C;
}
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
    assert(is_matrix(M));
    for(int i = 0; i < M->r; i++) {
        for(int j = 0; j < M->c; j++)
            printf("%.3f\t", M->A[i][j]);
        printf("\n");
    }
}

bool is_vector(Matrix V)
{
    assert(is_matrix(V));
    return V->r == 1 || V->c == 1;
}
Matrix vector_normalize(Matrix V)
{
    assert(is_vector(V));
    double sum = 0;
    if(V->r == 1) {
        for(int i = 0; i < V->c; i++)
            sum += V->A[0][i] * V->A[0][i];
    }
    else {
        for(int i = 0; i < V->r; i++)
            sum += V->A[i][0] * V->A[i][0];
    }
    Matrix N = matrix_multiply_scalar(V, sqrt(1 / sum));
    assert(is_vector(N));
    return N;
}
double vector_dot_product(Matrix V1, Matrix V2)
{
    assert(is_vector(V1) && is_vector(V2));
    double product = 0;
    if(V1->r == 1) {
        assert(V2->c == 1 && V1->c == V2->r);
        for(int i = 0; i < V1->c; i++)
            product += V1->A[0][i] * V2->A[i][0];
    }
    else {
        assert(V2->r == 1 && V1->r == V2->c);
        for(int i = 0; i < V1->r; i++)
            product += V1->A[i][0] * V2->A[0][i];
    }
    return product;
}

bool matrix_equals(Matrix M1, Matrix M2)
{
    assert(is_matrix(M1) && is_matrix(M2));
    if(!(M1->r == M2->r && M1->c == M2->c))
        return false;
    for(int i = 0; i < M1->r; i++)
        for(int j = 0; j < M1->c; j++)
            if(M1->A[i][j] != M2->A[i][j])
                return false;
    return true;
}
Matrix matrix_add(Matrix M1, Matrix M2)
{
    assert(is_matrix(M1) && is_matrix(M2));
    assert(M1->r == M2->r && M1->c == M2->c);
    int r = M1->r;
    int c = M1->c;
    Matrix M = matrix_new_empty(r, c);
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++) {
            assert(FLT_MAX - M2->A[i][j] >= M1->A[i][j]);
            M->A[i][j] = M1->A[i][j] + M2->A[i][j];
        }
    assert(is_matrix(M));
    return M;
}
Matrix matrix_subtract(Matrix M1, Matrix M2)
{
    assert(is_matrix(M1) && is_matrix(M2));
    assert(M1->r == M2->r && M1->c == M2->c);
    int r = M1->r;
    int c = M1->c;
    Matrix M = matrix_new_empty(r, c);
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++) {
            assert(FLT_MIN + M1->A[i][j] >= M2->A[i][j]);
            M->A[i][j] = M1->A[i][j] - M2->A[i][j];
        }
    assert(is_matrix(M));
    return M;
}
Matrix matrix_multiply(Matrix M1, Matrix M2)
{
    assert(is_matrix(M1) && is_matrix(M2));
    assert(M1->c == M2->r);
    int r = M1->r;
    int c = M2->c;
    int m = M1->c;
    Matrix M = matrix_new_empty(r, c);
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++) {
            for(int k = 0; k < m; k++) {
                assert((FLT_MAX - M->A[i][j]) / M1->A[i][k] >= M2->A[k][j]);
                M->A[i][j] += M1->A[i][k] * M2->A[k][j];
            }
        }
    assert(is_matrix(M));
    return M;
}
Matrix matrix_multiply_scalar(Matrix M, double n)
{
    assert(is_matrix(M));
    int r = M->r;
    int c = M->c;
    Matrix N = matrix_new_empty(r, c);
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++) {
                assert(FLT_MAX / n >= M->A[i][j]);
                N->A[i][j] = n * M->A[i][j];
            }
    assert(is_matrix(N));
    return N;
}
Matrix matrix_transpose(Matrix M)
{
    assert(is_matrix(M));
    int r = M->r;
    int c = M->c;
    Matrix N = matrix_new_empty(c, r);
    for(int i = 0; i < r; i++)
        for(int j = 0; j < c; j++)
            N->A[j][i] = M->A[i][j];
    assert(is_matrix(N));
    return N;
}
Matrix matrix_inverse(Matrix M)
{
    assert(is_square_matrix(M));
    int n = M->r;
    Matrix N = matrix_new_empty(n, n);
    // TODO: fill in actual code here
    assert(matrix_multiply(M, N));
    return N;
}