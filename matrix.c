#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include "matrix.h"

double* array_copy(double* A, int n)
{
    double* B = calloc(n, sizeof(double));
    for(int i = 0; i < n; i++)
        B[i] = A[i];
    return B;
}

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

bool is_normal_vector(double* v, int n)
{
    int sum = 0;
    for(int i = 0; i < n; i++)
        sum += v[i] * v[i];
    return sum == 1;
}
double* vector_normalize(double* v, int n)
{
    double* vn = calloc(n, sizeof(double));
    double sum = 0;
    for(int i = 0; i < n; i++)
        sum += v[i] * v[i];
    double norm = sqrt(sum);
    for(int i = 0; i < n; i++)
        vn[i] = v[i] / norm;
    return vn;
}
double vector_dot_product(double* v1, double* v2, int n)
{
    double product = 0;
    for(int i = 0; i < n; i++)
        product += v1[i] * v2[i];
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

bool matrix_is_symmetric(Matrix M)
{
    if(!is_square_matrix(M))
        return false;
    int n = M->r;
    for(int i = 0; i < n; i++)
        for(int j = i; j < n; j++)
            if(M->A[i][j] != M->A[j][i])
                return false;
    return true;
}
bool matrix_is_orthogonal(Matrix M)
{
    if(!is_square_matrix(M))
        return false;
    int n = M->c;
    for(int i = 0; i < n - 1; i++) {
        double* v1 = array_copy(M->A[i], n);
        if(!is_normal_vector(v1, n)) {
            free(v1);
            return false;
        }// fabs(vector_dot_product(v1, v2, n) - 1) > DBL_EPSILON
        for(int j = i + 1; j < n; j++) {
            double* v2 = array_copy(M->A[j], n);
            if(!is_normal_vector(v2, n) || 
                vector_dot_product(v1, v2, n) > DBL_EPSILON) {
                free(v1);
                free(v2);
                return false;
            }
            free(v2);
        }
        free(v1);
    }
    return true;
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
double* matrix_solve_system(Matrix A, double* b, int n)
{
    assert(is_square_matrix(A));
    assert(n == A->c);
    double* x = calloc(n, sizeof(double));
    double temp;
    
    for (int k = 0; k < n; k++) {
        int j = k;
        double max = A->A[j][j];
 
        for (int i = k + 1; i < n; i++)
            if ((temp = fabs(A->A[i][k])) > max) {
                j = i;
                max = temp;
            }
                
        matrix_swap_row(A, b, k, j);
 
        for (int i = k + 1; i < n; i++) {
            temp = A->A[i][k] / A->A[k][k];
            for (int j = k + 1; j < n; j++)
                A->A[i][j] -= temp * A->A[k][j];
            A->A[i][k] = 0;
            b[i] -= temp * b[k];
        }
    }
    for (int k = n - 1; k >= 0; k--) {
        temp = b[k];
        for (int j = n - 1; j > k; j--)
            temp -= x[j] * A->A[k][j];
        x[k] = temp / A->A[k][k];
    }
    assert(is_square_matrix(A));
    return x;
}
Matrix matrix_inverse(Matrix M)
{
    if(!is_square_matrix(M))
        return NULL;
    int n = M->r;
    Matrix N = matrix_new_empty(n, n);
    // TODO: fill in actual code here
    return N;
}
double matrix_determinant(Matrix M);