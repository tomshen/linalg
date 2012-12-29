// Basic operations
#include <stdlib.h>
#include <stdbool.h>
#include <assert.h>
#include <float.h>
#include <math.h>
#include "util.h"
#include "matrix.h"

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
        for(int j = 0; j < c; j++)
            for(int k = 0; k < m; k++)
                M->A[i][j] += M1->A[i][k] * M2->A[k][j];
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
Matrix matrix_cofactor(Matrix M, int i, int j)
{
    assert(is_matrix(M));
    assert(0 <= i && i < M->r);
    assert(0 <= j && j < M->c);
    Matrix C = matrix_new_empty(M->r - 1, M->c - 1);
    int r = 0;
    for(int k = 0; k < M->r; k++) {
        int c = 0;
        if(k != i) {
            for(int l = 0; l < M->c; l++) {
                if(l != j) {
                    C->A[r][c] = M->A[k][l];
                    c++;
                }
            }
            r++;
        }
    }
    assert(is_matrix(C));
    return C;
}
double matrix_determinant(Matrix M)
{
    assert(is_square_matrix(M));
    double temp;
    double det = 1;
    int n = M->c;

    if(n == 2) {
        return M->A[0][0] * M->A[1][1] - M->A[1][0] * M->A[0][1];
    }

    Matrix A = matrix_copy(M);
    for (int k = 0; k < n; k++) {
        int j = k;
        double max = A->A[j][j];
 
        for (int i = k + 1; i < n; i++)
            if ((temp = fabs(A->A[i][k])) > max) {
                j = i;
                max = temp;
            }
                
        matrix_swap_row(A, NULL, k, j);
        if(k != j)
            det *= -1;
 
        for (int i = k + 1; i < n; i++) {
            temp = A->A[i][k] / A->A[k][k];
            for (int j = k + 1; j < n; j++)
                A->A[i][j] -= temp * A->A[k][j];
            A->A[i][k] = 0;
        }
    }
    assert(is_square_matrix(A));
    for(int i = 0; i < n; i++)
        det *= A->A[i][i];
    free(A);
    return det;
}
Matrix matrix_inverse(Matrix M)
{
    if(!is_square_matrix(M))
        return NULL;
    int n = M->r;
    Matrix N = matrix_new_empty(n, n);
    double det = matrix_determinant(M);
    if(double_equals(det, 0))
        return NULL;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            Matrix C = matrix_cofactor(M, i, j);
            if((i + j) % 2 == 0)
                N->A[i][j] = matrix_determinant(C) / det;
            else
                N->A[i][j] = -1 * matrix_determinant(C) / det;
            matrix_free(C);
        }
    Matrix I = matrix_transpose(N);
    matrix_free(N);
    assert(is_square_matrix(I));
    return I;
}
Matrix* matrix_qr_decomposition(Matrix M)
{
    assert(is_matrix(M));
    int m = M->r;
    int n = M->c;
    Matrix Q = matrix_new_empty(m, n);
    for(int i = 0; i < n; i++) {
        double* u = matrix_column_vector(M, i);
        double* a = matrix_column_vector(M, i);
        for(int k = 0; k < i; k++) {
            double* e = matrix_column_vector(Q, k);
            double* p = vector_projection(e, a, m);
            free(e);
            for(int j = 0; j < m; j++)
                u[j] -= p[j];
            free(p);
        }
        free(a);
        double* e = vector_normalize(u, m);
        free(u);
        for(int l = 0; l < m; l++)
            Q->A[l][i] = e[l];
        free(e);
    }
    Matrix Qt = matrix_transpose(Q);
    Matrix R = matrix_multiply(Qt, M);
    matrix_free(Qt);

    Matrix* QR = calloc(2, sizeof(Matrix));
    QR[0] = Q;
    QR[1] = R;
    return QR;
}