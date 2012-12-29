#include <stdlib.h>
#include <stdio.h>
#include <stdbool.h>
#include <math.h>
#include <float.h>
#include <assert.h>
#include <limits.h>
#include "util.h"
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
bool is_symmetric_matrix(Matrix M)
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
bool is_orthogonal_matrix(Matrix M)
{
    if(!is_square_matrix(M))
        return false;
    int n = M->c;
    for(int i = 0; i < n - 1; i++) {
        double* v1 = array_copy(M->A[i], n);
        if(!is_normal_vector(v1, n)) {
            free(v1);
            return false;
        }
        for(int j = i + 1; j < n; j++) {
            double* v2 = array_copy(M->A[j], n);
            if(!(is_normal_vector(v2, n) && 
                double_equals(vector_dot_product(v1, v2, n), 0))) {
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
bool is_identity_matrix(Matrix M)
{
    if(!is_square_matrix(M)) {
        return false;
    }
    int n = M->r;
    for(int i = 0; i < n; i++)
        for(int j = 0; j < n; j++) {
            if(i == j && !double_equals(M->A[i][j], 1))
                return false;
            else if(i != j && !double_equals(M->A[i][j], 0))
                return false;
        }
    return true;
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
Matrix matrix_new_empty(int row, int col)
{
    assert(row > 0 && col > 0);
    Matrix M = malloc(sizeof(matrix));
    M->r = row;
    M->c = col;
    double** A = calloc(row, sizeof(double*));
    for(int i = 0; i < row; i++)
        A[i] = calloc(col, sizeof(double));
    M->A = A;
    assert(is_matrix(M));
    return M;
}
Matrix matrix_new_identity(int n)
{
    assert(n > 0);
    Matrix M = matrix_new_empty(n, n);
    for(int i = 0; i < n; i++)
        M->A[i][i] = 1;
    assert(is_square_matrix(M));
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
double* vector_projection(double* e, double* a, int n)
{
    double ea = vector_dot_product(e, a, n);
    double ee = vector_dot_product(e, e, n);
    double* p = calloc(n, sizeof(double));
    for(int i = 0; i < n; i++)
        p[i] = ea * e[i] / ee;
    return p;
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
double* matrix_column_vector(Matrix M, int c)
{
    assert(is_matrix(M));
    double* v = calloc(M->r, sizeof(double));
    for(int i = 0; i < M->r; i++)
        v[i] = M->A[i][c];
    return v;
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


double* solve_system(Matrix A, double* b, int n)
{
    assert(is_square_matrix(A));
    assert(n == A->r);
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
double* least_squares_regression(Matrix A, double* b)
{
    Matrix* QR = matrix_qr_decomposition(A);
    Matrix Qt = matrix_transpose(QR[0]);
    double* d = calloc(Qt->r, sizeof(double));
    for(int i = 0; i < Qt->r; i++)
        for(int j = 0; j < Qt->c; j++)
            d[i] += Qt->A[i][j] * b[j];
    free(Qt);
    matrix_print(QR[1]);
    double* x = solve_system(QR[1], d, Qt->r);
    matrix_free(QR[1]);
    free(d);
    free(QR);
    return x;
}