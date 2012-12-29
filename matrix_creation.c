// Creation
#include <stdlib.h>
#include <assert.h>
#include "matrix.h"

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