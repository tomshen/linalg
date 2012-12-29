// Matrix applications

#include <stdlib.h>
#include <assert.h>
#include <math.h>
#include "matrix.h"

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