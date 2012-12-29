// Classification
#include <stdlib.h>
#include <stdbool.h>
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