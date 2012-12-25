#include <stdlib.h>
#include <stdbool.h>
#include <float.h>

double* array_copy(double* A, int n)
{
    double* B = calloc(n, sizeof(double));
    for(int i = 0; i < n; i++)
        B[i] = A[i];
    return B;
}
bool double_equals(double d, int i)
{
    return d - i <= DBL_EPSILON;
}