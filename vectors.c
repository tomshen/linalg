// Vectors
#include <stdlib.h>
#include <math.h>
#include "matrix.h"

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