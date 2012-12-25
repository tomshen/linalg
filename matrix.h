#include <stdbool.h>

typedef struct matrix *Matrix;

typedef struct matrix {
    int r;
    int c;
    double** A;
} matrix;

bool is_matrix(Matrix M);
bool is_square_matrix(Matrix M);
Matrix matrix_new_empty(int row, int col);
Matrix matrix_new(double** array, int row, int col);
Matrix matrix_new_identity(int n);
Matrix matrix_copy(Matrix M);
void matrix_free(Matrix M);
void matrix_print(Matrix M);

bool is_normal_vector(double* v, int n);
double* vector_normalize(double* v, int n);
double vector_dot_product(double* v1, double* v2, int n);

bool matrix_equals(Matrix M1, Matrix M2);
Matrix matrix_add(Matrix M1, Matrix M2);
Matrix matrix_subtract(Matrix M1, Matrix M2);
Matrix matrix_multiply(Matrix M1, Matrix M2);
Matrix matrix_multiply_scalar(Matrix M, double n);

bool matrix_is_symmetric(Matrix M);
bool matrix_is_orthogonal(Matrix M);
bool matrix_is_identity(Matrix M);

Matrix matrix_transpose(Matrix M); // returns M^t
double* matrix_solve_system(Matrix A, double* b, int n); // returns solution
Matrix matrix_cofactor(Matrix M, int i, int j);
double matrix_determinant(Matrix M); // requires square matrix
Matrix matrix_inverse(Matrix M); // returns M^-1 or NULL if not possible