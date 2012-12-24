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
Matrix matrix_copy(Matrix M);
void matrix_free(Matrix M);
void matrix_print(Matrix M);

bool is_vector(Matrix V);
bool is_normal_vector(Matrix V);
Matrix vector_normalize(Matrix V);
double vector_dot_product(Matrix V1, Matrix V2);

bool matrix_equals(Matrix M1, Matrix M2);
Matrix matrix_add(Matrix M1, Matrix M2);
Matrix matrix_subtract(Matrix M1, Matrix M2);
Matrix matrix_multiply(Matrix M1, Matrix M2);
Matrix matrix_multiply_scalar(Matrix M, double n);

bool matrix_is_symmetric(Matrix M);
bool matrix_is_orthogonal(Matrix M);

Matrix matrix_transpose(Matrix M); // returns M^t
Matrix matrix_lu_decompose(Matrix M); // returns L
Matrix matrix_inverse(Matrix M); // returns M^-1 or NULL if not possible
double matrix_determinant(Matrix M); // requires square matrix