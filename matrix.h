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
Matrix vector_normalize(Matrix V);
double vector_dot_product(Matrix V1, Matrix V2);

bool matrix_equals(Matrix M1, Matrix M2);
Matrix matrix_add(Matrix M1, Matrix M2);
Matrix matrix_subtract(Matrix M1, Matrix M2);
Matrix matrix_multiply(Matrix M1, Matrix M2);
Matrix matrix_multiply_scalar(Matrix M, double n);
Matrix matrix_transpose(Matrix M);
Matrix matrix_inverse(Matrix M);