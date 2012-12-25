#include <stdbool.h>

/* A Matrix is a pointer to a struct. The struct contains a 2D-array that 
 * contains the actual entries, and ints representing the number of rows and 
 * columns. This is to avoid passing around those values separately, which 
 * would create quite convoluted function headers. */
typedef struct matrix *Matrix;

typedef struct matrix {
    int r;
    int c;
    double** A;
} matrix;


// Classification - to check if a matrix is of a specifc type
bool is_matrix(Matrix M);
bool is_square_matrix(Matrix M);
bool is_symmetric_matrix(Matrix M);
bool is_orthogonal_matrix(Matrix M);
bool is_identity_matrix(Matrix M);


// Creation - different ways to create a matrix
Matrix matrix_new(double** array, int row, int col);
Matrix matrix_new_empty(int row, int col);
Matrix matrix_new_identity(int n);
Matrix matrix_copy(Matrix M);


// Meta operations
void matrix_free(Matrix M);
void matrix_print(Matrix M);


// Vectors - common vector operations
bool is_normal_vector(double* v, int n);
double* vector_normalize(double* v, int n);
double vector_dot_product(double* v1, double* v2, int n);
double* vector_projection(double* v1, double* v2, int n);


// Arithmetic - basic operations with matrices
bool matrix_equals(Matrix M1, Matrix M2);
Matrix matrix_add(Matrix M1, Matrix M2);
Matrix matrix_subtract(Matrix M1, Matrix M2);
Matrix matrix_multiply(Matrix M1, Matrix M2);
Matrix matrix_multiply_scalar(Matrix M, double n);


// Useful operations - more complicated, common operations with matrices

// returns the transpose of M
Matrix matrix_transpose(Matrix M);

// returns the (i, j) cofactor of M
Matrix matrix_cofactor(Matrix M, int i, int j);

// requires square matrix
/* finds determinant by triangulizing M through gaussian elimination, taking 
 * the product of the diagonal entries, and then multiplying by -1^(number of
 * permutations) */
double matrix_determinant(Matrix M); 

/* finds inverse by finding the adjugate matrix of M (tranpose of matrix of 
 * cofactors), and then dividing by the determinant of M. If det(M) = 0, or if
 * the matrix is not square, then no inverse exists, and NULL is returned. */
Matrix matrix_inverse(Matrix M);

/* finds the QR decomposition of M through the Gram-Schmidt process, which 
 * relies on projection and normalization of the column vectors. returns
 * an array of matrices, with Q as the first and R as the second entry. */
Matrix* matrix_qr_decomposition(Matrix M);


// Applications - common applications of matrices in linear algebra

/* Given matrix A and vector b (and the length of b, int n), this function 
 * finds the solution vector x to the system Ax=b. */
double* solve_system(Matrix A, double* b, int n);

/* Given matrix A and vector b (and the length of b, int n), this function 
 * finds the most accurate coefficients for the overdetermined system Ax=b.
 * This function uses QR decomposition to accomplish this. */
double* least_squares_regression(Matrix A, double* b);