#pragma once

/** Only standardlibraries */
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <string.h>

/** This is the core-struct in this library. All matrix-operations are based on this Struct. */
struct matrix {
  size_t columns;
  size_t rows;
  size_t size;
  value *start;
  bool diagonals;
};

/* Matrix instead of struct matrix */
typedef struct matrix matrix;

#define DOUBLE

/** Setup for the preprocessor depending on mode */
#ifdef INT
typedef int value;
#define FORMAT_STRING "%i "
#define PRECISION 0
#endif

#ifdef FLOAT
typedef float value;
#define FORMAT_STRING "%f "
#define PRECISION 0.01
#endif

#ifdef DOUBLE
typedef double value;
#define FORMAT_STRING "%f "
#define PRECISION 0.0001
#endif


/** Create a matrix */
matrix* create_matrix(size_t row, size_t col);

/** Is normally not needed for this implementation but might be needed on others */
matrix* create_zero_matrix(size_t row,size_t col);

/** Creates a identity matrix */
matrix* create_identity_matrix(size_t row,size_t col);

/** Calculate the dot product value = a • b */
value dot_product(matrix* a, matrix* b);

/** Calculate the cross product c = a x b */
void cross_product(matrix* a, matrix* b, matrix* c);

/** Destroy a matrix */
void free_matrix(matrix* mat);

/** Checks if the position exists in the matrix */
bool check_boundaries(size_t row, size_t col, matrix* mat);

/** Insert a array into the matrix */
bool insert_array(value arr[], matrix* mat);

/** Returns true if matrices a and b look the same */
bool compare_matrices(matrix* a, matrix* b);

/** Return true if the matrix are the same */
bool is_matrix(matrix* a, matrix* b);

/** Insert a value into matrix */
bool insert_value(value insert,size_t row, size_t col, matrix* mat);

/** As insert_value without check */
void insert_value_without_check(value insert, size_t row, size_t col, matrix* mat);

/** Get a value from matrix */
value get_value(size_t row, size_t col, matrix* mat);

/** As get_value without check */
value get_value_without_check(size_t row, size_t col, matrix* mat);

/** Adds a and b into c */
bool add_matrices(matrix* a, matrix* b, matrix* c);

/** Subtract a and b into c. c=a-b */
bool subtract_matrices(matrix* a, matrix* b, matrix* c);

/** Multiply a and b into c. c=a*b */
bool multiply_matrices(matrix* a, matrix* b, matrix* c);

/** Multiply a and b into c using the naive algorithm. c=a*b */
bool multiply_matrices_naive(matrix* a, matrix* b, matrix* c);

/** Multiply a and b into c. Uses row-major optimization. c=a*b */
bool multiply_matrices_optimized(matrix* a, matrix* b, matrix* c);

/** Returns the determinant of matrix a */
value get_determinant(matrix* a);

/** Calculates the inverse of a and puts it into c */
bool get_inverse(matrix* a, matrix* c);

/** Calculates the invers of a 2x2 matrix and returns it in matrix b  */
bool get_inverse_of_2x2(matrix* a,matrix* b);

/** Crout algorithm to divide matrix a into l and u that holds a=lu */
bool crout(matrix* a, matrix* l, matrix* u);

/** If no solution can be found with solve_linear, this function finds the closest one */
void least_square(matrix* a, matrix* x, matrix* b);

/** Gauss eliminates the matrix a */
bool gauss_jordan(matrix* a);

/** Solves the system of linear equations using gauss jordan */
bool gauss_jordan_solver(matrix* a,matrix* x,matrix* b);

/** Returns the lowest of the two values */
value min(value a,value b);

/** Returns on which row the largest element in the column is after start */
size_t largest_element_in_column_index(size_t column,size_t start,matrix* a);

/** Returns on which row the smallest element in the column is after start */
size_t smallest_element_in_column_index(size_t column,size_t start,matrix* a);

/** Returns on which row the first nonezero element is in the column is after start returns -1
   if no nonezero element is found */
size_t first_nonezero_in_column_index(size_t column, size_t start, matrix* a);

/** Returns on which column the first nonezero element is in the column is after start returns -1
   if no nonezero element is found */
size_t first_nonezero_in_row_index(size_t row,size_t start, matrix* a);

/** Adds each element in row1 and row 2 and puts the result on row2 */
void add_rows(size_t row1, size_t row2, matrix* a);

/** Transposes matrix a into b */
bool transpose_matrix(matrix* a, matrix*b);

/** Return the sum of a row in matrix mat */
value sum_of_row(size_t row, matrix* mat);

/** Return the sum of a column in matrix mat */
value sum_of_column(size_t column, matrix* mat);

/** Return the product of a row in matrix mat */
value product_of_row(size_t row, matrix* mat);

/** Return the product of a column in matrix mat */
value product_of_column(size_t column, matrix* mat);

/** Multiplies matrix mat with scalar */
void multiply_matrix_with_scalar(value scal, matrix* mat);

/** Divides matrix mat with scalar */
void divide_matrix_with_scalar(value scal, matrix* mat);

/** Multiplies a row with a scalar */
void multiply_row_with_scalar(value scal, size_t row, matrix* mat);

/** Divides a row with a scalar */
void divide_row_with_scalar(value scal, size_t row, matrix* mat);

/** Multiplies a column with a scalar */
void multiply_column_with_scalar(value scal, size_t col, matrix* mat);

/** Divides a column with a scalar */
void divide_column_with_scalar(value scal, size_t col, matrix* mat);

/** Takes row vector from matrix a and puts it into b */
bool get_row_vector(size_t row, matrix* a, matrix* b);

/** Inserts row vector a into b:s row */
bool insert_row_vector(size_t row, matrix* a, matrix* b);

/** Switches rows in a */
bool switch_rows(size_t row1, size_t row2, matrix* a);

/** Takes column vector from matrix a and puts it into b */
bool get_column_vector(size_t column, matrix* a, matrix* b);

/** Inserts column vector a into matrix b at position column */
bool insert_column_vector(size_t column, matrix *a, matrix* b);

/** Get a sub matrix from a */
bool get_sub_matrix(size_t start_row, size_t end_row, size_t start_col, size_t end_col, matrix* a, matrix* b);

/** Inserts the submatrix defined by start_row,end_row,start_col,end_col and put it into matrix b */
bool insert_sub_matrix(size_t start_row, size_t end_row, size_t start_col, size_t end_col, matrix* b, matrix* a);

/** Copy and return new matrix. */
matrix* matrix_copy(matrix* source);

/** Copies all the data from matrix A into matrix B */
void matrix_copy_data(matrix* A, matrix* B);

/** Checks if all elements in a matrix is equal to zero */
bool is_zero_matrix(matrix* v);

/** Checks if all elements in a matrix is positive */
bool is_non_negative_matrix(matrix* v);

/** Checks if all elements along the diagonal in a symmetric matrix is positive */
bool is_non_negative_diagonal_matrix(matrix* A);

/** Takes the diagonal in a and puts it into b */
bool get_diagonal(matrix* a,matrix* b);

/** Returns a pointer to a matrix with the derivative of var if the a matrix second order coefficiants */
matrix* derivate_matrix_with_return(size_t var,matrix* a);

/** Transforms matrix to reduced row echelon form */
void transform_to_reduced_row_echelon_form(matrix* M);

/** Return true if b contains value a */
bool matrix_contains(value a,matrix* b);

/** Compare two element values */
int compare_elements(value a, value b);

/** Returns the absolute value of a */
value matlib_fabs(value a);