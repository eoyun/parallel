#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

typedef struct {
	int row, col;
	double *data;
} matrix;


matrix matrix_create (int row, int col);
matrix matrix_copy (matrix *matrix_);
void boundary_condition_D(matrix *matrix_);
void boundary_condition_N(matrix *matrix_);
double get_element (matrix *matrix_, int row, int col);
void matrix_set_element (matrix *matrix_,int row, int col, double data);
void matrix_delete (matrix *matrix_);
double accuracy (matrix *matrix_);
double accuracy_N (matrix *matrix_);
double accuracy_compare (matrix *M1,matrix *M2);
