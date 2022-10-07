#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>

typedef struct {
	int row, col;
	float *data;
} matrix;


matrix matrix_create (int row, int col);
matrix matrix_copy (matrix *matrix_);
void boundary_condition_D(matrix *matrix_);
void boundary_condition_N(matrix *matrix_);
float get_element (matrix *matrix_, int row, int col);
void matrix_set_element (matrix *matrix_,int row, int col, float data);
void matrix_delete (matrix *matrix_);
float accuracy (matrix *matrix_);
float accuracy_compare (matrix *M1,matrix *M2);
