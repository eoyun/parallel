#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>


matrix matrix_create (int row, int col){
	matrix matrix_;
	matrix_.row = row;
	matrix_.col = col;
	matrix_.data =(float*) calloc((row+1) * (col+1), sizeof(float));
	return matrix_;

}

matrix matrix_copy (matrix *matrix_){
	matrix M1 = matrix_create(matrix_->row,matrix_->col);
	M1.data = matrix_->data;
	return M1;
}
void boundary_condition_D(matrix *matrix_){
	int i;
	for (i=0;i<matrix_->row+1;i++){ matrix_set_element(matrix_,i,0,0); matrix_set_element(matrix_,i,matrix_->col,0);}
	for (i=1;i<matrix_->col;i++){ matrix_set_element(matrix_,0,i,0); matrix_set_element(matrix_,matrix_->row,i,0);}

}
void boundary_condition_N(matrix *matrix_){
	int i;
	for (i=0;i<matrix_->row+1;i++){ 
		matrix_set_element(matrix_,i,0,get_element(matrix_,i,1)); 
		matrix_set_element(matrix_,i,matrix_->col,get_element(matrix_,i,matrix_->col-1));
	}
	for (i=1;i<matrix_->col;i++){ matrix_set_element(matrix_,0,i,0); matrix_set_element(matrix_,matrix_->row,i,0);}
}

float get_element (matrix *matrix_, int row, int col){
	return matrix_->data[row*matrix_->col+col];
	//return 3;
}

void matrix_set_element (matrix *matrix_,int row, int col, float data){
	matrix_->data[row*matrix_->col+col] = data;
}

void matrix_delete (matrix *matrix_){
	free(matrix_->data);
}

float accuracy (matrix *matrix_){
	int i,j;
	float u_n, u_a;
	float PI = 3.141592;
	float x,y;
	float sum=0;
	float acc;

	for (i=0;i<matrix_->row+1;i++){
		for (j=0;j<matrix_->col+1;j++){
			u_n = get_element(matrix_,i,j);
			x = (float) i/matrix_->row;
			y = (float) j/matrix_->col;
			u_a = -sinf(x * PI) * cosf(y * PI)/2/PI/PI;
			sum += (u_n - u_a) * (u_n - u_a);
			//printf("x is %f, y is %f, analytic sol is %f\n",x,y,u_a);
		}
	}
	acc = pow(sum/(matrix_->row+1)/(matrix_->col+1),0.5);
	return acc;
}
float accuracy_compare (matrix *M1, matrix *M2){
	int i,j;
	float u_1, u_2;
	float PI = 3.141592;
	//float x,y;
	float sum=0;
	float acc;

	for (i=0;i<M1->row+1;i++){
		for (j=0;j<M2->col+1;j++){
			u_1 = get_element(M1,i,j);
			u_2 = get_element(M2,i,j);
			sum += (u_1 - u_2) * (u_1 - u_2);
			//printf("x is %f, y is %f, analytic sol is %f\n",x,y,u_a);
		}
	}
	acc = pow(sum/(M1->row+1)/(M1->col+1),0.5);
	return acc;
}
//matrix matrix_add(matrix *M1, matrix *M2){

//}

//matrix matrix_multiply(matrix *M1, matrix *M2){
//	matrix M3;
	
//}
/*
int main(){
	int c = 3;
	int d = 4;
	matrix a = matrix_create(c,d);
	//a.data[(1)*d+d-1]=1;
	//a.data[1]=0;
	matrix_set_element(&a,1,2,3.0);
	matrix_set_element(&a,1,3,4.5);
	printf("element is %f\n",matrix_get_element(&a,1,2));
	printf("element is %f\n",matrix_get_element(&a,1,3));
	return 0;
}*/
