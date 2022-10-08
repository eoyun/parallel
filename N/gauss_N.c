#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>


int main(){
	int Nbin = 100;
	int Nbin_x = Nbin;
	int Nbin_y = Nbin+2;
	int i,j,iter;
	matrix U = matrix_create(Nbin_x,Nbin_y);
	//matrix U_tmp = matrix_create(Nbin,Nbin);
	matrix F = matrix_create(Nbin_x,Nbin_y);
	float f_ele;
	float x,y;
	float u_tmp, u1, u2, u3, u4;
	const float PI = 3.141592;	
	FILE *fp;
	float acc;
	char txtname[100];
	
	//initial condition ; all zeros
	for(i=0;i<Nbin_x+1;i++){
		for(j=0;j<Nbin_y+1;j++){
			matrix_set_element(&U,i,j,i*j);
			x = (float) i/Nbin;	
			y = (float) (j-1)/Nbin;	
			f_ele = sinf(x * PI)*cosf(y * PI);
		//	printf("x is %f, y is %f ; f is %f\n",x,y,f_ele);
			matrix_set_element(&F,i,j,f_ele);	
		}
	}
	boundary_condition_N(&U);
	//for (iter=0;iter<5;iter++){
	acc = accuracy_N(&U);
	iter = 0;
	printf("hello?\n");
	while(fabs(acc)>0.0001){
		sprintf(txtname,"gauss%d.txt",iter);
		//fp = fopen(txtname,"wt");
		
	/*	for(i=0;i<Nbin+1;i++){
			for(j=0;j<Nbin+1;j++){
				u_tmp = get_element(&U,i,j);
				matrix_set_element(&U_tmp,i,j,u_tmp);
			}
		}*/
		for(i=1;i<Nbin_x;i++){
			for(j=1;j<Nbin_y;j++){
				//if (i==0) u1 = get_element(&U,i,j); 
				u1 = get_element(&U,i-1,j);
				//if (j==0) u2 = get_element(&U,i,j);
				u2 = get_element(&U,i,j-1);
				//if (i==Nbin) u3 = get_element(&U,i,j);
				u3 = get_element(&U,i+1,j);
				//if (j==Nbin) u4 = get_element(&U,i,j);
				u4 = get_element(&U,i,j+1);
				f_ele = get_element(&F,i,j);
				u_tmp = (u1 + u2 + u3 + u4 - f_ele/Nbin_x/Nbin_y)/4;
				
				matrix_set_element(&U,i,j,u_tmp);
				//fprintf(fp,"%f ",u_tmp);
			}
			//fprintf(fp,"\n");
		}
		//acc_com = accuracy_compare(&U,&U_tmp);
		printf("iteration number is %d, accuracy is %f\n",iter,acc);
		boundary_condition_N(&U);
		acc = accuracy_N(&U);

		iter+=1;
		if (iter>100000){printf("something wrong\n"); break;}
		//printf("this is dbg line; end of loop not close file\n",i,j);
		//fclose(fp);
	}
	//printf("this is dbg line; end of loop\n",i,j);
	matrix_delete(&U);			
	//matrix_delete(&U_tmp);			
	matrix_delete(&F);			
	return 0;
}
