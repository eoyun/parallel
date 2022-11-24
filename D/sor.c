#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>


int main(){
	int Nbin = 100;
	int i,j,iter;
	matrix U = matrix_create(Nbin,Nbin);
	//matrix U_tmp = matrix_create(Nbin,Nbin);
	matrix F = matrix_create(Nbin,Nbin);
	float f_ele;
	float x,y;
	float u_tmp, u1, u2, u3, u4, u;
	const float PI = 3.141592;	
	//FILE *fp;
	//char txtname[100];
	float w=1.8;
	float acc;
	//initial condition ; all zeros
	for(i=0;i<Nbin+1;i++){
		for(j=0;j<Nbin+1;j++){
			matrix_set_element(&U,i,j,1);
			x = (float) i/Nbin;	
			y = (float) j/Nbin;	
			f_ele = sinf(x * PI)*cosf(y * PI);
		//	printf("x is %f, y is %f ; f is %f\n",x,y,f_ele);
			matrix_set_element(&F,i,j,f_ele);	
		}
	}
	boundary_condition_D(&U);
	acc = accuracy(&U);
	iter=0;
	while (fabs(acc)>0.000001){
	//for (iter=0;iter<5;iter++){
		//sprintf(txtname,"gauss%d.txt",iter);
		//fp = fopen(txtname,"wt");
		
	/*	for(i=0;i<Nbin+1;i++){
			for(j=0;j<Nbin+1;j++){
				u_tmp = get_element(&U,i,j);
				matrix_set_element(&U_tmp,i,j,u_tmp);
			}
		}*/
		for(i=1;i<Nbin;i++){
			for(j=1;j<Nbin;j++){
				//if (i==0) u1 = get_element(&U,i,j); 
				u1 = get_element(&U,i-1,j);
				//if (j==0) u2 = get_element(&U,i,j);
				u2 = get_element(&U,i,j-1);
				//if (i==Nbin) u3 = get_element(&U,i,j);
				u3 = get_element(&U,i+1,j);
				//if (j==Nbin) u4 = get_element(&U,i,j);
				u4 = get_element(&U,i,j+1);
				u = get_element(&U,i,j);
				f_ele = get_element(&F,i,j);
				u_tmp = u + w * (u1 + u2 + u3 + u4 - f_ele/Nbin/Nbin - u * 4)/4;
				
				matrix_set_element(&U,i,j,u_tmp);
				//fprintf(fp,"%f ",u_tmp);
			}
			//fprintf(fp,"\n");
		}
		printf("iteration number is %d, accuracy is %f\n",iter,acc);
		boundary_condition_D(&U);
		acc = accuracy(&U);
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
