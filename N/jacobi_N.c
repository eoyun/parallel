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
	matrix U_tmp = matrix_create(Nbin_x,Nbin_y);
	matrix F = matrix_create(Nbin_x,Nbin_y);
	double f_ele;
	double x,y;
	double u_tmp, u1, u2, u3, u4;
	const double PI = 3.141592;
	double acc, acc_com;	
	//FILE *fp;
	//char txtname[100];
	//initial condition ; all zeros
	or(i=0;i<Nbin_x+1;i++){
		for(j=0;j<Nbin_y+1;j++){
			matrix_set_element(&U,i,j,i*j);
			x = (double) i/Nbin;	
			y = (double) (j-1)/Nbin;	
			f_ele = sinf(x * PI)*cosf(y * PI);
		//	printf("x is %f, y is %f ; f is %f\n",x,y,f_ele);
			matrix_set_element(&F,i,j,f_ele);	
		}
	}
	boundary_condition_N(&U);
	acc = accuracy(&U);
	iter =0;
	while (fabs(acc)>0.0005){
	//for (iter=0;iter<800;iter++){
		//sprintf(txtname,"jacobi%d.txt",iter);
		//fp = fopen(txtname,"wt");
		
		for(i=0;i<Nbin_x+1;i++){
			for(j=0;j<Nbin_y+1;j++){
				u_tmp = get_element(&U,i,j);
				matrix_set_element(&U_tmp,i,j,u_tmp);
			}
		}
		for(i=1;i<Nbin_x;i++){
			for(j=1;j<Nbin_y;j++){
				//if (i==0) u1 = get_element(&U_tmp,i,j); 
				u1 = get_element(&U_tmp,i-1,j);
				//if (j==0) u2 = get_element(&U_tmp,i,j);
				u2 = get_element(&U_tmp,i,j-1);
				//if (i==Nbin) u3 = get_element(&U_tmp,i,j);
				u3 = get_element(&U_tmp,i+1,j);
				//if (j==Nbin) u4 = get_element(&U_tmp,i,j);
				u4 = get_element(&U_tmp,i,j+1);
				f_ele = get_element(&F,i,j);
				u_tmp = (u1 + u2 + u3 + u4 - f_ele/Nbin_x/Nbin_y)/4;
				
				matrix_set_element(&U,i,j,u_tmp);
				//fprintf(fp,"%f ",u_tmp);
			}
			//fprintf(fp,"\n");
		}
		boundary_condition_N(&U);
		acc_com = accuracy_compare(&U,&U_tmp);
		printf("iteration number is %d, accuracy is %f, com accuracy is %f\n",iter,acc,acc_com);
		if (iter==100000) break;
		acc = accuracy(&U);
		iter+=1;
		//printf("this is dbg line; end of loop not close file\n",i,j);
		//fclose(fp);
	}
	//printf("this is dbg line; end of loop\n",i,j);
	matrix_delete(&U);			
	matrix_delete(&U_tmp);			
	matrix_delete(&F);			
	return 0;
}
