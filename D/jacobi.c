#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>


int main(){
	int Nbin = 40;
	double h = (double)1/Nbin;
	int i,j,iter;
	matrix U = matrix_create(Nbin,Nbin);
	matrix U_tmp = matrix_create(Nbin,Nbin);
	matrix F = matrix_create(Nbin,Nbin);
	double f_ele;
	double x,y;
	double u_tmp, u1, u2, u3, u4;
	const double PI = 3.141592;
	double acc, acc_com;	
	//FILE *fp;
	//char txtname[100];
	//initial condition ; all zeros
	for(i=0;i<Nbin+1;i++){
		for(j=0;j<Nbin+1;j++){
			matrix_set_element(&U,i,j,i*j);
			x = (double) i/Nbin;	
			y = (double) j/Nbin;
			//printf("%d ",i*100+j);
			//if (i==0 || i==Nbin||j==0||j==Nbin) matrix_set_element(&U,i,j,-sin(x*M_PI)*cos(y*M_PI)/M_PI/M_PI/2);	
			f_ele = sin(x * M_PI)*cos(y * M_PI);
		//	printf("x is %f, y is %f ; f is %f\n",x,y,f_ele);
			matrix_set_element(&F,i,j,f_ele);	
		}
		//printf("\n");
	}
	boundary_condition_D(&U);
	for (i=0;i<Nbin+1;i++){
		for (j=0;j<Nbin+1;j++){
			printf("%f ",get_element(&U,i,j));
		} 
		printf("\n");
	}
	acc = accuracy(&U);
	iter =0;
	//while (fabs(acc)>0.00001){
	while (iter<10000){
	//for (iter=0;iter<800;iter++){
		//sprintf(txtname,"jacobi%d.txt",iter);
		//fp = fopen(txtname,"wt");
		
		for(i=0;i<Nbin+1;i++){
			for(j=0;j<Nbin+1;j++){
				u_tmp = get_element(&U,i,j);
				matrix_set_element(&U_tmp,i,j,u_tmp);
			}
		}
		for(i=1;i<Nbin;i++){
			for(j=1;j<Nbin;j++){
				//if (i==0) u1 = get_element(&U_tmp,i,j); 
				u1 = get_element(&U_tmp,i-1,j);
				//if (j==0) u2 = get_element(&U_tmp,i,j);
				u2 = get_element(&U_tmp,i,j-1);
				//if (i==Nbin) u3 = get_element(&U_tmp,i,j);
				u3 = get_element(&U_tmp,i+1,j);
				//if (j==Nbin) u4 = get_element(&U_tmp,i,j);
				u4 = get_element(&U_tmp,i,j+1);
				f_ele = get_element(&F,i,j);
				u_tmp = (u1 + u2 + u3 + u4 - f_ele*h*h)/4;
				
				matrix_set_element(&U,i,j,u_tmp);
				//fprintf(fp,"%f ",u_tmp);
			}
			//fprintf(fp,"\n");
		}
		for (i=0;i<Nbin+1;i++){
			for (j=0;j<Nbin+1;j++){
				printf("%f ",get_element(&U,i,j));
			} 
			printf("\n");
		}
		//boundary_condition_D(&U);
		acc_com = accuracy_compare(&U,&U_tmp);
		acc = accuracy(&U);
		printf("iteration number is %d, accuracy is %f, com accuracy is %f\n",iter,acc,acc_com);
		if (iter==100000) break;
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
