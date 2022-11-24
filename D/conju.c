#include <stdio.h>
#include <malloc.h>
#include <stdlib.h>
#include "matrix.h"
#include <math.h>


int main(){
	int Nbin = 100;
	int i,j;
	matrix X = matrix_create(Nbin,Nbin);
	matrix R = matrix_create(Nbin,Nbin);
	matrix R_tmp = matrix_create(Nbin,Nbin);
	matrix P = matrix_create(Nbin,Nbin);
	double alpha;
	double alpha_deno=0.;
	double alpha_numer=0.;
	double beta;
	double beta_deno=0.;
	double beta_numer=0.;
	double f_ele;
	double x,y;
	const double PI = 3.141592;	
	matrix F = matrix_create(Nbin,Nbin);
	double x_tmp, x1, x2, x3, x4, x0;
	double r, p;
	double p_square=0.;
	int iter=0;
	double acc,acc_no_com;
	// define X0 and b
	for(i=0;i<Nbin+1;i++){
		for(j=0;j<Nbin+1;j++){
			matrix_set_element(&X,i,j,i*j);
			x = (double) i/Nbin;	
			y = (double) j/Nbin;	
			f_ele = sin(x * M_PI)*cos(y * M_PI);
		//	printf("x is %f, y is %f ; f is %f\n",x,y,f_ele);
			matrix_set_element(&F,i,j,f_ele);	
		}
	}
	boundary_condition_D(&X);
 	// define P0 and R0
	for(i=1;i<Nbin;i++){
        	for(j=1;j<Nbin;j++){
                	//if (i==0) x1 = get_element(&X,i,j);
                	if (i==0) x1 = 0;
                        else x1 = get_element(&X,i-1,j);
                        //if (j==0) x2 = get_element(&X,i,j);
                        if (j==0) x2 = 0;
                        else x2 = get_element(&X,i,j-1);
                        //if (i==Nbin) x3 = get_element(&X,i,j);
                        if (i==Nbin) x3 = 0;
                        else x3 = get_element(&X,i+1,j);
			//if (j==Nbin) x4 = get_element(&X,i,j);
                        if (j==Nbin) x4 = 0;
                        else x4 = get_element(&X,i,j+1);
			x0 = get_element(&X,i,j);
                        f_ele = get_element(&F,i,j);
                        r = (4*x0 - x1 - x2 - x3 - x4)*Nbin*Nbin + f_ele;
                        //r = (4*x0 - x1 - x2 - x3 - x4) + f_ele;
			
                        matrix_set_element(&R,i,j,r);
                        matrix_set_element(&P,i,j,r);
                }
        }
	// alpha
	for (i=1;i<Nbin;i++){
		for (j=1;j<Nbin;j++){
			alpha_numer += get_element(&R,i,j)*get_element(&R,i,j);
                	//if (i==0) x1 = get_element(&P,i,j);
                	if (i==0) x1 = 0;
                        else x1 = get_element(&P,i-1,j);
                        //if (j==0) x2 = get_element(&P,i,j);
                        if (j==0) x2 = 0;
                        else x2 = get_element(&P,i,j-1);
                        //if (i==Nbin) x3 = get_element(&P,i,j);
                        if (i==Nbin) x3 = 0;
                        else x3 = get_element(&P,i+1,j);
                        //if (j==Nbin) x4 = get_element(&P,i,j);
                        if (j==Nbin) x4 = 0;
                        else x4 = get_element(&P,i,j+1);
			x0 = get_element(&P,i,j);
			p_square = (x1+x2+x3+x4-4*x0)*Nbin*Nbin;
			//p_square = (x1+x2+x3+x4-4*x0);
			alpha_deno += get_element(&P,i,j)*p_square;
		}
	}
	alpha = alpha_numer/alpha_deno;
	// iteration start 1.x 2.r 3.p 4.alpha & beta
	acc=accuracy(&X);
	while(acc>0.000001){
		for (i=1;i<Nbin;i++){
			for (j=1;j<Nbin;j++){
				x_tmp = get_element(&X,i,j) + alpha * get_element(&P,i,j);
				matrix_set_element(&X,i,j,x_tmp);
				printf("%f ",x_tmp);
									
			}
			printf("\n");
		}
		boundary_condition_D(&X);
	       	beta_numer =0.;	
	       	beta_deno =0.;	
		for (i=1;i<Nbin;i++){
			for (j=1;j<Nbin;j++){
				matrix_set_element(&R_tmp,i,j,get_element(&R,i,j));	
		        	//if (i==0) x1 = get_element(&X,i,j);
		        	if (i==0) x1 = 0;
		                else x1 = get_element(&P,i-1,j);
		                //if (j==0) x2 = get_element(&X,i,j);
		                if (j==0) x2 = 0;
		                else x2 = get_element(&P,i,j-1);
		                //if (i==Nbin) x3 = get_element(&X,i,j);
		                if (i==Nbin) x3 = 0;
		                else x3 = get_element(&P,i+1,j);
		                //if (j==Nbin) x4 = get_element(&X,i,j);
		                if (j==Nbin) x4 = 0;
		                else x4 = get_element(&P,i,j+1);
				x0 = get_element(&P,i,j);
		                //f_ele = get_element(&F,i,j);
		                r = alpha*(4*x0 - x1 - x2 - x3 - x4)*Nbin*Nbin + get_element(&R,i,j);
		                //r = (4*x0 - x1 - x2 - x3 - x4) + f_ele;
				matrix_set_element(&R,i,j,r);
				//printf("%f\n",r);
				beta_deno += get_element(&R_tmp,i,j) * get_element(&R_tmp,i,j);
				beta_numer += r * r;
			}
		}
		beta = beta_numer/beta_deno;
		alpha_numer = 0.;
		alpha_deno = 0.;
		for (i=1;i<Nbin;i++){
			for (j=1;j<Nbin;j++){
				p = get_element(&R,i,j) + beta * get_element(&P,i,j);
				matrix_set_element(&P,i,j,p);
			}
		}	
		for (i=1;i<Nbin;i++){
			for (j=1;j<Nbin;j++){
				alpha_numer += get_element(&R,i,j)*get_element(&R,i,j);
		        	//if (i==0) x1 = get_element(&P,i,j);
		        	if (i==0) x1 = 0;
		                else x1 = get_element(&P,i-1,j);
		                //if (j==0) x2 = get_element(&P,i,j);
		                if (j==0) x2 = 0;
		                else x2 = get_element(&P,i,j-1);
		                //if (i==Nbin) x3 = get_element(&P,i,j);
		                if (i==Nbin) x3 = 0;
		                else x3 = get_element(&P,i+1,j);
		                //if (j==Nbin) x4 = get_element(&P,i,j);
		                if (j==Nbin) x4 = 0;
		                else x4 = get_element(&P,i,j+1);
				x0 = get_element(&P,i,j);
				p_square = (x1+x2+x3+x4-4*x0)*Nbin*Nbin;
				//p_square = (x1+x2+x3+x4-4*x0);
				alpha_deno += get_element(&P,i,j)*p_square;
				
			}
		}
		alpha = alpha_numer/alpha_deno;
		acc= accuracy(&X);
		printf("iteration number is %d, accuracy is %f,alpha is %f, beta is %f\n",iter,acc,alpha,beta);
		//printf("alpha numer is %f, alpha deno is %f, beta numer is %f, beta deno is %f\n",alpha_numer,alpha_deno,beta_numer,beta_deno);
		iter+=1;
		//if (iter ==10000) break;	
		//acc=accuracy_compare(&R,&R_tmp);
	}
	//iteration end;




	matrix_delete(&X);			
	matrix_delete(&R);			
	matrix_delete(&R_tmp);			
	matrix_delete(&P);			
	return 0;
}
