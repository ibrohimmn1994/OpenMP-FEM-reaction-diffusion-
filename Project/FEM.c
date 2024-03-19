/*
 * FEM.c
 *
 *  Created on: Mar 15, 2022
 *      Author: Ibrohim
 */

 
#include <stdio.h>
#include <stdlib.h>
#include <stdbool.h>
#include <math.h>
#include "graphics.h"
#include <sys/time.h>
#include <omp.h>

static double get_wall_seconds(){
	struct timeval tv;
	gettimeofday(&tv, NULL);
	double seconds = tv.tv_sec + (double)tv.tv_usec / 1000000;
	return seconds;
}
//#############################################
// General purpose functions:
void print_matrix(int n, double** matrix){
	printf("\n");
	for (int i=0; i<n; i++){
		for (int j=0; j<n; j++){
			printf(" %f ", matrix[i][j] );
		}
		printf("\n");
	}
}
//----------------------------------------
void print_array(int n, double* array){
	printf("\n");
	for (int i=0; i<n; i++){
		printf(" %f ", array[i]);
	}
	printf("\n");
}
//----------------------------------------
void free_matrix(int n, double **mat){
        int i;
	for( i=0; i<n; i++ ){
		free(mat[i]);
	}
       
	free(mat);
	
}
//----------------------------------------
double **build_matrix(int n){
	double **matrix = (double**)malloc(n*sizeof(double*));
	int i;
   
	for ( i=0; i<n; i++ ){
		matrix[i] = (double *)malloc(n*sizeof(double));
	}
	return matrix;
}
//----------------------------------------
double *build_array(int n){
	double *array = (double*)malloc(n*sizeof(double));
	return array;
}
//----------------------------------------
double norm(int n, double *array){
	double sum = 0;
	int i;
	for ( i=0; i<n; i++){
		sum += array[i]*array[i];
	}
	return sqrt(sum);
}

//############################################################################################
// CSR functions and assemblance

#define PI (3.141592653589793)
// This struct is to store the value array, column array and index array of a particular matrix in CSR format together
typedef struct {
	double *V;
	double *C ;
	double *I ;
}assemble;
//----------------------------------------
// A slopy function that convert from CSR to regular format and it helps to check if the assemble matrices were done properly
// V stands for value array, C for columns array, I for the offset or index arrays, sizeof_V is the length of V array
// It works for the type of matrices that is used in this probelm. No confidance if it works for any sparce matrix.
double ** CSR_to_normal(int n, double *V, double *C, double *I, double sizeof_V){
	double **M = build_matrix(n);
	int vn = sizeof_V;
	int i,j;
	for (i=0; i<n; i++){
		for (j=0; j<n; j++){
			M[i][j] = 0;
		}
	}
	int ii=0;
    int t=1;
	for (i=0; i<vn; i++){
		M[ii][(int)C[i]] = V[i];
		if (i+1 == I[t]){
			ii++;
			t++;
		}
	}
	return M;
}
//----------------------------------------
// Function for matrix array multiplication in CSR format
// Referance for this function:
//https://www.it.uu.se/education/phd_studies/phd_courses/pasc/lecture-1
void CSR_matrix_array(int n, double *V, double *C, double *I, double *arr, double *result){
	register int i,j;
	register double d0;
        #pragma omp parallel for private(i,j,d0) num_threads(5)	
	for ( i=0; i<n; ++i ){
		d0 = 0.0;
		result[i] = 0;
	     for ( j=I[i]; j<I[i+1]; ++j){
	    	 d0 += V[j] * arr[(int)C[j]];
        	}
	     result[i] = d0;
	}
	
}
//######################################################################################
// Functions for assemblance
// Function to assemble the advection matrix directly into CSR format. This assemblance algorithms works for this specific problem

void Advection_CSR(int n, double **V, double **C, double **I){
	int N=(n-1)*2;
	*V = realloc(*V,N*sizeof(double));
	*C = realloc(*C,N*sizeof(double));
	*I = realloc(*I,(n+1)*sizeof(double));
	register int r;// = N/2;
	register int l;// = N-1;
	register int i;
        #pragma omp parallel for private(i,r,l) 
	for (i=N-1; i>0; i-=2){
		r= (N/2)-( (N-1-i)/2 ); // equivalent to r--   //in-dependent
		l = (N-1)-(N-i)+1;  // equivalent to l -=2
		(*V)[i]  = -0.5;
		(*V)[i-1]= 0.5;
		(*C)[i]  = r-1;
		(*C)[i-1]= r;
		(*I)[r] = l;
		//r--;         //dependent
		//l -=2;
	}
	(*I)[0]=0;
	(*I)[1]=1;
	(*I)[n]=N;
}
//----------------------------------------
// Function to assemble the diffusion matrix directly into CSR format
void Diffusion_CSR(int n, double h, double **V, double **C, double **I){
	int N = (n-1)*2+n;

	*V = realloc(*V,N*sizeof(double));
	*C = realloc(*C,N*sizeof(double));
	*I = realloc(*I,(n+1)*sizeof(double));
	register int i;
	register int r;// = n-1;
	register int l;// = N-2;

        #pragma omp parallel for private(i,r,l) 
	for (i=N-1; i>0; i-=3){
		r= (n-1)-( (N-1-i)/3 ); // equivalent to r--
	        l = (N-2)-(N-i)+1;  // equivalent to l -=3
		(*V)[i]   = 2.0/h;
		(*V)[i-1] = -1.0/h;
		(*V)[i-2] = -1.0/h;
		(*C)[i]   = r;
		(*C)[i-1] = r-1;
		(*C)[i-2] = r;
		(*I)[r] = l;
	//	r--;
	//	l-=3;
	}
	(*I)[0]=0;
	(*V)[0]=2/h;
	(*C)[0]=0;
	(*I)[n]=N;
}
//----------------------------------------
// Function to assemble the mass matrix directly into CSR format
void Mass_CSR(int n, double h, double **V, double **C, double **I){
	int N = (n-1)*2+n;
	*V = realloc(*V,N*sizeof(double));
	*C = realloc(*C,N*sizeof(double));
	*I = realloc(*I,(n+1)*sizeof(double));
	register int i;
	register int r;// = n-1;
	register int l;// = N-2;
        #pragma omp parallel for private(i,r,l) 
	for (i=N-1; i>0; i-=3){
		r= (n-1)-( (N-1-i)/3 ); // equivalent to r--
		l = (N-2)-(N-i)+1;  // equivalent to l -=3
		(*V)[i]   = 2*h/3;
		(*V)[i-1] = h/6;
		(*V)[i-2] = h/6;
		(*C)[i]   = r;
		(*C)[i-1] = r-1;
		(*C)[i-2] = r;
		(*I)[r] = l;

	}
	(*I)[0]=0;
	(*V)[0]=2*h/3;
	(*C)[0]=0;
	(*I)[n]=N;
}
//#######################################################################################
// Conjugate gradient solver:
void CG(int n, assemble A, double *b, double *x){
     double tol = 0.0001;

	double *r = build_array(n);
	double *r_old = build_array(n);
	double *p = build_array(n);
	double *w = build_array(n);
	double alpha,B;
	double normr,normp,normr_old;

	normr_old = 0;
	int i;
	for (i=0; i<n; i++){
		x[i] = 0;
		p[i]=b[i];
		r[i]=b[i];
	    normr_old += r[i]*r[i];
	}
	int k;

	for (k=0; k<n; k++){
                
	        	
                
		
		CSR_matrix_array(n,A.V,A.C,A.I,p,w); 
                #pragma omp parallel num_threads(5) 
	        {
                #pragma omp single
		{
                normp = 0;
		}
	       
                #pragma omp for private(i) reduction(+:normp)
                for (i=0; i<n; i++){
                	normp += p[i] *w[i];
                } 
                #pragma omp single
		{
		alpha = normr_old/normp;
		}
                
                #pragma omp for private(i)
		for ( i=0; i<n; i++ ){
			    r_old[i] = r[i];
				x[i] = x[i] + alpha*p[i];
				r[i] = r[i] - alpha*w[i];
		}

                #pragma omp single
		{
		normr = 0;		
		}
		//----
	
                #pragma omp for private(i) reduction(+:normr) 
		
		for (i=0; i<n; i++){
			normr += r[i]*r[i];
		}
                #pragma omp single
		{
		B = normr/normr_old;             // (r'*r)/(r_old'*r_old)
		normr_old = normr;
		}
                #pragma omp for private(i) 
		for (i=0; i<n; i++){
			p[i] = r[i] + B*p[i];
		}
               
		normr_old = normr;

		}// <<-- end of parallel region
		if (sqrt(normr) < tol){break;}
               
		
	}
	free(r);free(r_old);free(w);free(p);
}
//---------------------------------------------
// Conjugate gradient algorithmically optimized
void CG2(int n, assemble A, double *b, double *x){
	double tol = 0.0001;

	double *r = build_array(n);
	double *r_old = build_array(n);
	double *p = build_array(n);
	double *w = build_array(n);
	double alpha,B;
	double normr,normp,normr_old;

	normr_old = 0;
	int i;
	//#pragma omp parallel for private(i)
	for (i=0; i<n; i++){
		x[i] = 0;
		p[i]=b[i];
		r[i]=b[i];
	    normr_old += r[i]*r[i];
	}
        
	normp = 0;
	int k;
	for (k=0; k<n; k++){
		CSR_matrix_array(n,A.V,A.C,A.I,p,w);
                #pragma omp parallel private(alpha,B) 
		{
		
                #pragma omp single nowait
		{
                normr = 0;
		}
                #pragma omp for private(i) reduction(+:normp)
                for (i=0; i<n; i++){
			normp += p[i] *w[i];
                }
                #pragma omp single
                {
        	alpha = normr_old*1/normp;
	        }
               #pragma omp for private(i)
	       for ( i=0; i<n; i++ ){

        	        r_old[i] = r[i];
			x[i] = x[i] + alpha*p[i];
			r[i] = r[i] - alpha*w[i];
        	}
                #pragma omp single nowait
	       {
        	normp = 0;
	       }
                #pragma omp for private(i) reduction(+:normr) 
	        for (i=0; i<n; i++){
		       normr += r[i]*r[i];
	        }
                #pragma omp single
		{
          	B = normr*1/normr_old; 
                }                 		// (r'*r)/(r_old'*r_old)
                #pragma omp for private(i) 
	        for (i=0; i<n; i++){
		p[i] = r[i] + B*p[i];
	        }
                #pragma omp single nowait
		{
                normr_old = normr;
		}

		}// <<--  end of parallel region
		 //Check convergance
	        if (sqrt(normr) < tol){break;}
        	
	}
	free(r);free(r_old);free(w);free(p);
}
//#####################################################################################
// Right hand side function in the AX=b
void RHS_function(int n, double *u, assemble A, assemble S, double h, double *result){
	double epsilon = h/2;
	double *beta1 = build_array(n);
	double *beta2 = build_array(n);
	double *beta1_extra = build_array(n);
        int i;
	
	#pragma omp parallel for private(i) num_threads(3)
	for (i=0; i<n; i++){
		beta1[i] =  -0.5*u[i]*u[i];
		
	 }
        
	omp_set_nested(1);
         #pragma omp sections 
	{
        #pragma omp section

		{ CSR_matrix_array(n,S.V,S.C,S.I,u,beta2);}
		
	#pragma omp section
		{ CSR_matrix_array(n,A.V,A.C,A.I,beta1,beta1_extra);}
		
	}
	
         #pragma omp parallel for private(i)
         for (i=0; i<n; i++){
        	result[i] = beta1_extra[i]- epsilon*beta2[i];
        }

	free(beta1);free(beta2);free(beta1_extra);
}

// FE- simulation for the one-dimensional advection difussion equation
//  To use the algorithmically iotimized conjugate gradient just switch GC to GC2
void FEM(int n, int command){
	// Define variables
	double b = 2*PI;
	double a = 0;
	double h = (b-a)/(n-1);           // mesh size
	double dt = 0.01*h;                // time increment
	double T_final = 2.0;              // final time for the simulatiom
	int t_steps = round(T_final/dt);   // number of time steps
	//------------------------
	int N = n-2;
	double *x        =build_array(n);             //grid
	double *ID       =build_array(n);            // Initial_data
	double *u0_extra = build_array(N);
	double *u_new    = build_array(N);
	double *u0       = build_array(N);
         
	// Here we get the grid and the initial data
	int i;
	
	#pragma omp parallel for private(i)
	for ( i=0; i<n; i++ ){
		x[i] = a+i*h;             // independent loops
		ID[i] = sin(x[i]);
		//x[i] = x[i-1] + h;         // dependent loops
		//ID[i-1] = sin(x[i-1]);
	}

	// here we omit the boundary as the solution is fixed there and later we append the boundary to the final solution
	// we make two copied of the initial data without the boundary, the extra copy needed in the computation
	
	//#pragma omp parallel for private(i)
	for ( i=0; i< N; i++ ){
		u0[i] = ID[i+1];
		u0_extra[i]=u0[i];
	}
	//------------------------
	// assemble matrices

        assemble A; assemble S; assemble M;
       // Different matrices will have different V,C,I lengths in CSR fromat which are constructed inside function
      // and it was not possible to pass uninicialized args => they are all initialized by malloc of one double before passing them
       A.V=build_array(1); A.C=build_array(1); A.I=build_array(1);
       S.V=build_array(1); S.C=build_array(1); S.I=build_array(1);
       M.V=build_array(1); M.C=build_array(1); M.I=build_array(1);

       // Assemble matrices
       Advection_CSR(N,&A.V,&A.C,&A.I);
       Diffusion_CSR(N,h,&S.V,&S.C,&S.I);
       Mass_CSR(N,h,&M.V,&M.C,&M.I);

       // F is where we store the RHS and K is for RK4 for time discretization
       double *K=build_array(N);
       double *F=build_array(N);
    
       // Start the simulation
       int k;
       double t;
       
       
       for ( k=1; k<t_steps+1; k++){
        	t = dt*k;  // time point
       
   	        // Here are the 4 steps of RK4
        	//  the loops between the steps are to append k1, k2, k3, k4 to the u_new
	       //-------------------------------------------
	       // (1)
                
        	RHS_function(N,u0,A,S,h,F);	
	        
             	CG(N,M,F,K);
        
        
	       // here we append k1/6 to u_new and simultaniously apply u0 + k1/2 --> u0
	       //  that is used in the next step. The rest of the steps follow the same manner
               // #pragma omp parallel for private(i) 
     	       for( i=0; i<N; i++){
     	        	u_new[i] = dt/6*K[i];
     	                u0[i] = u0_extra[i] + 0.5*dt*K[i];
            	}
	
	
             	//---------------------------------------------
          	// (2)

        	RHS_function(N,u0,A,S,h,F);
        	CG(N,M,F,K);

               // #pragma omp parallel for private(i)
        	for(i=0; i<N; i++){
           	     u_new[i] += 2*dt/6*K[i];
    	             u0[i] = u0_extra[i] + 0.5*dt*K[i];
    	        }
             	//--------------------------------------------
        	//(3)

    	        RHS_function(N,u0,A,S,h,F);
          	CG(N,M,F,K);

               // #pragma omp parallel for private(i)
    	       for(i=0; i<N; i++){
    	             u_new[i] += 2*dt/6*K[i];
    	             u0[i] = u0_extra[i] + dt*K[i];
    	        }
        	//--------------------------------------------
          	//(4)
          	RHS_function(N,u0,A,S,h,F);
            	CG(N,M,F,K);

           	for(i=0; i<N; i++){
    	        	u_new[i] += dt/6*K[i] +u0_extra[i];  //<<-- u_new = u0 + (k1+2*k2+2*k3+k4)/6  which is the solution
    	                u0[i] = u_new[i];
    	                u0_extra[i] = u_new[i];
    	         }
	
	
        }
    
        // Add the boundary to the solution
        //-------------------------------------------
        u_new = (double*)realloc(u_new,n*sizeof(double));
        
        for(i=n-1; i>0; i--){
    	     u_new[i] = u_new[i-1];
        }
        u_new[0] = ID[0]; u_new[n-1] = ID[n-1];
	
        //---------------------------------------------
         printf("At time point %f the solution is:\n",t);
        //uncoment to print the final solution:
	//for n = 11
	//The solution should be this for refereance
	//0.000000  0.188349  0.364064  0.481468  0.391387  -0.000000  -0.391387  -0.481468  -0.364064  -0.188349  -0.000000
       //   print_array(n,u_new);

        // To view the plot of the final solution and verify if the simulation went correct
    
          if (command){
        	Draw(n,x,u_new);
          }
	  


          free(x);free(ID);free(u0);free(u0_extra);free(K);free(F);
          free(M.V);free(M.C);free(M.I);
	  free(A.V);free(A.C);free(A.I);
	  free(S.V);free(S.C);free(S.I);
}

//##############################################################

// the program is executed with ./FEM N1 N2
// N1 is the size of the mesh, preferably an odd number
// N2 is command to get the plot if 0 => no plot, if 1 => plot of the final solution
int main(int argc, char *argv[]){
	if (argc != 3){printf("The program has two input N1: mesh-size, N2:command for the plot either 0 or one");}
	else{
	int n = atoi(argv[1]);
	int command = atoi(argv[2]);

	//---------------------------------------
	//double time = get_wall_seconds();
	double time = omp_get_wtime();
	FEM(n,command);
	time = omp_get_wtime() - time;
	//time = get_wall_seconds() - time;
	printf("\n Time of the simulation is: %f\n", time);
	}
	//--------------------------------------

	printf("\nThe end\n");
	return 0;
}

