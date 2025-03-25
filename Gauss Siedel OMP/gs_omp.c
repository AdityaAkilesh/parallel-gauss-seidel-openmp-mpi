#include <stdio.h>
#include <math.h>
#include <stdlib.h>
#include <sys/time.h>
#include <omp.h>
#include <unistd.h>

/*Global Variables declaration*/
#define tolerance 10e-7
#define maximum_iterations 100
double **A, **B, difference;

void print(double **A, int n) /*This function used print the matrixes used for debugging*/
{
    for(int i=0; i<n+1; i++){
        for(int j=0; j<n+1; j++){
            printf("%f", A[i][j]);
        }
        printf("\n");
    }
    printf("\n");
}

void initialize_matrix(double **A, int n) //initiate the matrix with top and left with and remaining with zeros
 {
    int i, j ;
    for(j=0; j<n+1; j++){
        A[0][j]=1.0;        //allocating ones

    }
    for (i=1;i<n+1;i++){
      A[i][0]=1.0;
      for (j=1;j<n+1;j++) A[i][j]=0.0; //remaining elements zeros 
    }     
   
}

void serial_solve(double **B, int n ){  // this Fuction used solve system sequentially gauss seidel serial method 
    printf("\n**************Serial Solver*************"); 
    int convergence =0; 
    double difference, temp; // to know difference and temp used for swaping 
    int iterations=0; // to get ietations count 
    int curnt_iterations; // variable current iteration 
    int i,j; 
    for(curnt_iterations=1; curnt_iterations<=maximum_iterations; curnt_iterations++){
        difference=0.0;
        for(i=1; i<=n; i++){
            for(j=1; j<=n; j++){
                temp = B[i][j];
                B[i][j] = 0.2*(B[i][j] + B[i][j-1] + B[i-1][j] + B[i][j+1] + B[i+1][j]); //calulatinh the values of neighbors and updating 
                difference += fabs(B[i][j] - temp); 
            }
        }
        iterations++;
        printf("Difference after %3d iterations: %f\n", iterations, difference);
        if (difference/((double)n*(double)n) < tolerance) // Checking the if the difference is less than tolerance 
        {
            printf("\nConvergence achieved %d iterations....Now exiting\n\n", iterations);
            return;
        }        
    }
    printf("\n\nMaximum Iterations Reached...exit out Exiting\n\n"); // if covergence not reached then maximum iterations reached exit out 
}

void parallel_red_black_solve(double **A, int n, int thread_count){
    printf("\n*************Red Black Solver*************");
    int convergence=0;
    int curnt_iterations;      
    int iterations=0;
    double difference, temp;
    int i,j;
    /* varibles  decalred for calculating*/
    for(curnt_iterations=1; curnt_iterations<=maximum_iterations; curnt_iterations++)
    {
        difference=0;
        #pragma omp parallel num_threads(thread_count) private(temp, i, j) reduction(+:difference) // initiate the openmp and thread count reduction 
        {
            #pragma omp for 
            for(i=1; i<=n; i++)
            {
                for(j=1; j<=n; j++) 
                {                                     //check the codition odd and even 
                    if((i+j)%2==1)      //checking odd or black 
                    {
                        temp = A[i][j];
                        A[i][j] = 0.2*(A[i][j] + A[i][j-1] + A[i-1][j] + A[i][j+1] + A[i+1][j]); // update the values of neighbors 
                        difference += fabs(A[i][j] - temp);
                    }
                }
            }
            #pragma omp barrier
        }
        #pragma omp parallel num_threads(thread_count) private(temp, i, j) reduction(+:difference)
        {
            #pragma omp for
            for(i=1; i<=n; i++){
                for(j=1; j<=n; j++){
                    if((i+j)%2==0){             //check the even or reb ordering 
                        temp = A[i][j];
                        A[i][j] = 0.2*(A[i][j] + A[i][j-1] + A[i-1][j] + A[i][j+1] + A[i+1][j]);
                        difference += fabs(A[i][j] - temp);
                        
                    }
                }
            }
            #pragma omp barrier            
        }
        iterations++;
        printf("difference reduced to each after %d iterations: %f\n", iterations, difference);
        if(difference/((double)n*(double)n)<tolerance)
        {
            printf("\nConvergence  after %d iterations.....\n\n", iterations);
            return;
        }

       
    }
    printf("\nMaximum Iterations Reached......exit out ..exiting");

}

int main(int argc, char *argv[]){
    int i; 
     

    if (argc != 3) {
        printf("Usage: %s <thread_count> <N>\n", argv[0]);
        exit(1);
    }
    const int thread_count = strtol(argv[1], NULL, 10); //thread count input from command line as an arguement 
    const int N = strtol(argv[2],NULL, 10); // Also, matrix size input from commad line as an arguement 
    printf("thread_count = %d\n", thread_count);
    printf("N = %d\n", N);


    double start_time;
    double end_time; //declaration of start time and end time for time calculation 

    
    A = (double **)malloc((N+2) * sizeof(double *)); //allocation of memory of matrix A
    B = (double **)malloc((N+2) * sizeof(double *)); //allocation of memory of matrix  B

    for(i=0; i<N+2; i++)
    {

        A[i] = (double*)malloc((N+2)*sizeof(double));
        B[i] = (double*)malloc((N+2)*sizeof(double));    
    }
    srand(time(NULL)); 
    initialize_matrix(B, N); // initiate matrix B
    start_time = omp_get_wtime();
    serial_solve(B,N); // call serial function
    end_time = omp_get_wtime();

    double total_time_serial;
    total_time_serial = end_time - start_time;
    printf("\ntotal Execution time Serial or Sequential execution: %f\n\n", total_time_serial);

    initialize_matrix(A,N); // initate the matrix A 

    start_time = omp_get_wtime();
    parallel_red_black_solve(A, N, thread_count); //call the function parallel red black 
    end_time = omp_get_wtime();

    double total_time_parallel;
    total_time_parallel = end_time-start_time;

    printf("\ntotal execution time parallel execution:%f\n\n", total_time_parallel);
    
	FILE *f;    //declaration for accessin the file 
	if (access("data.csv", F_OK) == -1) 
    {
 	f = fopen("data.csv", "a");
	fprintf(f, "thread_count;MatrixSize;Total_time_serial;total_time_parallel\n");
	}
	else 
    {
		f = fopen("data.csv", "a");
	}

	fprintf(f, "%d;%d;%f;%f\n", thread_count, N, total_time_serial, total_time_parallel);
	fclose(f);


}