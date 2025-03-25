#include <stdio.h>
#include <stdlib.h>
#include <unistd.h>
#include <math.h>
#include "mpi.h"

/*Global Variables declaration*/
#define maximum_iterations 500
#define tolerance 10e-7
#define MAX 10

/*global declaration for functions*/

float generate_random(int max); // This function generate float numbers with value of max
int find_rows_with_max_values(int number_nodes, int n); // This function calculates maximu how many rows are given to each node
int find_node_position(int nodex_index, int n, int max_num_rows); // This function gets the position from which elements are gonna sent or received
int find_node_elements(int node_index, int n, int max_num_rows); // This calculates how many elements are going to a given node
 
void initialize_master_matrix(float **A, int n, int m);  //initialize 2D matrx in the master node

void initialize_node_matrix(float **A, int number_elements); // intialize 2D matrix in slave nodes with submatrixes



void Gauss_seidel_solver (float **A, int n, int number_elements);

int main(int argc, char *argv[]){
    int comm_sz, my_rank;
    int n, i;
    float *A, *b;

    /*initilaize MPI and take the input from command line arguement*/

    MPI_Init(&argc, &argv);
	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	MPI_Comm_size(MPI_COMM_WORLD, &comm_sz);
    if (argc < 2) {
		if (my_rank == 0) {
			printf("please enter matrix size: matrix_size \n");
			printf("\t matrix_size shoul be power of 2 (e.g. : 16, 1024)\n");
				
		}

		MPI_Finalize();
		exit(1);
	}
    
    n = atoi(argv[1]);

    if(my_rank==0){
        printf("Matrix size = %d\n", n);
    }
     

    int max_num_rows = find_rows_with_max_values(comm_sz, n);

    //Calculate array conatining the offset and numner of elements per node
    int nodes_offsets[comm_sz];
	int nodes_elements[comm_sz];
	for (i = 0; i < comm_sz; i++) {
		nodes_offsets[i] = find_node_position(i, n, max_num_rows);
		nodes_elements[i] = find_node_elements(i, n, max_num_rows);
    }
    // assign a variable to store the node local number of elements

    int number_elements = nodes_elements[my_rank];

    //Initialize the 2D total matrix in the master node
    if(my_rank == 0) {
		
			initialize_master_matrix(&A, n, n);		
    }
    //simlarly, allocate node matrix where the reciving rows are computed
    //stored in number_elements
    initialize_node_matrix(&b, number_elements);

    // Use Scatter function to broadcast the values on all processors
    MPI_Scatterv(A, nodes_elements, nodes_offsets, MPI_FLOAT, b, number_elements, MPI_FLOAT, 0,MPI_COMM_WORLD);  

    double start_time = MPI_Wtime();

    //*************************Gauss Seidel****************
    Gauss_seidel_solver(&b, n, number_elements);

    double end_time = MPI_Wtime();
    //Gather all the values
    MPI_Gatherv(b, number_elements, MPI_FLOAT, A, nodes_elements, nodes_offsets, MPI_FLOAT, 0, MPI_COMM_WORLD);
    if(my_rank==0){
        printf("Total execution time = %lf seconds\n", end_time - start_time);
        free(A);
    }	
	if (my_rank == 0) {
		double total_time_elapsed = end_time - start_time; 
		printf("Total time: %f\n", total_time_elapsed);

		FILE *f;
		if (access("data.csv", F_OK) == -1) {
 			f = fopen("data.csv", "a");
			fprintf(f, "Processes;Size;Total-time;\n");
		}
		else {
			f = fopen("data.csv", "a");
		}

		fprintf(f, "%d;%d;%f;\n", comm_sz, n, total_time_elapsed);
		fclose(f);
	}
	//free(A);
	free(b);

	MPI_Finalize();
	return 0;
 }

void Gauss_seidel_solver(float **A, int n, int number_elements) {

	
	int convergence = 0, iterations = 0, my_rank;
    float difference = 0, temp;

  	while (!convergence && ( iterations< maximum_iterations)) {
  		difference = 0;

  		// the first row or the last row are not solved these are updated by each ieteration
  		//  it starts at "n" and it goes up to "number_elements - 2n")
  		for (int i = n; i < number_elements -(2*n); i++) {
  			// (So the first and last positions of "rows" are skipped)
  			if ((i % n == 0) || ((i+1) % n == 0)) {
				continue;
			}

  			int north = i - n;
  			int south = i + n;
  			int west = i - 1;
  			int east = i + 1;

  			temp = (*A)[i];
			(*A)[i] = 0.25 * ((*A)[i] + (*A)[east] + (*A)[north] + (*A)[west] + (*A)[south]);
			difference += abs((*A)[i] - temp);
      	}

		if (difference/(n*n) < tolerance) {
			convergence = 1;
		}
		iterations ++;
	}


	MPI_Comm_rank(MPI_COMM_WORLD, &my_rank);
	if (convergence) {
		printf("Node %d: Gauss Seidel converged after %d iterations\n", my_rank, iterations);
	}
	else {
		printf("Node %d: Gauss Seidel not converged after %d iterations\n", my_rank, iterations);
	}
}

 // Generate random numbers
float generate_random(int max) {
	return ((float)rand() / (float)(RAND_MAX)) * max;
}

// maximum How many rows given to each node
int find_rows_with_max_values(int num_nodes, int n) {
	return (int)(ceil((n-2) / num_nodes) + 2);
}

// Find the position from which elements are gonna sent or received
int find_node_position(int node_index, int n, int max_num_rows) {
	return node_index * n * (max_num_rows-2);
}

//  how many elements are going to a given node
int find_node_elements(int node_index, int n, int max_num_rows) {

	int node_start_index = find_node_position(node_index, n, max_num_rows);
	int node_elements = max_num_rows * n;

	// Case1 : which the node receive the full set of elements
	if (node_start_index + node_elements <= (n*n)) {
		return node_elements;
	}

	// Case2: of the last node, which could get less elements
	else {
		return (n*n) - node_start_index;
	}
}

// Initialize 2D matrix in the master node
void initialize_master_matrix(float **A, int n, int m){

	*A = (float *) malloc(n * m * sizeof(float));
	for (int i = 0; i < (n*m); i++) {
		(*A)[i] = generate_random(MAX);
	}
}




// Initialize 2D matrix in the slaves nodes
void initialize_node_matrix(float **A, int number_elements) {
	*A = (float *) malloc(number_elements * sizeof(float));
}

