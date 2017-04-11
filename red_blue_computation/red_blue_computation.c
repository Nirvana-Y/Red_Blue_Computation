#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

// declera the functions
void allocate_memory(int l, int w, int **grid_flat, int ***grid);
void init_grid(int n, int ***grid);
void print_grid(int n, int ***grid);
int* distribute_row_for_processes(int numprocs, int n, int t);


int main(int argc, char **argv) {
	int *grid_flat;	// one-dimension version of grid
	int **grid;	// two-dimension grid
	int n, t, c, max_iters;	// n-cell grid size, t-tile grid size, c-terminating threshold, max_iters- maximum number of iterations
	int n_itrs = 0;	// the iteration times
	float threshold; // the threshold shown in percentage
	int redcount = 0, bluecount = 0; // the count of red and blue
	int myid;
	int numprocs;
	MPI_Status status;
	int i, j;

	MPI_Init(&argc, &argv);
	MPI_Comm_size(MPI_COMM_WORLD, &numprocs);
	MPI_Comm_rank(MPI_COMM_WORLD, &myid);

	// check the command line arguments
	if (myid == 0) {
		if (argc != 5) {
			printf("Wrong number of arguments.\n");
			printf("Please enter the command in the following format:\n");
			printf("mpirun -np [proc num] red_blue_computation [cell grid size] [tile grid size] [terminating threshold] [maximum number of iterations]\n");
			printf("Note: tile grid size should divides cell grid size; process number should smaller than tile grid size.\n");
			return 0;
		}
	}

	// parse the command line arguments
	n = atoi(argv[1]);
	t = atoi(argv[2]);
	c = atoi(argv[3]);
	max_iters = atoi(argv[4]);
	threshold = c / 100.0;

	if (myid == 0) {
		if ((n % t != 0) || ((numprocs - 1) > (n / t))) {
			printf("Illegal arguments.\n");
			printf("Please enter the command in the following format:\n");
			printf("mpirun -np [proc num] red_blue_computation [cell grid size] [tile grid size] [terminating threshold] [maximum number of iterations]\n");
			printf("Note: tile grid size should divides cell grid size; process number should smaller than tile row number.\n");
			return 0;
		}
	}

	int *row = distribute_row_for_processes(numprocs, n, t);

	if (myid == 0) {
		allocate_memory(n, n, &grid_flat, &grid);
		init_grid(n, &grid);
		print_grid(n, &grid);
	
		if (numprocs == 1) {
			return 0;
		}
		else {	
			for (i = 1; i < numprocs; i++) {
				int index = 0;
				for (j = 0; j < i; j++) {
					index = index + row[j];
				}
				index = index - row[i - 1];

				MPI_Send(&grid_flat[index * n], row[i - 1] * n, MPI_INT, i, 1, MPI_COMM_WORLD);
			}
			for (i = 1; i < numprocs; i++) {
				int index = 0;
				for (j = 0; j < i; j++) {
					index = index + row[j];
				}
				index = index - row[i - 1];
				
				MPI_Recv(&grid_flat[index * n], row[i - 1] * n, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			}
			//print_grid(n, &grid);
		}
	}
	else {
		allocate_memory(n, row[myid - 1] + 2, &grid_flat, &grid);
		MPI_Recv(&grid_flat[n], row[myid - 1] * n, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		
		// create ghost line for the easier manipulation of the grid
		MPI_Sendrecv(&grid_flat[(row[myid - 1] - 1) * n + n], n, MPI_INT, myid % (numprocs - 1) + 1, 1, &grid_flat[0], n, MPI_INT, (myid - 2 + (numprocs - 1)) % (numprocs - 1) + 1, 1, MPI_COMM_WORLD, &status);
		MPI_Sendrecv(&grid_flat[n], n, MPI_INT, (myid - 2 + (numprocs - 1)) % (numprocs - 1) + 1, 2, &grid_flat[row[myid - 1] * n + n], n, MPI_INT, myid % (numprocs - 1) + 1, 2, MPI_COMM_WORLD, &status);
		
		printf("\n");
		for (i = 0; i < row[myid - 1] + 2; i++) {
			for (j = 0; j < n;j++) {
				printf("%d", grid[i][j]);
			}
			printf("\n");
		}

		MPI_Send(&grid_flat[n], row[myid - 1] * n, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

	MPI_Finalize();
	return 0;
}

// memory allocation for grid
void allocate_memory(int l, int w, int **grid_flat, int ***grid) {
	int count = l * w;
	*grid_flat = (int *)malloc(sizeof(int) * count);
	*grid = (int **)malloc(sizeof(int *) * w);
	int i;

	for (i = 0; i < w; i++) {
		(*grid)[i] = &((*grid_flat)[i * l]);
	}	
}

// initialize grid
void init_grid(int n, int ***grid) {
	time_t s;
	srand((unsigned)time(&s));
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			(*grid)[i][j] = rand() % 3;
		}
	}
}

// print grid
void print_grid(int n, int ***grid) {
	int i, j;

	for (i = 0; i < n; i++) {
		for (j = 0; j < n; j++) {
			printf("%d", (*grid)[i][j]);
		}
		printf("\n");
	}
	printf("\n");
}

// calculate the row distributed to each process
int* distribute_row_for_processes(int numprocs, int n, int t) {
	int *row = (int *)malloc(sizeof(int) * (numprocs - 1)); // array 'row' store the number of row distributed to each process
	int tile_row_count = n / t;
	int base = tile_row_count / (numprocs - 1);
	int remain = tile_row_count % (numprocs - 1);
	int i;

	// calculate the tile row distributed to each process
	for (i = 0; i < (numprocs - 1); i++) {
		row[i] = base;
	}
	for (i = 0; i < remain; i++) {
		row[i] = row[i] + 1;
	}

	// calculate the row distributed to each process
	for (i = 0; i < (numprocs - 1); i++) {
		row[i] = row[i] * t;
	}

	return row;
}