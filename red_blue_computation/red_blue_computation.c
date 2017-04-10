#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

int main(int argc, char **argv) {
	int **grid;	// two-dimension grid
	int *grid_flat;	// one-dimension version of grid
	int n, t, c, max_iters;	// n-cell grid size, t-tile grid size, c-terminating threshold, max_iters- maximum number of iterations
	int n_itrs = 0;	// the iteration times
	int redcount = 0, bluecount = 0; // the count of red and blue
	int myid;
	int numprocs;
	MPI_Status status;
	time_t s;
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

	if (myid == 0) {
		if ((n % t != 0) || (numprocs > t)) {
			printf("Illegal arguments.\n");
			printf("Please enter the command in the following format:\n");
			printf("mpirun -np [proc num] red_blue_computation [cell grid size] [tile grid size] [terminating threshold] [maximum number of iterations]\n");
			printf("Note: tile grid size should divides cell grid size; process number should smaller than tile grid size.\n");
			return 0;
		}
	}

	MPI_Finalize();
	return 0;
}

