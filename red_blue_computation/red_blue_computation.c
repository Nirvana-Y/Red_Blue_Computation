#include <mpi.h>
#include <stdlib.h>
#include <stdio.h>
#include <time.h>

#define RED 1
#define BLUE 2
#define BOTH 3
#define PARALLEL 1
#define SEQUENTIAL 2

// declera the functions
void allocate_memory(int l, int w, int **grid_flat, int ***grid);
void init_grid(int n, int ***grid);
void print_grid(int n, int ***grid);
int* distribute_row_for_processes(int numprocs, int n, int t);
void do_red(int l, int w, int ***grid);
void do_blue(int l, int w, int ***grid);
int red_blue_computation(float **red_blue_array, int ***grid, int tile_number, int n, int t, float threshold, int type);
void print_computation_result(float **red_blue_array, int tile_number);
void sequential_computation(int **grid_flat, int ***grid, int tile_number, int n, int t, float threshold, int max_iters);



int main(int argc, char **argv) {
	int *grid_flat;	// one-dimension version of grid
	int **grid;	// two-dimension grid
	int *grid_flat_copy; // the copy of one-dimension version of grid
	int **grid_copy; // the copy of the initial two-dimension grid
	int n, t, c, max_iters;	// n-cell grid size, t-tile grid size, c-terminating threshold, max_iters- maximum number of iterations
	int n_itrs = 0;	// the iteration times
	int finished_flag = 0; // the flag show whether the iteration finished or not
	int finished_flag_p = 0; //	the finished flag for processes
	float threshold; // the threshold shown in percentage
	int tile_number; // the number of tile in the grid
	int redcount = 0, bluecount = 0; // the count of red and blue
	int type; // the algorithm type, parallel or sequential
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

	int *row;
	if (numprocs != 1) {
		row = distribute_row_for_processes(numprocs, n, t);
	}

	if (myid == 0) {
		tile_number = (n * n) / (t * t);
		float *red_blue_array = (float *)malloc(sizeof(float) * tile_number * 3);	// store the red and blue ratio in grid

		allocate_memory(n, n, &grid_flat, &grid);
		allocate_memory(n, n, &grid_flat_copy, &grid_copy);

		init_grid(n, &grid);
		memcpy(grid_flat_copy, grid_flat, sizeof(int) * n * n);

		printf("The initial grid: \n");
		print_grid(n, &grid);	
	
		if (numprocs == 1) {
			sequential_computation(&grid_flat, &grid, tile_number, n, t, threshold, max_iters);
			goto EXIT;
		}
		else {	
			// send the sub-grid to corresponding processes
			for (i = 1; i < numprocs; i++) {
				int index = 0;
				for (j = 0; j < i; j++) {
					index = index + row[j];
				}
				index = index - row[i - 1];

				MPI_Send(&grid_flat[index * n], row[i - 1] * n, MPI_INT, i, 1, MPI_COMM_WORLD);
			}
			
			// terminate when the computation in other processes terminate
			while (!finished_flag && n_itrs < max_iters) {
				// receive the iteration number from process 1
				MPI_Recv(&n_itrs, 1, MPI_INT, 1, 1, MPI_COMM_WORLD, &status);
				MPI_Allreduce(&finished_flag_p, &finished_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
			}

			// receive the final grid and the computation result from other processes
			for (i = 1; i < numprocs; i++) {
				int index = 0;
				for (j = 0; j < i; j++) {
					index = index + row[j];
				}
				index = index - row[i - 1];

				MPI_Recv(&red_blue_array[index * n / (t * t) * 3], (n * row[i - 1]) / (t * t) * 3, MPI_FLOAT, i, 1, MPI_COMM_WORLD, &status);
				MPI_Recv(&grid_flat[index * n], row[i - 1] * n, MPI_INT, i, 1, MPI_COMM_WORLD, &status);
			}
			
			printf("The parallel computation result: \n");
			printf("After %d interations, the final grid: \n", n_itrs);
			print_grid(n, &grid);
			print_computation_result(&red_blue_array, tile_number);

			sequential_computation(&grid_flat_copy, &grid_copy, tile_number, n, t, threshold, max_iters);
		}
	}
	else {
		tile_number = (n * row[myid - 1]) / (t * t);
		type = PARALLEL;
		float *red_blue_array = (float *)malloc(sizeof(float) * tile_number * 3);	// store the red and blue ratio in each tile grid

		allocate_memory(n, row[myid - 1] + 2, &grid_flat, &grid);
		MPI_Recv(&grid_flat[n], row[myid - 1] * n, MPI_INT, 0, 1, MPI_COMM_WORLD, &status);
		while (!finished_flag && n_itrs < max_iters) {
			n_itrs = n_itrs + 1; // renew the iteration number

			// send the iteration number to process 0
			if (myid == 1) {
				MPI_Send(&n_itrs, 1, MPI_INT, 0, 1, MPI_COMM_WORLD);
			}

			// create ghost line for the easier manipulation of the grid
			MPI_Sendrecv(&grid_flat[(row[myid - 1] - 1) * n + n], n, MPI_INT, myid % (numprocs - 1) + 1, 1, &grid_flat[0], n, MPI_INT, (myid - 2 + (numprocs - 1)) % (numprocs - 1) + 1, 1, MPI_COMM_WORLD, &status);
			MPI_Sendrecv(&grid_flat[n], n, MPI_INT, (myid - 2 + (numprocs - 1)) % (numprocs - 1) + 1, 2, &grid_flat[row[myid - 1] * n + n], n, MPI_INT, myid % (numprocs - 1) + 1, 2, MPI_COMM_WORLD, &status);

			do_red(n, row[myid - 1] + 2, &grid);
			do_blue(n, row[myid - 1] + 2, &grid);

			finished_flag_p = red_blue_computation(&red_blue_array, &grid, tile_number, n, t, threshold, type);
			MPI_Allreduce(&finished_flag_p, &finished_flag, 1, MPI_INT, MPI_SUM, MPI_COMM_WORLD);
		}
		MPI_Send(&red_blue_array[0], tile_number * 3, MPI_FLOAT, 0, 1, MPI_COMM_WORLD);
		MPI_Send(&grid_flat[n], row[myid - 1] * n, MPI_INT, 0, 1, MPI_COMM_WORLD);
	}

EXIT:
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

// red color movement
void do_red(int l, int w, int ***grid) {
	int i, j;
	for (i = 0; i < w; i++) {
		//the first column.
		if ((*grid)[i][0] == 1 && (*grid)[i][1] == 0) {
			(*grid)[i][0] = 4;
			(*grid)[i][1] = 3;
		}
		//the rest column.
		for (j = 1; j < l; j++) {
			if ((*grid)[i][j] == 1 && ((*grid)[i][(j + 1) % l] == 0)) {
				(*grid)[i][j] = 0;
				(*grid)[i][(j + 1) % l] = 3;
			}
			else if ((*grid)[i][j] == 3)
				(*grid)[i][j] = 1;
		}
		//cast back to the changed colours.
		if ((*grid)[i][0] == 3)
			(*grid)[i][0] = 1;
		else if ((*grid)[i][0] == 4)
			(*grid)[i][0] = 0;
	}
}

// blue color movement
void do_blue(int l, int w, int ***grid) {
	int i, j;
	for (j = 0; j < l; j++) {
		if ((*grid)[0][j] == 2 && (*grid)[1][j] == 0) {
			(*grid)[0][j] = 4;
			(*grid)[1][j] = 3;
		}
		for (i = 1; i < w; i++) {
			if ((*grid)[i][j] == 2 && (*grid)[(i + 1) % w][j] == 0) {
				(*grid)[i][j] = 0;
				(*grid)[(i + 1) % w][j] = 3;
			}
			else if ((*grid)[i][j] == 3)
				(*grid)[i][j] = 2;
		}
		if ((*grid)[0][j] == 3)
			(*grid)[0][j] = 2;
		else if ((*grid)[0][j] == 4)
			(*grid)[0][j] = 0;
	}
}

// check whether there is any tile exceeding the threshold
int red_blue_computation(float **red_blue_array, int ***grid, int tile_number, int n, int t, float threshold, int type) {
	float tile_count = t * t;	// the number of cells in each tile
	int tile_row, tile_column; // the index of tile in tile grid
	int redcount = 0, bluecount = 0;
	float red_ratio, blue_ratio;
	int finished_flag = 0;
	int i, j, k;

	for (i = 0; i < tile_number; i++) {
		tile_row = i / (n / t);
		tile_column = i % (n / t);
		for (j = t * tile_row; j < t * tile_row + t; j++) {
			for (k = t * tile_column; k < t * tile_column + t; k++) {
				if (type == PARALLEL) {
					if ((*grid)[j + 1][k] == 1) {
						redcount = redcount + 1;
					}
					if ((*grid)[j + 1][k] == 2) {
						bluecount = bluecount + 1;
					}
				}
				else {
					if ((*grid)[j][k] == 1) {
						redcount = redcount + 1;
					}
					if ((*grid)[j][k] == 2) {
						bluecount = bluecount + 1;
					}
				}
			}
		}

		red_ratio = redcount / tile_count;
		blue_ratio = bluecount / tile_count;
		red_ratio = (int)(100.0 * red_ratio + 0.5) / 100.0;
		blue_ratio = (int)(100.0 * blue_ratio + 0.5) / 100.0;
		(*red_blue_array)[3 * i + 1] = red_ratio;
		(*red_blue_array)[3 * i + 2] = blue_ratio;

		if (red_ratio > threshold) {
			(*red_blue_array)[3 * i] = RED;	
			finished_flag = 1;
		}

		if (blue_ratio > threshold) {
			(*red_blue_array)[3 * i] = BLUE;
			finished_flag = 1;
		}

		if (blue_ratio > threshold && red_ratio > threshold) {
			(*red_blue_array)[3 * i] = BOTH;
			finished_flag = 1;
		}		

		redcount = 0; 
		bluecount = 0;
	}
	return finished_flag;
}

void print_computation_result(float **red_blue_array, int tile_number) {
	int exceed_tile = 0;
	int i;

	for (i = 0; i < tile_number; i++) {
		if ((*red_blue_array)[3 * i] == RED) {
			printf("In tile %d, the red color exceed the threshold with the ratio %.2f.\n", i, (*red_blue_array)[3 * i + 1]);
			exceed_tile = exceed_tile + 1;
		}

		if ((*red_blue_array)[3 * i] == BLUE) {
			printf("In tile %d, the blue color exceed the threshold with the ratio %.2f.\n", i, (*red_blue_array)[3 * i + 2]);
			exceed_tile = exceed_tile + 1;
		}

		if ((*red_blue_array)[3 * i] == BOTH) {
			printf("In tile %d, the red color exceed the threshold with the ratio %.2f and the blue color exceed the threshold with the ratio %.2f.\n", i, (*red_blue_array)[3 * i + 1], (*red_blue_array)[3 * i + 2]);
			exceed_tile = exceed_tile + 1;
		}
	}

	if (exceed_tile == 0) {
		printf("There is no tile containning color exceeding threshold.\n");
		printf("The computation terminated becausu the maximum iteration number has been reached.\n");
	}
	printf("\n");
}

void sequential_computation(int **grid_flat, int ***grid, int tile_number, int n, int t, float threshold, int max_iters) {
	int finished_flag = 0;
	int n_itrs = 0;
	int type = SEQUENTIAL;
	float *red_blue_array = (float *)malloc(sizeof(float) * tile_number * 3);	// store the red and blue ratio in grid

	while (!finished_flag && n_itrs < max_iters) {
		n_itrs = n_itrs + 1; // renew the iteration number

		do_red(n, n, grid);
		do_blue(n, n, grid);

		finished_flag = red_blue_computation(&red_blue_array, grid, tile_number, n, t, threshold, type);
	}

	printf("The sequential computation result: \n");
	printf("After %d interations, the final grid: \n", n_itrs);
	print_grid(n, grid);
	print_computation_result(&red_blue_array, tile_number);
}