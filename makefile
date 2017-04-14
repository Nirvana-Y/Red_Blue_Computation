red_blue_computation: red_blue_computation.c
	mpicc -o $@ $^

clean:
	rm red_blue_computation

