CC = mpicc
FLAGS = -c -O3
LINK_FLAGS = #-o3
OBJS = main.o mpi_comm.o sor.o utils.o
EXECUTABLE_NAME = parallel-sor


# Compile
all: $(OBJS)
	$(CC) -O3 $(OBJS) $(LINK_FLAGS) -o $(EXECUTABLE_NAME) -lm
	rm $(OBJS)

main.o: main.c
	$(CC) $(FLAGS) $(LINK_FLAGS) main.c

mpi_comm.o: mpi_comm.c
	$(CC) $(FLAGS) $(LINK_FLAGS) mpi_comm.c

sor.o: sor.c
	$(CC) $(FLAGS) $(LINK_FLAGS) sor.c

utils.o: utils.c
	$(CC) $(FLAGS) $(LINK_FLAGS) utils.c


# Clean-up
clean:
	rm -f $(EXECUTABLE_NAME) $(OBJS)
