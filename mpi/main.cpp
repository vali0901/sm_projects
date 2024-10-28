#include <mpi.h>
#include <iostream>

int NUMBERS = 10;

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    int numbers_per_proc = NUMBERS / world_size;
    int start = world_rank * numbers_per_proc + 1;
    int end = start + numbers_per_proc;

    int local_sum = 0;
    for (int i = start; i < end; ++i) {
        local_sum += i;
    }

    int global_sum = 0;
    MPI_Reduce(&local_sum, &global_sum, 1, MPI_INT, MPI_SUM, 0, MPI_COMM_WORLD);

    if (world_rank == 0) {
        std::cout << "Sum of the first 10 numbers is " << global_sum << std::endl;
    }

    MPI_Finalize();
    return 0;
}