#include <mpi.h>
#include <iostream>

#include "utils.h"
#include "messages.h"
#include "task_pool.h"
#include "worker.h"


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if(world_rank == MASTER) {
        struct input_gen input(0, input_type::VERY_SMALL);
        struct system **systems = input.get_input();
        TaskPool task_pool(input.no_systems, world_size - 1);
        task_pool.load_jacobi_tasks(systems);
        task_pool.loop();

        #ifdef DEBUG
        log_solution_correctness(systems, input.no_systems, "hybrid mpi openmp");
        #endif

        for(int i = 0; i < input.no_systems; i++)
            delete systems[i];
        delete[] systems;
    } else {
        Worker worker(world_rank);
        printf("Worker %d started\n", world_rank);
        worker.loop();
    }

    MPI_Finalize();
    return 0;
}