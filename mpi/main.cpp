#include <mpi.h>
#include <iostream>
#include "../utils/utils.h"

#define EQUATION_SOLVER = 1
#define SYSTEM_SOLVER = 2


#define RANK_ROOT 0

struct Task {
    int type = 0;
};

struct Equation_solver_Task : Task {
    float *A_row;
    float *x;
    float b;

    int N;
    int index;
};

struct System_solver_task : Task {
    float **A;
    float *x;
    float *b;
    int N;

    System_solver_task(float **A, float *x, float *b, int N) 
        : A(A), x(x), b(b), N(N) {
    };
};

enum Channels {
    REQUEST_WORK = 1,
    ADD_TASK_TO_POOL = 2,
    SEND_SYSTEM_SOLVED = 3,
    SEND_EQUATION_SOLVED = 4,
    CHECK_RESULTS = 5
};

void process_run() {

};

void handle_add_task_to_pool() {

};

void handle_send_work() {

}

void handle_got_system_result() {

}

void handle_got_equation_result() {
    
}

void handle_check_if_iteration_done() {
    
}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank != RANK_ROOT) {
        process_run();
        MPI_Finalize();
        return 0;
    } 

    input_gen input_generator(1, input_type::VERY_SMALL);
    auto systems = input_generator.get_input();
    int systems_count = input_generator.no_systems;

    auto tasks = std::vector<Task*>();

    for (int i = 0; i < systems_count; i++) {
        Task* task = new System_solver_task(systems[i]->A, systems[i]->x, systems[i]->b, systems[i]->N);
        tasks.push_back(task);
    }

    MPI_Request get_work_request = MPI_REQUEST_NULL;
    MPI_Request add_task_to_pool = MPI_REQUEST_NULL;
    MPI_Request sent_system_result = MPI_REQUEST_NULL;
    MPI_Request sent_equation_result = MPI_REQUEST_NULL;
    MPI_Request check_result_of_iteration = MPI_REQUEST_NULL;

    int dummy;
    int recv_task_equation_index;
    int system_solved_index;
    int equation_solved_index;
    int operation_completed = 0;
    MPI_Status status;
    while (true) {
        if (get_work_request == MPI_REQUEST_NULL)
            MPI_Irecv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_WORK, MPI_COMM_WORLD, &get_work_request);

        if (add_task_to_pool == MPI_REQUEST_NULL)
            MPI_Irecv(&recv_task_equation_index, 1, MPI_INT, MPI_ANY_SOURCE, ADD_TASK_TO_POOL, MPI_COMM_WORLD, &add_task_to_pool);

        if (sent_system_result == MPI_REQUEST_NULL)
            MPI_Irecv(&system_solved_index, 1, MPI_INT, MPI_ANY_SOURCE, SEND_SYSTEM_SOLVED, MPI_COMM_WORLD, &sent_system_result);

        if (sent_equation_result == MPI_REQUEST_NULL)
            MPI_Irecv(&equation_solved_index, 1, MPI_INT, MPI_ANY_SOURCE, SEND_SYSTEM_SOLVED, MPI_COMM_WORLD, &sent_equation_result);


        MPI_Test(&get_work_request, &operation_completed, &status);
        if (operation_completed)
            handle_send_work();
        
        MPI_Test(&add_task_to_pool, &operation_completed, &status);
        if (operation_completed)
            handle_add_task_to_pool();

        MPI_Test(&sent_system_result, &operation_completed, &status);
        if (operation_completed)
            handle_got_system_result();
        
        MPI_Test(&sent_equation_result, &operation_completed, &status);
        if (operation_completed)
            handle_got_equation_result();
        
        MPI_Test(&check_result_of_iteration, &operation_completed, &status);
        if (operation_completed)
            handle_check_if_iteration_done();

    }

    MPI_Finalize();
    return 0;
}