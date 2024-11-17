#include <mpi.h>
#include <iostream>
#include "../utils/utils.h"

#define EQUATION_SOLVER = 1
#define SYSTEM_SOLVER = 2

bool introduced_all_tasks = false;

bool operations_finished() {
    return !introduced_all_tasks;
}

#define RANK_ROOT 0
struct Task {
    
};

struct Equation_solver_task : Task {
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

auto equation_tasks = std::queue<Equation_solver_task*>();
auto system_tasks = std::queue<System_solver_task*>();
int world_size;

enum Channels {
    REQUEST_WORK = 1,
    ADD_TASK_TO_POOL = 2,
    SEND_SYSTEM_SOLVED = 3,
    SEND_EQUATION_SOLVED = 4,
    CHECK_RESULTS = 5,
    SEND_EQUATION_TASK = 6,
    SEND_SYSTEM_TASK = 7
};

Equation_solver_task read_equation_task(int N) {

    Equation_solver_task task;
    task.N = N;
    task.A_row = new float[N];
    task.b = 0;
    task.x = new float[N];
    task.index = 0;

    MPI_Recv(&(task.index), 1, MPI_INT, RANK_ROOT, SEND_EQUATION_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&(task.A_row), N, MPI_FLOAT, RANK_ROOT, SEND_EQUATION_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&(task.x), N, MPI_FLOAT, RANK_ROOT, SEND_EQUATION_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&(task.b), 1, MPI_FLOAT, RANK_ROOT, SEND_EQUATION_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    return task;
}

System_solver_task read_system_solver_task(int N) {
    float **A = new float*[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new float[N];
    }
    float *b = new float[N];
    float *x = new float[N];
    
    System_solver_task task(A, x, b, N);

    for (int i = 0; i < N; i++) {
        MPI_Recv(task.A[i], N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Recv(task.x, N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(task.b, N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return task;
};

float solve_equation_task(const Equation_solver_task &task) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    std::cout << "[" << world_rank << " worker] " << "solved equation\n";
};

void solve_system_task(const System_solver_task &task) {
    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);
    std::cout << "[" << world_rank << " worker] " << "solved system\n";
};

void send_equation_result(float result) {
    MPI_Send(&result, 1, MPI_FLOAT, MPI_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);
}

void worker_process_run() {
    int N;
    int dummy;
    MPI_Status status;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    while (true) {
        MPI_Send(&dummy, 1, MPI_INT, RANK_ROOT, REQUEST_WORK, MPI_COMM_WORLD);

        MPI_Recv(&N, 1, MPI_INT, RANK_ROOT, MPI_ANY_TAG, MPI_COMM_WORLD, &status);

        if (N == -1) {
            return;
        }
        int tag = status.MPI_TAG;
        if (tag == SEND_EQUATION_TASK) {
            auto task = read_equation_task(N);
            auto result = solve_equation_task(task);
            send_equation_result(result);
            std::cout << "[" << rank << " worker] has sent equation result" << std::endl;

            // TODO: will move the frees to the corresponding destructor of task
            free(task.A_row);
            free(task.x);
        } else {
            auto task = read_system_solver_task(N);
            solve_system_task(task);

            // TODO: will move the frees to the corresponding destructor of task
            for (int i = 0; i < task.N; i++)
                free(task.A[i]);
            free(task.A);
            free(task.b);
            free(task.x);
        }
    }
};

void handle_add_task_to_pool(MPI_Status status) {

};

void send_equation_task(int destination, Equation_solver_task *task, Channels ch) {
    MPI_Send(&(task->N), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(&(task->index), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(task->A_row, task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(task->x, task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(&(task->b), 1, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);
}

void send_system_task(int destination, System_solver_task *task, Channels ch) {
    MPI_Send(&(task->N), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    for (int i = 0; i < task->N; i++)
        MPI_Send(&(task->A[i]), task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(&(task->x), task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);
    MPI_Send(&(task->b), task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

}

bool handle_send_work(MPI_Status status) {
    int destination = status.MPI_SOURCE;

    // we prioritize the tasks that solve equations - solves the problem of starvation    
    if (!equation_tasks.empty()) {
        send_equation_task(destination, equation_tasks.front(), SEND_EQUATION_TASK);
        equation_tasks.pop();
        return true;
    }
    
    if (!system_tasks.empty()) {
        send_system_task(destination, system_tasks.front(), SEND_SYSTEM_TASK);
        system_tasks.pop();
        return true;
    }

    // TODO CHECK IF THERE IS NO WORK IN PROGRESS
    // TODO if work available return false
    if (!operations_finished()) {
        std::cout << "not all ops finished\n";
        return false;
    }

    int close = -1;
    MPI_Send(&close, 1, MPI_INT, destination, SEND_EQUATION_TASK, MPI_COMM_WORLD);
    return true;
}

void handle_got_system_result(MPI_Status status) {

}

void handle_got_equation_result(MPI_Status status) {
    
}

void handle_check_if_iteration_done(MPI_Status status) {

}


int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank != RANK_ROOT) {
        worker_process_run();
        std::cout << "[" << world_rank << " worker] finished work" << std::endl;

        MPI_Finalize();
        return 0;
    } 

    input_gen input_generator(2, input_type::SMALL_SMALL);
    auto systems = input_generator.get_input();
    int systems_count = input_generator.no_systems;
    std::cout << "before\n";

    for (int i = 0; i < systems_count; i++) {
        std::cout << systems[i] << std::endl;
        // auto task = new System_solver_task(systems[i]->A, systems[i]->x, systems[i]->b, systems[i]->N);
        // system_tasks.push(task);
    }
    std::cout << "after\n";

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
        if (operation_completed) {
            if (handle_send_work(status))
                get_work_request = MPI_REQUEST_NULL;
        }
        
        MPI_Test(&add_task_to_pool, &operation_completed, &status);
        if (operation_completed) {
            handle_add_task_to_pool(status);
            add_task_to_pool = MPI_REQUEST_NULL;
        }

        MPI_Test(&sent_system_result, &operation_completed, &status);
        if (operation_completed) {
            handle_got_system_result(status);
            sent_system_result = MPI_REQUEST_NULL;
        }
        
        MPI_Test(&sent_equation_result, &operation_completed, &status);
        if (operation_completed) {
            handle_got_equation_result(status);
            sent_equation_result = MPI_REQUEST_NULL;
            
        }
        
        MPI_Test(&check_result_of_iteration, &operation_completed, &status);
        if (operation_completed) {
            handle_check_if_iteration_done(status);
            check_result_of_iteration = MPI_REQUEST_NULL;
        }

    }
    std::cout << "[ROOT] finished work" << std::endl;
    MPI_Finalize();
    return 0;
}