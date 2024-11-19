#include <mpi.h>
#include <iostream>
#include "../utils/utils.h"

#define EQUATION_SOLVER = 1
#define SYSTEM_SOLVER = 2

#define RANK_ROOT 0
#define EPS 0.005

struct Equation_solver_task {
    float *A_row;
    float *x;
    float b;

    int N;
    int x_index;
    int system_ID;

    int dealloc_responsability = false;

    Equation_solver_task(int N) {
        this->N = N;
        A_row = new float[N];
        b = 0;
        x = new float[N];
        x_index = 0;
        system_ID = 0;
        dealloc_responsability = true;
    }

    Equation_solver_task(int N, int x_index, float b, float *x, float *A_row, int system_ID) {
        this->N = N;
        this->x_index = x_index;
        this->b = b;
        this->x = x;
        this->A_row = A_row;
        this->system_ID = system_ID;
    }

    ~Equation_solver_task() {
        if (dealloc_responsability) {
            delete[] A_row;
            delete[] x;
        }
    }
};

struct System_solver_task {
    float **A;
    float *x;
    float *b;
    int N;
    int ID;

    bool solved = false;
    int equations_solved = 0;
    int assigned_process_rank = -1;
    
    bool dealloc_responsability = false;

    System_solver_task(float **A, float *x, float *b, int N, int ID) 
        : A(A), x(x), b(b), N(N), ID(ID) {
    };

    System_solver_task(int N)
    {
        this->N = N;
        A = new float*[N];
        for (int i = 0; i < N; ++i) {
            A[i] = new float[N];
        }
        b = new float[N];
        x = new float[N];
        ID = -1;
        dealloc_responsability = true;
    }

    ~System_solver_task() 
    {
        if (dealloc_responsability) {
            for (int i = 0; i < N; i++)
                delete[] A[i];

            delete[] A;
            delete[] b;
            delete[] x;
        }
    }
};

auto equation_tasks = std::queue<Equation_solver_task*>();
auto system_tasks = std::queue<System_solver_task*>();
auto system_tasks_registry = std::vector<System_solver_task*>();

long unsigned int count_systems_solved = 0;

enum Channels {
    REQUEST_WORK = 1,
    ADD_TASK_TO_POOL = 2,
    SEND_SYSTEM_SOLVED = 3,
    SEND_EQUATION_SOLVED = 4,
    GET_SYSTEM_NEW_RESULTS = 5,

    SEND_EQUATION_TASK = 6,
    SEND_SYSTEM_TASK = 7
};

enum Get_work_return {
    CLOSE = -1,
    RETRY = -2,
    SUCCESS = 0
};

bool operations_finished() {
    return count_systems_solved == system_tasks_registry.size();
}

Equation_solver_task* read_equation_task(int N, int source, int ch) {

    Equation_solver_task *task = new Equation_solver_task(N);

    MPI_Recv(&(task->system_ID), 1, MPI_INT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&(task->x_index), 1, MPI_INT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(task->A_row, N, MPI_FLOAT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(task->x, N, MPI_FLOAT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&(task->b), 1, MPI_FLOAT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return task;
}

System_solver_task* read_system_task(int N) {
    System_solver_task *task = new System_solver_task(N);
    
    MPI_Recv(&(task->ID), 1, MPI_INT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < N; i++) {
        MPI_Recv(task->A[i], N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Recv(task->x, N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(task->b, N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return task;
};

float solve_equation_task(const Equation_solver_task *task) {
    float result = task->b;

    for (int j = 0; j < task->N; ++j) {
        if (task->x_index == j)
            continue;

        result -= task->A_row[j] * task->x[j];
    }
    result /= task->A_row[task->x_index];

    return result;
};

void get_system_latest_iteration(int N, float* result) {
    MPI_Recv(result, N, MPI_FLOAT, RANK_ROOT, GET_SYSTEM_NEW_RESULTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void send_equation_task(int destination, Equation_solver_task *task, Channels ch);

void solve_system_task(const System_solver_task &task) {
    int max_iter = task.N * 100;
    float *x_new = new float[task.N];
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);

    for (int iter = 0; iter < max_iter; ++iter) {
        
        double error = error_MSE(task.x, x_new, task.N);
        for (int i = 0; i < task.N; ++i) {
            auto equation_task = new Equation_solver_task(task.N, i, task.b[i], task.x, task.A[i], task.ID);

            send_equation_task(RANK_ROOT, equation_task, ADD_TASK_TO_POOL);
            delete equation_task;
        }
        get_system_latest_iteration(task.N, x_new);

        memcpy(task.x, x_new, task.N * sizeof(float));

        if (error < EPS) {
            break;
        }
    }
    delete[] x_new;
};

void send_equation_result(float result, int system_ID, int x_index) {
    int dummy = -1;
    MPI_Send(&dummy, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);

    MPI_Send(&system_ID, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);

    MPI_Send(&x_index, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);

    MPI_Send(&result, 1, MPI_FLOAT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);
};

void send_system_result(System_solver_task task) {
    MPI_Send(&(task.ID), 1, MPI_INT, RANK_ROOT, SEND_SYSTEM_SOLVED, MPI_COMM_WORLD);
};

void worker_process_run() {
    int N;
    int dummy;
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    Channels tag;
    Get_work_return op_result;

    while (true) {
        MPI_Send(&dummy, 1, MPI_INT, RANK_ROOT, REQUEST_WORK, MPI_COMM_WORLD);

        MPI_Recv(&op_result, 1, MPI_INT, RANK_ROOT, REQUEST_WORK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (op_result == Get_work_return::CLOSE) {
            return;
        }
        if (op_result == Get_work_return::RETRY) {
            continue;
        }
        MPI_Recv(&tag, 1, MPI_INT, RANK_ROOT, REQUEST_WORK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

        if (tag == SEND_EQUATION_TASK) {
            MPI_Recv(&N, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

            auto task = read_equation_task(N, RANK_ROOT, SEND_EQUATION_TASK);
            auto result = solve_equation_task(task);
            send_equation_result(result, task->system_ID, task->x_index);

            // TODO: will move the frees to the corresponding destructor of task
            delete task;
        } else {
            if (tag == SEND_SYSTEM_TASK) {
                MPI_Recv(&N, 1, MPI_INT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                auto task = read_system_task(N);
                solve_system_task(*task);
                send_system_result(*task);

            }
        }
    }
};

void handle_add_task_to_pool(MPI_Status status, int N) {
    // std::cout << "[ROOT] received a task to add in the pool from " << status.MPI_SOURCE << std::endl;
    auto task = read_equation_task(N, status.MPI_SOURCE, ADD_TASK_TO_POOL);
    // std::cout << "[ROOT] has added the task in the pool\n";
    equation_tasks.push(task);
};

void send_equation_task(int destination, Equation_solver_task *task, Channels ch) {
    //std::cout << "send_equation_task A = " << std::endl;
    for (int i = 0; i < task->N; i++) {
        //std::cout << task->A_row[i] << ' ';
    }
    //std::cout << '\n';
    MPI_Send(&(task->N), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(&(task->system_ID), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(&(task->x_index), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(task->A_row, task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(task->x, task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(&(task->b), 1, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);
}

void send_system_task(int destination, System_solver_task *task, Channels ch) {
    //std::cout << "send_system_task dest " << destination << std::endl;

    MPI_Send(&(task->N), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(&(task->ID), 1, MPI_INT, destination, ch, MPI_COMM_WORLD);

    for (int i = 0; i < task->N; i++)
        MPI_Send(task->A[i], task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

    MPI_Send(task->x, task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);
    MPI_Send(task->b, task->N, MPI_FLOAT, destination, ch, MPI_COMM_WORLD);

}

bool handle_send_work(MPI_Status status, int &count_workers_closed) {
    int destination = status.MPI_SOURCE;

    // we prioritize the tasks that solve equations - solves the problem of starvation    
    if (!equation_tasks.empty()) {
        Get_work_return result = Get_work_return::SUCCESS;
        MPI_Send(&result, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);

        Channels tag = SEND_EQUATION_TASK;
        MPI_Send(&tag, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);
        
        send_equation_task(destination, equation_tasks.front(), SEND_EQUATION_TASK);
        equation_tasks.pop();
        return true;
    }
    
    if (!system_tasks.empty()) {
        Get_work_return result = Get_work_return::SUCCESS;
        MPI_Send(&result, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);

        Channels tag = SEND_SYSTEM_TASK;
        MPI_Send(&tag, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);

        send_system_task(destination, system_tasks.front(), SEND_SYSTEM_TASK);
        system_tasks.front()->assigned_process_rank = destination;
        
        system_tasks.pop();
        return true;
    }

    if (!operations_finished())
        return false;

    int close = Get_work_return::CLOSE;
    MPI_Send(&close, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);
    count_workers_closed++;
    return true;
}

void send_updated_values_for_equation(int destination, float *data, int N) {
    MPI_Send(data, N, MPI_FLOAT, destination, GET_SYSTEM_NEW_RESULTS, MPI_COMM_WORLD);
}


void handle_got_equation_result(MPI_Status status) {
    int x_index;
    float result;
    int system_ID;

    MPI_Recv(&system_ID, 1, MPI_INT, status.MPI_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&x_index, 1, MPI_INT, status.MPI_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    MPI_Recv(&result, 1, MPI_FLOAT, status.MPI_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (const auto& task : system_tasks_registry) {
        if (task->ID == system_ID) {
            task->x[x_index] = result;
            task->equations_solved++;
            if (task->equations_solved == task->N) {
                task->equations_solved = 0;;
                send_updated_values_for_equation(task->assigned_process_rank, task->x, task->N);
            }

            return;
        }
    }

}

void handle_got_system_result(int system_ID) {
    for(const auto& task : system_tasks_registry) {
        if (task->ID == system_ID) {
            task->solved = true;
            count_systems_solved++;

            std::cout << "THE RESULT FOR " << system_ID << " IS \n";
            for (int i = 0; i < task->N; i++) {
                std::cout << task->x[i] << ' ';
            }
            std::cout << '\n';
            return;
        }
    }
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);;

    int world_size;
    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank != RANK_ROOT) {
        worker_process_run();
        MPI_Finalize();
        return 0;
    } 

    input_gen input_generator(2, input_type::SMALL_SMALL);
    auto systems = input_generator.get_input();
    int systems_count = input_generator.no_systems;

    for (int i = 0; i < systems_count; i++) {
        auto task = new System_solver_task(systems[i]->A, systems[i]->x, systems[i]->b, systems[i]->N, i);
        system_tasks.push(task);
        system_tasks_registry.push_back(task);
    }

    MPI_Request get_work_request = MPI_REQUEST_NULL;
    MPI_Request add_task_to_pool = MPI_REQUEST_NULL;
    MPI_Request worker_sent_system_result = MPI_REQUEST_NULL;
    MPI_Request worker_sent_equation_result = MPI_REQUEST_NULL;

    int dummy;
    int recv_equation_task_N;
    int system_solved_index;
    int equation_solved_dummy = -1;
    int operation_completed = 0;

    MPI_Status status;
    
    int count_workers_closed = 0;
    bool have_worker_ready = false;
    MPI_Status status_saved;

    while (count_workers_closed < (world_size - 1)) {
        if (have_worker_ready && operations_finished()) {
            handle_send_work(status_saved, count_workers_closed);
            have_worker_ready = false;
        }

        if (!have_worker_ready || operations_finished())
            if (get_work_request == MPI_REQUEST_NULL)
                MPI_Irecv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_WORK, MPI_COMM_WORLD, &get_work_request);

        if (add_task_to_pool == MPI_REQUEST_NULL)
            MPI_Irecv(&recv_equation_task_N, 1, MPI_INT, MPI_ANY_SOURCE, ADD_TASK_TO_POOL, MPI_COMM_WORLD, &add_task_to_pool);

        if (worker_sent_system_result == MPI_REQUEST_NULL)
            MPI_Irecv(&system_solved_index, 1, MPI_INT, MPI_ANY_SOURCE, SEND_SYSTEM_SOLVED, MPI_COMM_WORLD, &worker_sent_system_result);

        if (worker_sent_equation_result == MPI_REQUEST_NULL)
            MPI_Irecv(&equation_solved_dummy, 1, MPI_INT, MPI_ANY_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, &worker_sent_equation_result);


        if (!have_worker_ready || operations_finished()) {
            MPI_Test(&get_work_request, &operation_completed, &status);
            if (operation_completed) {
                if (!handle_send_work(status, count_workers_closed)) {
                    status_saved = status;
                    have_worker_ready = true;
                }
            }
        } else {
            if (system_tasks.size() != 0 || equation_tasks.size() != 0) {
                handle_send_work(status_saved, count_workers_closed);
                have_worker_ready = false;
            }
        }

        MPI_Test(&add_task_to_pool, &operation_completed, &status);
        if (operation_completed) {
            handle_add_task_to_pool(status, recv_equation_task_N);
        }

        MPI_Test(&worker_sent_system_result, &operation_completed, &status);
        if (operation_completed) {
            handle_got_system_result(system_solved_index);
        }
        
        MPI_Test(&worker_sent_equation_result, &operation_completed, &status);
        if (operation_completed) {
            handle_got_equation_result(status);
            worker_sent_equation_result = MPI_REQUEST_NULL;
        }
    }

    for (const auto& task : system_tasks_registry) {
        delete task;
    }

    MPI_Finalize();
    return 0;
}