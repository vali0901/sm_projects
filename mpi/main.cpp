#include <mpi.h>
#include <iostream>
#include "../utils/utils.h"

#define EQUATION_SOLVER = 1
#define SYSTEM_SOLVER = 2

#define RANK_ROOT 0
#define EPS 0.005
struct Task {
    
};

struct Equation_solver_task : Task {
    float *A_row;
    float *x;
    float b;

    int N;
    int x_index;
    int system_ID;
};

struct System_solver_task : Task {
    float **A;
    float *x;
    float *b;
    int N;
    int ID;

    bool solved = false;
    int equations_solved = 0;
    int assigned_process_rank = -1;

    System_solver_task(float **A, float *x, float *b, int N, int ID) 
        : A(A), x(x), b(b), N(N), ID(ID) {
    };
};

auto equation_tasks = std::queue<Equation_solver_task*>();
auto system_tasks = std::queue<System_solver_task*>();
auto system_tasks_registry = std::vector<System_solver_task*>();

int world_size;
int count_systems_solved = 0;

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

    Equation_solver_task *task = new Equation_solver_task();
    task->N = N;
    task->A_row = new float[N];
    task->b = 0;
    task->x = new float[N];
    task->x_index = 0;
    task->system_ID = 0;

    MPI_Recv(&(task->system_ID), 1, MPI_INT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&(task->x_index), 1, MPI_INT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(task->A_row, N, MPI_FLOAT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(task->x, N, MPI_FLOAT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&(task->b), 1, MPI_FLOAT, source, ch, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    // //std::cout << ">>>> SYSTEM_ID = " << task->system_ID << std::endl;
    return task;
}

System_solver_task read_system_task(int N) {
    float **A = new float*[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new float[N];
    }
    float *b = new float[N];
    float *x = new float[N];
    int ID = -1;
    
    System_solver_task task(A, x, b, N, ID);
    
    MPI_Recv(&(task.ID), 1, MPI_INT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    for (int i = 0; i < N; i++) {
        MPI_Recv(task.A[i], N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Recv(task.x, N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(task.b, N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    return task;
};

float solve_equation_task(const Equation_solver_task *task) {
    // //std::cout << "eq = \n A= ";
    // for (int i = 0; i < task->N; i++)
    //     //std::cout << task->A_row[i] << ' ';
    // //std::cout << "\n b= " << task->b << "\nx = ";

    // for (int i = 0; i < task->N; i++)
    //     //std::cout << task->x[i] << ' ';
    // //std::cout << "\n";

    float result = task->b;

    for (int j = 0; j < task->N; ++j) {
        if (task->x_index == j)
            continue;

        result -= task->A_row[j] * task->x[j];
    }
    result /= task->A_row[task->x_index];

    return result;
};

void get_system_latest_iteration(int ID, int N, float* result) {
    //std::cout << "before getting the latest" << N << '\n';
    MPI_Recv(result, N, MPI_FLOAT, RANK_ROOT, GET_SYSTEM_NEW_RESULTS, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    //std::cout << "after getting the latest\n";

}

// TODO: workaround to not move the function upward, add definition in header file 
void send_equation_task(int destination, Equation_solver_task *task, Channels ch);

void solve_system_task(const System_solver_task &task) {
    int max_iter = task.N * 10;
    int needed_iter = max_iter;
    float *x_new = new float[task.N];
    int rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &rank);
    //std::cout << "[" << rank << " worker] " << "started working on system\n";
    

    for (int iter = 0; iter < max_iter; ++iter) {
        
        double error = error_MSE(task.x, x_new, task.N);
        for (int i = 0; i < task.N; ++i) {
            auto equation_task = new Equation_solver_task();
            equation_task->N = task.N;
            equation_task->x_index = i;
            equation_task->b = task.b[i];
            equation_task->x = task.x;
            equation_task->A_row = task.A[i];
            equation_task->system_ID = task.ID;

            // std::cout << "[" << rank << " worker] " << "before sending the task\n";
            send_equation_task(RANK_ROOT, equation_task, ADD_TASK_TO_POOL);
            // std::cout << "[" << rank << " worker] " << "has sent a task\n";
            delete equation_task;
        }
        // std::cout << "[" << rank << " worker] " << "TRIES to get latest\n";
        get_system_latest_iteration(task.ID, task.N, x_new);
        // std::cout << "[" << rank << " worker] " << "got latest latest\n";

        memcpy(task.x, x_new, task.N * sizeof(float));

        if (error < EPS) {
            needed_iter = iter;
            break;
        }
    }
    delete[] x_new;
};

void send_equation_result(float result, int system_ID, int x_index) {
    // if (system_ID)
    // std::cout << "[=====] SEND equation result systemID = " << system_ID  << " x_index = " << x_index << " result = " << result << std::endl << std::endl;
    int dummy = -1;
    MPI_Send(&dummy, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);

    MPI_Send(&system_ID, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);

    MPI_Send(&x_index, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);

    MPI_Send(&result, 1, MPI_FLOAT, RANK_ROOT, SEND_EQUATION_SOLVED, MPI_COMM_WORLD);
};

void send_system_result(System_solver_task task) {
    MPI_Send(&(task.ID), 1, MPI_INT, RANK_ROOT, SEND_SYSTEM_SOLVED, MPI_COMM_WORLD);

    // MPI_Send(task.x, task.N, MPI_FLOAT, RANK_ROOT, SEND_SYSTEM_SOLVED, MPI_COMM_WORLD);
};

void worker_process_run() {
    int N;
    int dummy;
    MPI_Status status;
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
            // //std::cout << "[" << rank << " worker] " << "got a retry\n";
            continue;
        }
        //std::cout << "here0?\n";

        MPI_Recv(&tag, 1, MPI_INT, RANK_ROOT, REQUEST_WORK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
        //std::cout << "here1?\n";

        if (tag == SEND_EQUATION_TASK) {
            MPI_Recv(&N, 1, MPI_INT, RANK_ROOT, SEND_EQUATION_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
            //std::cout << "here?\n";
            auto task = read_equation_task(N, RANK_ROOT, SEND_EQUATION_TASK);
            auto result = solve_equation_task(task);
            send_equation_result(result, task->system_ID, task->x_index);
            //std::cout << "[" << rank << " worker] " << "solved equation\n";

            // TODO: will move the frees to the corresponding destructor of task
            delete[] task->A_row;
            delete[] task->x;
            delete task;
        } else {
            if (tag == SEND_SYSTEM_TASK) {
                // std::cout << "[" << rank << " worker] " << "got system to solve\n";

                MPI_Recv(&N, 1, MPI_INT, RANK_ROOT, SEND_SYSTEM_TASK, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

                auto task = read_system_task(N);
                // std::cout << "[" << rank << " worker] " << "BEFORE" << task.ID << "\n";
                solve_system_task(task);
                // std::cout << "[" << rank << " worker] " << "WANTS to send system solved\n";
                
                // std::cout << "[" << rank << " worker] " << "AFTER" << task.ID << "\n";
                send_system_result(task);

                // std::cout << "[" << rank << " worker] " << "solved system\n";

                // TODO: will move the frees to the corresponding destructor of task
                for (int i = 0; i < task.N; i++)
                    delete[] task.A[i];

                delete[] task.A;
                delete[] task.b;
                delete[] task.x;
            } else {
                //std::cout <<"[ERROR] got data where i shouldnt\n";
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

void handle_send_work(MPI_Status status, int &count_workers_closed) {
    int destination = status.MPI_SOURCE;
    // we prioritize the tasks that solve equations - solves the problem of starvation    
    if (!equation_tasks.empty()) {
        Get_work_return result = Get_work_return::SUCCESS;
        MPI_Send(&result, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);

        Channels tag = SEND_EQUATION_TASK;
        MPI_Send(&tag, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);
        
        send_equation_task(destination, equation_tasks.front(), SEND_EQUATION_TASK);
        equation_tasks.pop();
        return;
    }
    
    if (!system_tasks.empty()) {
        // std::cout << "[ROOT] tries to send system task\n";
        Get_work_return result = Get_work_return::SUCCESS;
        MPI_Send(&result, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);

        Channels tag = SEND_SYSTEM_TASK;
        MPI_Send(&tag, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);

        send_system_task(destination, system_tasks.front(), SEND_SYSTEM_TASK);
        system_tasks.front()->assigned_process_rank = destination;
        system_tasks.pop();
        return;
    }

    // TODO CHECK IF THERE IS NO WORK IN PROGRESS
    // TODO if work available return false
    if (!operations_finished()) {
        // //std::cout << "not all ops finished\n";
        int retry = Get_work_return::RETRY;
        MPI_Send(&retry, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);

        return;
    }

    int close = Get_work_return::CLOSE;
    MPI_Send(&close, 1, MPI_INT, destination, REQUEST_WORK, MPI_COMM_WORLD);
    count_workers_closed++;
}

void send_updated_values_for_equation(int destination, float *data, int N) {
    //std::cout << "send_updated_values dest= " << destination << std::endl;
    MPI_Send(data, N, MPI_FLOAT, destination, GET_SYSTEM_NEW_RESULTS, MPI_COMM_WORLD);
}


void handle_got_equation_result(MPI_Status status) {
    //systemID, x_index, result
    int x_index;
    float result;
    int system_ID;
    // TODO remove the comment, is commented because we already got the sys id
    MPI_Recv(&system_ID, 1, MPI_INT, status.MPI_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    MPI_Recv(&x_index, 1, MPI_INT, status.MPI_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    
    MPI_Recv(&result, 1, MPI_FLOAT, status.MPI_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    // std::cout << "[Error]: got EQUATION result existent system ID " << system_ID << " x_index = " << x_index << " result = " << result<< "\n";

    for (const auto& task : system_tasks_registry) {
        if (task->ID == system_ID) {
            if (task->solved) {
                std::cout << "[ERROR] task already solved\n";
            }
            task->x[x_index] = result;
            task->equations_solved++;
            if (task->equations_solved == task->N) {
                task->equations_solved = 0;;
                send_updated_values_for_equation(task->assigned_process_rank, task->x, task->N);
            }

            //std::cout << "[ ROOT ] got equation result sys_ID = " << system_ID << " x_index = " << x_index << " result = " << result << std::endl;
            return;
        }
    }

}

void handle_got_system_result(MPI_Status status, int system_ID) {
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

    //std::cout << "[ERROR] got the result of a non valid system_id = " << system_ID << std::endl;
}

int main(int argc, char** argv) {
    MPI_Init(&argc, &argv);;

    MPI_Comm_size(MPI_COMM_WORLD, &world_size);

    int world_rank;
    MPI_Comm_rank(MPI_COMM_WORLD, &world_rank);

    if (world_rank != RANK_ROOT) {
        worker_process_run();
        // std::cout << "[" << world_rank << " worker] finished work" << std::endl;

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
        // for (int m = 0; m < systems[i]->N; m++) {
        //     for (int n = 0; n < systems[i]->N; n++)
        //         std::cout << systems[i]->A[m][n] << ' ';
        //     std::cout << '\n';
        // }
        // std::cout << '\n';
        // for (int j = 0; j < systems[i]->N; j++) {
        //     std::cout << systems[i]->b[j] << ' ';
        // }
        // std::cout << '\n';
        
        // TODO remove break, this only generates one system
        // break;
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
    while (count_workers_closed < (world_size - 1)) {
        // std::cout << count_systems_solved << ' ' << count_systems_solved << std::endl;

        if (get_work_request == MPI_REQUEST_NULL)
            MPI_Irecv(&dummy, 1, MPI_INT, MPI_ANY_SOURCE, REQUEST_WORK, MPI_COMM_WORLD, &get_work_request);

        if (add_task_to_pool == MPI_REQUEST_NULL)
            MPI_Irecv(&recv_equation_task_N, 1, MPI_INT, MPI_ANY_SOURCE, ADD_TASK_TO_POOL, MPI_COMM_WORLD, &add_task_to_pool);

        if (worker_sent_system_result == MPI_REQUEST_NULL)
            MPI_Irecv(&system_solved_index, 1, MPI_INT, MPI_ANY_SOURCE, SEND_SYSTEM_SOLVED, MPI_COMM_WORLD, &worker_sent_system_result);

        if (worker_sent_equation_result == MPI_REQUEST_NULL)
            MPI_Irecv(&equation_solved_dummy, 1, MPI_INT, MPI_ANY_SOURCE, SEND_EQUATION_SOLVED, MPI_COMM_WORLD, &worker_sent_equation_result);


        MPI_Test(&get_work_request, &operation_completed, &status);
        if (operation_completed) {
            handle_send_work(status, count_workers_closed);
        }

        MPI_Test(&add_task_to_pool, &operation_completed, &status);
        if (operation_completed) {
            handle_add_task_to_pool(status, recv_equation_task_N);
        }

        MPI_Test(&worker_sent_system_result, &operation_completed, &status);
        if (operation_completed) {
            handle_got_system_result(status, system_solved_index);
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
    //std::cout << "[ROOT] finished work" << std::endl;
    MPI_Finalize();
    return 0;
}