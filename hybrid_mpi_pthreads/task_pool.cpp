#include "task_pool.h"

TaskPool::TaskPool(int no_tasks, int no_workers) {
    this->no_tasks = no_tasks;
    this->no_workers = no_workers;
    this->no_workers_active = no_workers;
    this->task_consumed = new bool[no_tasks];
    memset(this->task_consumed, 0, no_tasks * sizeof(bool));
    this->worker_active = new bool[no_workers];
    memset(this->worker_active, 1, no_workers * sizeof(bool));
    available_task = 0;
}

TaskPool::~TaskPool() {
    delete[] this->task_consumed;
    for(auto task : this->jacobi_tasks)
        delete task;
    this->jacobi_tasks.clear();
    delete[] this->worker_active;
}

void TaskPool::load_jacobi_tasks(struct system **systems) {
    for(int i = 0; i < this->no_tasks; i++) {
        struct JacobiTask *task = new JacobiTask();
        task->task_id = i;
        task->sys = systems[i];
        this->jacobi_tasks.push_back(task);
    }
}

void TaskPool::mpi_send_task(int worker_id) {
    int task_id = this->available_task;

    if(task_id >= no_tasks) {
        int data = SHUT_DOWN;
        MPI_Send(&data, 1, MPI_INT, worker_id, SIGNAL, MPI_COMM_WORLD);
        worker_active[worker_id] = false;
        no_workers_active--;
        return;
    }

    struct JacobiTask *task = this->jacobi_tasks[task_id];
    MPI_Send(&task->task_id, 1, MPI_INT, worker_id, DATA, MPI_COMM_WORLD);
    MPI_Send(&task->sys->N, 1, MPI_INT, worker_id, DATA, MPI_COMM_WORLD);
    for(int i = 0; i < task->sys->N; i++)
        MPI_Send(task->sys->A[i], task->sys->N, MPI_FLOAT, worker_id, DATA, MPI_COMM_WORLD);
    MPI_Send(task->sys->b, task->sys->N, MPI_FLOAT, worker_id, DATA, MPI_COMM_WORLD);
    this->task_consumed[task_id] = true;

    this->available_task++;
}

void TaskPool::mpi_receive_task(int worker_id) {
    int task_id;
    MPI_Recv(&task_id, 1, MPI_INT, worker_id, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    if(task_id < 0 || task_id >= no_tasks) {
        return;
    }
    struct JacobiTask *task = this->jacobi_tasks[task_id];
    MPI_Recv(task->sys->x, task->sys->N, MPI_FLOAT, worker_id, RESULT, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void TaskPool::loop() {
    while(true) {
        int worker_id;
        int signal;
        MPI_Status status;
        MPI_Recv(&signal, 1, MPI_INT, MPI_ANY_SOURCE, SIGNAL, MPI_COMM_WORLD, &status);
        switch (signal)
        {
        case REQUEST_WORK:
            worker_id = status.MPI_SOURCE;
            mpi_send_task(worker_id);
            break;
        case WORK_DONE:
            worker_id = status.MPI_SOURCE;
            mpi_receive_task(worker_id);
            break;
        default:
            break;
        }
        if(no_workers_active == 0)
            break;
    }
}

