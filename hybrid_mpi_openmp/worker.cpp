#include "worker.h"

Worker::Worker(int id) : id(id) {
    task = new JacobiTask();
}
Worker::~Worker() {
    delete task;
}

void Worker::mpi_receive_work(int task_id) {
    int N;
    task->task_id = task_id;

    MPI_Recv(&N, 1, MPI_INT, MASTER, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    task->sys = new struct system(N, false);

    for(int i = 0; i < task->sys->N; i++){
        // printf("%d\n", i);
        MPI_Recv(task->sys->A[i], task->sys->N, MPI_FLOAT, MASTER, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Recv(task->sys->b, task->sys->N, MPI_FLOAT, MASTER, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void Worker::work() {
    auto tol = 1e-6;
    auto sys = task->sys;
    float **A = sys->A;
    float *b = sys->b;
    float *x = sys->x;
    int N = sys->N;
    int max_iter = N * 100;

    float *x_new = new float[N];
    int needed_iter = max_iter;

    for (int iter = 0; iter < max_iter; ++iter) {
        #pragma omp for
        for (int i = 0; i < N; ++i) {
            x_new[i] = b[i];
            for (int j = 0; j < N; ++j) {
                if (i == j)
                    continue;

                x_new[i] -= A[i][j] * x[j];
            }
            x_new[i] /= A[i][i];
        }

        /* Need to choose what error funcion will we use */
        double error = error_MSE(x, x_new, N);
        // double error = error_RIN(x, x_new, N);

        memcpy(x, x_new, N * sizeof(float));

        if (error < tol) {
            needed_iter = iter;
            break;
        }
    }

    delete[] x_new;
    return;
}

void Worker::mpi_send_result() {
    MPI_Send(&task->task_id, 1, MPI_INT, MASTER, RESULT, MPI_COMM_WORLD);
    MPI_Send(task->sys->x, task->sys->N, MPI_FLOAT, MASTER, RESULT, MPI_COMM_WORLD);

    delete task->sys;
}

void Worker::loop() {
    int recv_work_signal = REQUEST_WORK;
    int done_signal = WORK_DONE;
    int data_recv;
    MPI_Status status;
    while(true) {
        MPI_Send(&recv_work_signal, 1, MPI_INT, MASTER, SIGNAL, MPI_COMM_WORLD);

        MPI_Recv(&data_recv, 1, MPI_INT, MASTER, MPI_ANY_TAG, MPI_COMM_WORLD, &status);
        switch (status.MPI_TAG)
        {
        case SIGNAL:
            if(data_recv == SHUT_DOWN) {
                printf("Worker %d shutting down\n", id);
                return;
            }
            continue;
            break;
        case DATA:
            mpi_receive_work(data_recv);
            work();
            MPI_Send(&done_signal, 1, MPI_INT, MASTER, SIGNAL, MPI_COMM_WORLD);
            mpi_send_result();
            break;
        default:
            break;
        }
    }
}