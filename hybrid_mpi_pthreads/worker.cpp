#include "worker.h"

Worker::Worker(int id, int no_threads) : id(id) {
    this->no_threads = no_threads;
    task = new JacobiTask();
    threads = new pthread_t *[no_threads];
    semaphores = new sem_t *[no_threads];
    master_barrier = new pthread_barrier_t;
    worker_barrier = new pthread_barrier_t;
    data = new Worker::thread_data *[no_threads];

    should_stop = false;

    pthread_barrier_init(master_barrier, NULL, no_threads + 1);
    pthread_barrier_init(worker_barrier, NULL, no_threads);

    for (int i = 0; i < no_threads; i++) {
        threads[i] = new pthread_t;
        semaphores[i] = new sem_t;

        sem_init(semaphores[i], 0, 0);

        struct thread_data *data = new thread_data;
        this->data[i] = data;

        data->thread_id = i;
        data->semaphore = semaphores[i];
        data->master_barrier = master_barrier;
        data->worker_barrier = worker_barrier;
        data->sys = &this->task->sys;
        data->should_stop = &this->should_stop;
        data->x_new = &x_new;
        data->system_solved = &system_solved;
        data->no_threads = no_threads;

        pthread_create(threads[i], NULL, Worker::thread_work, data);
    }
}

Worker::~Worker() {
    delete task;
    for (int i = 0; i < no_threads; i++) {
        sem_destroy(semaphores[i]);
        pthread_join(*threads[i], NULL);
        delete semaphores[i];
        delete threads[i];
        delete data[i];
    }

    delete[] semaphores;
    delete[] threads;
    delete[] data;
    delete master_barrier;
    delete worker_barrier;
}

void Worker::mpi_receive_work(int task_id) {
    int N;
    task->task_id = task_id;

    MPI_Recv(&N, 1, MPI_INT, MASTER, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);

    task->sys = new struct system(N, false);
    x_new = new float[N];
    system_solved = false;

    for(int i = 0; i < task->sys->N; i++){
        MPI_Recv(task->sys->A[i], task->sys->N, MPI_FLOAT, MASTER, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
    }

    MPI_Recv(task->sys->b, task->sys->N, MPI_FLOAT, MASTER, DATA, MPI_COMM_WORLD, MPI_STATUS_IGNORE);
}

void Worker::work() {
    for(int i = 0 ; i < no_threads; i++) {

        sem_post(semaphores[i]);
    }

    pthread_barrier_wait(master_barrier);
}

void Worker::mpi_send_result() {
    MPI_Send(&task->task_id, 1, MPI_INT, MASTER, RESULT, MPI_COMM_WORLD);
    MPI_Send(task->sys->x, task->sys->N, MPI_FLOAT, MASTER, RESULT, MPI_COMM_WORLD);

    delete task->sys;
    delete[] x_new;
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
                this->should_stop = true;
                for(int i = 0; i < no_threads; i++)
                    sem_post(semaphores[i]);

                for(int i = 0; i < no_threads; i++)
                    pthread_join(*threads[i], NULL);

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

void *Worker::thread_work(void *arg) {
    struct thread_data *data = (struct thread_data *)arg;
    int thread_id = data->thread_id;
    sem_t *semaphore = data->semaphore;
    pthread_barrier_t *master_barrier = data->master_barrier;
    pthread_barrier_t *worker_barrier = data->worker_barrier;
    struct system **sys = data->sys;

    float tol = 1e-6;

    while(true) {
        sem_wait(semaphore);
        if(*data->should_stop) {
            break;
        }

        float **A = (*sys)->A;
        float *b = (*sys)->b;
        float *x = (*sys)->x;
        int N = (*sys)->N;
        int max_iter = N * 100;
        int needed_iter = max_iter;
        float *x_new = *data->x_new;

        int start = thread_id * N / data->no_threads;
        int end = (thread_id + 1) * N / data->no_threads;
        if(N < end) end = N;

        for(int iter = 0; iter < max_iter && !*data->system_solved; ++iter) {
            for(int i = start; i < end; ++i)  {
                x_new[i] = b[i];
                for(int j = 0; j < N; ++j) {
                    if(i == j)
                        continue;
                    x_new[i] -= A[i][j] * x[j];
                }
                x_new[i] /= A[i][i];
            }

            pthread_barrier_wait(worker_barrier);

            if(thread_id == 0) {
                double error = error_MSE(x, x_new, N);
                // double error = error_RIN(x, x_new, N);

                memcpy(x, x_new, N * sizeof(float));

                if(error < tol)
                    *data->system_solved = true;
            }

            pthread_barrier_wait(worker_barrier);
        }
        pthread_barrier_wait(master_barrier);
    }

    return NULL;
}