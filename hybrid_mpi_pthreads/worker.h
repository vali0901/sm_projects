#include <mpi.h>
#include <pthread.h>
#include <semaphore.h>

#include "utils.h"
#include "messages.h"

#define NO_THREADS 4

class Worker {
public:
    Worker(int id, int no_threads = NO_THREADS);
    ~Worker();

    void loop();
    void mpi_receive_work(int task_id);
    void work();
    void mpi_send_result();

private:
    struct thread_data {
        int thread_id;
        int no_threads;

        struct system **sys;
        float **x_new;
        bool *system_solved;

        sem_t *semaphore;
        pthread_barrier_t *master_barrier;
        pthread_barrier_t *worker_barrier;

        bool *should_stop;
    };

    int no_threads;
    int id;
    struct JacobiTask *task;

    pthread_t **threads;

    sem_t **semaphores;
    pthread_barrier_t *master_barrier;
    pthread_barrier_t *worker_barrier;
    struct thread_data **data;

    bool should_stop;
    float *x_new;
    bool system_solved;

    static void *thread_work(void *arg);
};
