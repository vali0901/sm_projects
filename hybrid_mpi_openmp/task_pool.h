#include <mpi.h>

#include "utils.h"
#include "messages.h"

class TaskPool {
public:
    TaskPool(int no_tasks, int no_workers);
    ~TaskPool();

    void load_jacobi_tasks(struct system **systems);
    void loop();


private:
    int no_tasks;
    int no_workers;
    bool *task_consumed;
    bool *worker_active;
    int no_workers_active;
    int available_task;
    std::vector<struct JacobiTask*> jacobi_tasks;

    void mpi_send_task(int worker_id);
    void mpi_receive_task(int worker_id);
};
