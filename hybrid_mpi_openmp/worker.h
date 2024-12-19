#include <mpi.h>
#include <omp.h>

#include "utils.h"
#include "messages.h"

class Worker {
public:
    Worker(int id);
    ~Worker();

    void loop();
    void mpi_receive_work(int task_id);
    void work();
    void mpi_send_result();

private:
    struct JacobiTask *task;
    int id;
};
