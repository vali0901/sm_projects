# SM_PROJECTS

This is the repository containing implementations of different
ways to parallelise the problem of solving N equation systems
using Jacobi Method. Each subfolder represents the specific
implementation containing source files and a Makefile (`make` for
building, `make run` for running, `make deep_clean` for clean up).

### Common elements
There is a utils folder containing some common structures used,
the most important one would be `input_gen` which is responsible
for generating different size inputs (using `input_type` enum
which for now only has 5 options, more will come for the
measurements stage) filled with random numbers based on a seed
(for tests consistency). The input is represented as an array of
systems, each system containing A matrix, b vector, and (unknown)
x vector. There are also two useful functions, one that solves a
system using linear time Jacobi-method and another one that
verifies if the system was solved correctly (computing the matrix
product Ax and comparing it to b).

### Implementations
For mpi, pthreads and openmp implementations, we chose to parallelise
the problem using a bag of tasks (where a task is solving a system).
The main reason behind this is that different systems can have
different sizes and complexities, so we can't just divide the systems
equally between workers. Similarly, if we think about mpi, in practice
the nodes do not have the same computational power, so for time
efficiency implementing a bag of tasks is the best approach. We do not
use any prioritizing mechanism, so the tasks are solved in the order
of generation (although in a real scenario a prioritizing mechanism
and a load balancer would be necessary).
For the hybrid implementations, the bag of tasks remains the same, but
now each task is solved dividing equally the work between the threads
(as we now know that the threads have, theoretically, the same
computational power). For these implementation, the bag of tasks code
remains the same (as the one from mpi implementation) but the worker's
work() funcion is changed respectively.
