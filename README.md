# SM_PROJECTS

This is the repository containing implementations of different
ways to parallelise the problem of solving N equation systems
using Jacobi Method. Each subfolder represents the specific
implementation containing source files, a Makefile (`make` for
building, `make run` for running, `make deep_clean` for clean up)
and a README detailing the implementation.

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
