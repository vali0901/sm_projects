#include <iostream>
#include <omp.h>
#include "utils.h"

int main(int argc, char *argv[]) {
    
    struct input_gen input(0, input_type::SMALL_SMALL);
    struct system **systems = input.get_input();

    #pragma omp parallel
    {
        #pragma omp single
        {
            for (int i = 0; i < input.no_systems; i++) {
                #pragma omp task
                {
                    simple_jacobi(systems[i]);
                }
            }
            #pragma omp taskwait
        }
    }

    #ifdef DEBUG
        log_solution_correctness(systems, input.no_systems, "openmp");
    #endif

    for (int i = 0; i < input.no_systems; i++)
        delete systems[i];
       
    return 0;
}
