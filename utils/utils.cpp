#include "utils.h"


system::system(int N, bool fill) : N(N) {
    A = new float*[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new float[N];
    }
    b = new float[N];
    x = new float[N];
    if(fill)
        fill_data();
}

system::~system() {
    for (int i = 0; i < N; ++i) {
        delete[] A[i];
    }
    delete[] A;
    delete[] b;
    delete[] x;
}

void system::print() {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j)
            std::cout << A[i][j] << " ";
        std::cout << '\n';
    }

    for (int i = 0; i < N; ++i)
        std::cout << x[i] << " ";

    std::cout << '\n';

    for (int i = 0; i < N; ++i)
        std::cout << b[i] << " ";

    std::cout << '\n';
}

void system::fill_data() {
    std::uniform_real_distribution<float> distribution(-25.0, 25.0);
    std::default_random_engine generator;

    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = distribution(generator);
        }
        b[i] = distribution(generator);
        x[i] = 0;
    }

    for (int i = 0; i < N; ++i) {
        int sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += abs(A[i][j]);
        }
        A[i][i] = sum;
    }
}

input_gen::input_gen(int seed, input_type type) : seed(seed), type(type) {
    srand(seed);
    int system_size = 0;
    switch (type) {
        case VERY_SMALL:
            no_systems = 1;
            system_size = 10;

            break;
        case SMALL_SMALL:
            no_systems = 100;
            system_size = 100;

            break;
        case SMALL_LARGE:
            no_systems = 100;
            system_size = 1000;

            break;
        case LARGE_SMALL:
            no_systems = 1000;
            system_size = 100;

            break;
        case LARGE_LARGE:
            no_systems = 1000;
            system_size = 1000;

            break;
    }

    systems = new struct system*[no_systems];
    for (int i = 0; i < no_systems; i++)
        systems[i] = new struct system(system_size);
}

struct system **input_gen::get_input() {
    return systems;
}

float error_MSE(float *x, float *y, int N) {
    float error = 0;
    for (int i = 0; i < N; ++i) {
        error += (x[i] - y[i]) * (x[i] - y[i]);
    }
    return error / N;
}

/* Relative Error using Infinity Norm */
float error_RIN(float *x, float *y, int N) {
    float relative_norm = 0;
    float norm = 0;
    for (int i = 0; i < N; ++i) {
        float temp = abs(x[i] - y[i]);
        if (temp > relative_norm) {
            relative_norm = temp;
        }
        if (abs(x[i]) > norm) {
            norm = abs(x[i]);
        }
    }
    return relative_norm / norm;
}

int simple_jacobi(struct system *sys, double tol) {
    float **A = sys->A;
    float *b = sys->b;
    float *x = sys->x;
    int N = sys->N;
    int max_iter = N * 100;

    float *x_new = new float[N];
    int needed_iter = max_iter;

    for (int iter = 0; iter < max_iter; ++iter) {
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
    return needed_iter;
}

bool system_solved(struct system *sys, double tol) {
    float **A = sys->A;
    float *b = sys->b;
    float *x = sys->x;
    int N = sys->N;

    for (int i = 0; i < N; ++i) {
        float sum = 0;
        for (int j = 0; j < N; ++j)
            sum += A[i][j] * x[j];
        if (abs(sum - b[i]) > tol)
            return false;
    }

    return true;
}

void log_solution_correctness(struct system **systems, int N, std::string solver_name) {
    std::cout << "Solver: " << solver_name << '\n';
    bool is_correct = true;
    for (int i = 0; i < N; i++) {
        struct system* sys = systems[i];
        if (!system_solved(sys)) {
            std::cout << "system " + i << "is not solved correctly" << std::endl;
            is_correct = false;
        }
    }

    if (is_correct)
        std::cout << "Solved systems correctly\n\n";
}