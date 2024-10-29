#include "utils.h"
#include <string.h>

system::system(int N) : N(N) {
    A = new int*[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new int[N];
    }
    b = new int[N];
    x = new float[N];
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

void system::fill_data() {
    for (int i = 0; i < N; ++i) {
        for (int j = 0; j < N; ++j) {
            A[i][j] = rand() % 51 - 25;
        }
        b[i] = rand() % 101 - 50;
        x[i] = 0;
    }

    for (int i = 0; i < N; ++i) {
        int sum = 0;
        for (int j = 0; j < N; ++j) {
            sum += A[i][j];
        }
        A[i][i] = sum;
    }
}

input_gen::input_gen(int seed, input_type type) : seed(seed), type(type) {
    srand(seed);
    switch (type) {
        case VERY_SMALL:
            no_systems = 1;
            systems = new struct system*[1]{new struct system(10)};
            break;
        case SMALL_SMALL:
            no_systems = 100;
            systems = new struct system*[100]{new struct system(100)};
            break;
        case SMALL_LARGE:
            no_systems = 100;
            systems = new struct system*[100]{new struct system(1000)};
            break;
        case LARGE_SMALL:
            no_systems = 1000;
            systems = new struct system*[1000]{new struct system(100)};
            break;
        case LARGE_LARGE:
            no_systems = 1000;
            systems = new struct system*[1000]{new struct system(1000)};
            break;
    }
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

void simple_jacobi(int **A, int *b, float *x, int N, int max_iter, double tol) {
    float *x_new = new float[N];
    for (int i = 0; i < N; ++i) {
        x[i] = 0;
    }

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

        memcpy(x, x_new, N * sizeof(int));

        /* Need to choose what error funcion will we use */

        double error = error_MSE(x, x_new, N);
        // double error = error_RIN(x, x_new, N);

        if (error < tol)
            break;
    }

    delete[] x_new;
}
