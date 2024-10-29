#include "utils.h"

system::system(int N) : N(N) {
    A = new int*[N];
    for (int i = 0; i < N; ++i) {
        A[i] = new int[N];
    }
    b = new int[N];
    x = new int[N];
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

