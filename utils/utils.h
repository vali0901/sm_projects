#pragma once
#include <stdlib.h>

struct system {
    int **A;
    int *b;
    int *x;
    int N;

    system(int N);

    ~system();

private:
    void fill_data();
};

/*
first dimension: the number of systems
second dimension: the size of the systems
*/
enum input_type {
    VERY_SMALL,
    SMALL_SMALL,
    SMALL_LARGE,
    LARGE_SMALL,
    LARGE_LARGE
};

struct input_gen {
    int no_systems;

    input_gen(int seed, input_type type);
    struct system **get_input();

private:
    int seed;
    input_type type;
    struct system **systems;
};
