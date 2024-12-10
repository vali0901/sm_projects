#pragma once
#include <stdlib.h>
#include <string.h>
#include <iostream>
#include <bits/stdc++.h>

struct system {
    float **A;
    float *b;
    float *x;
    int N;

    system(int N, bool fill = true);
    void print();
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

int simple_jacobi(struct system *sys, double tol = 1e-6);
bool system_solved(struct system *sys, double tol = 0.1);
