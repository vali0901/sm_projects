#include <iostream>
#include <pthread.h>
#include "../utils/utils.h"

#define DEBUG 1

pthread_mutex_t mutex;
int N = 10;

struct system** systems;
int count_systems;
int next_system_unsolved;

void* thread_run(void *args) {
    while (true) {
        pthread_mutex_lock(&mutex);
        if (next_system_unsolved >= count_systems) {

            pthread_mutex_unlock(&mutex);
            return nullptr;
        }
        struct system* current_system = systems[next_system_unsolved++];
    
        pthread_mutex_unlock(&mutex);

        simple_jacobi(current_system);
    }
    return nullptr;
};

int main() {
    pthread_t threads[N];
    int thread_ids[N];
    pthread_mutex_init(&mutex, nullptr);

    struct input_gen input(0, input_type::VERY_SMALL);
    systems = input.get_input();
    count_systems = input.no_systems;


    for (int i = 0; i < N; ++i) {
        thread_ids[i] = i;
        pthread_create(&threads[i], nullptr, thread_run, &thread_ids[i]);
    }

    thread_run(nullptr);

    for (int i = 0; i < N; ++i) {
        pthread_join(threads[i], nullptr);
    }

    #ifdef DEBUG
    std::cout << "Below you should see if any system wasn't solve successfully\n";
    for (int i = 0; i < count_systems; i++) {
        struct system* sys = systems[i];
        if (!system_solved(sys))
            std::cout << "system " + i << "is not solved" << std::endl;
    }
    #endif

    pthread_mutex_destroy(&mutex);

    return 0;
}