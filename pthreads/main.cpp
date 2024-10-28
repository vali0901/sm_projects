#include <iostream>
#include <pthread.h>

const int N = 4;
const int NUMBERS = 10;

int sum = 0;
pthread_mutex_t mutex;

void* compute_sum(void* arg) {
    int thread_id = *(int*)arg;
    int local_sum = 0;

    for (int i = thread_id; i < NUMBERS; i += N) {
        local_sum += i + 1;
    }

    pthread_mutex_lock(&mutex);
    sum += local_sum;
    pthread_mutex_unlock(&mutex);

    return nullptr;
}

int main() {
    pthread_t threads[N];
    int thread_ids[N];
    pthread_mutex_init(&mutex, nullptr);

    for (int i = 0; i < N; ++i) {
        thread_ids[i] = i;
        pthread_create(&threads[i], nullptr, compute_sum, &thread_ids[i]);
    }

    for (int i = 0; i < N; ++i) {
        pthread_join(threads[i], nullptr);
    }

    pthread_mutex_destroy(&mutex);

    std::cout << "Sum of first " << NUMBERS << " numbers is: " << sum << std::endl;

    return 0;
}