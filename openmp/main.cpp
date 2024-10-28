#include <iostream>
#include <omp.h>

int main(int argc, char *argv[]) {
    int sum = 0;

    #pragma omp parallel for reduction(+:sum)
    for (int i = 1; i <= 10; ++i) {
        sum += i;
    }

    std::cout << "Sum = " << sum << "\n";
    return 0;
}
