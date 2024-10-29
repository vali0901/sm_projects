import sys
import numpy as np

# Reads a 10x10 matrix, 10x1 vector x and 10x1 vector b from stdin
# and prints the result of A*x, the vector b and the difference between
# A*x and b

def read_matrix(rows, cols):
    matrix = []
    for _ in range(rows):
        row = list(map(float, input().strip().split()))
        matrix.append(row)
    return np.array(matrix)

def read_vector(size):
    vector = list(map(float, input().strip().split()))
    return np.array(vector)

def main():
    print("Enter the 10x10 matrix A:")
    A = read_matrix(10, 10)

    print("Enter the 10x1 vector x:")
    x = read_vector(10)

    print("Enter the 10x1 vector b:")
    b = read_vector(10)

    Ax = np.dot(A, x)

    print("Result of A*x:")
    print(Ax)

    print("Vector b:")
    print(b)

    print("Differences (A*x - b):")
    print(Ax - b)

if __name__ == "__main__":
    main()