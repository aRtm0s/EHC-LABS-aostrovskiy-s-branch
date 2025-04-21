#include <vector>
#include <iostream>
#include <random>
#include <string>
#include <chrono>
#include <omp.h>
#include <immintrin.h>

typedef double type;

class Gauss {
    type* matrix;    
    type* free_term;
    std::vector<type> variables;
    int n, m;

public:
    Gauss(int _n, int _m, type** _coefs, type* _free_terms)
        : n(_n), m(_m), variables(_n) {
        
        matrix = (type*)_mm_malloc(n*m*sizeof(type), 64);
        free_term = (type*)_mm_malloc(m*sizeof(type), 64);
        
        for (int i = 0; i < m; ++i) {
            for (int j = 0; j < n; ++j) {
                matrix[i*n + j] = _coefs[i][j];
            }
            free_term[i] = _free_terms[i];
        }
    }

    ~Gauss() {
        _mm_free(matrix);
        _mm_free(free_term);
    }

    double solve() {
        std::chrono::time_point start = std::chrono::high_resolution_clock::now();
        
        omp_set_num_threads(6);
        forward();
        backward();
        get_answer();
        
        std::chrono::time_point end = std::chrono::high_resolution_clock::now();
        return std::chrono::duration<double>(end - start).count();
    }

    void forward() {
        for (int i = 0; i < m; ++i) {
            type pivot = matrix[i*n + i];
            
            #pragma omp simd aligned(matrix, free_term:64)
            for (int j = i; j < n; ++j) {
                matrix[i*n + j] /= pivot;
            }
            free_term[i] /= pivot;

            #pragma omp parallel for schedule(dynamic) //nowait
            for (int j = i + 1; j < m; ++j) {
                type factor = matrix[j*n + i];
                #pragma omp simd aligned(matrix:64)
                for (int k = i; k < n; ++k) {
                    matrix[j*n + k] -= matrix[i*n + k] * factor;
                }
                free_term[j] -= free_term[i] * factor;
            }
        }
    }

    void backward() {
        for (int i = m - 1; i >= 0; --i) {
            type pivot = matrix[i*n + i];
            
            #pragma omp simd aligned(matrix, free_term:64)
            for (int j = i; j < n; ++j) {
                matrix[i*n + j] /= pivot;
            }
            free_term[i] /= pivot;

            #pragma omp parallel for schedule(dynamic) //nowait
            for (int j = i - 1; j >= 0; --j) {
                type factor = matrix[j*n + i];
                #pragma omp simd aligned(matrix:64)
                for (int k = i; k < n; ++k) {
                    matrix[j*n + k] -= matrix[i*n + k] * factor;
                }
                free_term[j] -= free_term[i] * factor;
            }
        }
    }

    void get_answer() {
        #pragma omp simd
        for (int i = 0; i < n; ++i) {
            variables[i] = free_term[i];
        }
    }

    type& operator[](int ind) {
        return variables[ind];
    }
};

bool correctness_test_run() {
    type line_1[] = {1, 2, 3};
    type line_2[] = {3, 5, 7};
    type line_3[] = {1, 3, 4};

    type free_terms[] = {3, 0, 1};
    type answer[] = {-4, -13, 11};

    type* coefs[] = {line_1, line_2, line_3};

    Gauss method(3, 3, coefs, free_terms);
    double test_time = method.solve();

    for (int i = 0; i < 3; ++i) {
        if (std::abs(answer[i] - method[i]) > 1e-9) {
            std::cout << "Small test failed. Calculations are incorrect.\n";
            return false;
        }
    }
    std::cout << "Small test passed.\n";
    return true;
}

bool performance_run(int task_size) {
    type** coefs = new type*[task_size];
    type* free_term = new type[task_size];

    std::mt19937 generator;
    std::uniform_real_distribution<type> dist(-10.0, 10.0);

    for (int i = 0; i < task_size; ++i) {
        coefs[i] = new type[task_size];
        for (int j = 0; j < task_size; ++j) {
            coefs[i][j] = dist(generator);
            if (i == j) coefs[i][j] += 100.0;
        }
        free_term[i] = dist(generator);
    }

    Gauss method(task_size, task_size, coefs, free_term);
    double seq_time = method.solve();
    std::cout << "Elapsed sequential Gaussian time: " << seq_time << " seconds\n";

    for (int i = 0; i < task_size; ++i) {
        delete[] coefs[i];
    }
    delete[] coefs;
    delete[] free_term;
    return true;
}

int main(int argc, char** argv) {
    int task_size = 1000;

    if (argc > 1) {
        task_size = std::stoi(argv[1]);
        std::cout << "Task size set to: " << task_size << std::endl;
    } else {
        std::cout << "Using default task size: " << task_size << std::endl;
    }

    correctness_test_run();
    performance_run(task_size);
    return 0;
}