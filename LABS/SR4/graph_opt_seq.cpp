#include <cstdlib>
#include <iostream>
#include <vector>
#include <chrono>
#include <string>
#include <cassert>
#include <algorithm>

#define MAX_DIST 10
#define BLOCK_SIZE 64

class Graph {
private:
    std::vector<uint64_t> m_adj_matr;
    size_t m_num_verts;

public:
    Graph(const Graph& other) 
        : m_num_verts(other.m_num_verts),
          m_adj_matr(other.m_adj_matr) 
    {}

    Graph(size_t num_verts) 
        : m_num_verts(num_verts),
          m_adj_matr(num_verts * num_verts, MAX_DIST)
    {
        for(size_t i = 0; i < m_num_verts; ++i) {
            m_adj_matr[i*m_num_verts + i] = 0;
        }
    }

    void add_edge(size_t src, size_t dst, uint64_t dist) {
        assert(src < m_num_verts && dst < m_num_verts);
        m_adj_matr[src*m_num_verts + dst] = dist;
    }

    void print_adj_matr() const {
        std::cout << "\n";
        for(size_t i = 0; i < m_num_verts; ++i) {
            for(size_t j = 0; j < m_num_verts; ++j) {
                uint64_t val = m_adj_matr[i*m_num_verts + j];
                if(val < MAX_DIST) std::cout << val << ' ';
                else std::cout << "- ";
            }
            std::cout << "\n";
        }
    }

    Graph get_shortest_paths() const {
        Graph res(*this);
        const size_t n = m_num_verts;
        
        for(size_t kk = 0; kk < n; kk += BLOCK_SIZE) {
            const size_t k_end = std::min(kk + BLOCK_SIZE, n);
            
            for(size_t k = kk; k < k_end; ++k) {
                for(size_t i = 0; i < n; ++i) {
                    const uint64_t ik = res.m_adj_matr[i*n + k];
                    
                    for(size_t j = 0; j < n; ++j) {
                        const uint64_t sum = ik + res.m_adj_matr[k*n + j];
                        if(res.m_adj_matr[i*n + j] > sum) {
                            res.m_adj_matr[i*n + j] = sum;
                        }
                    }
                }
            }
        }
        
        return res;
    }

    size_t get_graph_size() const { return m_num_verts; }
};

void fill_rand_graph(Graph& graph) {
    size_t size = graph.get_graph_size();
    for(size_t i = 0; i < size; ++i) {
        for(size_t j = 0; j < size; ++j) {
            if(i != j) {
                graph.add_edge(i, j, (rand() % (MAX_DIST-1)) + 1);
            }
        }
    }
}

int main(int argc, char** argv)
{
    size_t vertex = 0;
    if (argc != 2) {
        std::cerr << "Usage: graph.exe vertex_num [default is 400].\n";
        vertex = 400;
    } else {
        vertex = std::stoll(argv[1]);
    }

    srand(17);
    Graph gr1(3);
    fill_rand_graph(gr1);
    gr1.print_adj_matr();
    std::cout << "\nTest run started\n";
    Graph res_gr(gr1.get_shortest_paths());

    std::cout << "Test results:\n";
    
    res_gr.print_adj_matr();

    Graph gr_test(vertex);
    fill_rand_graph(gr_test);
    std::cout << "\nPerf run started\n";
    std::chrono::time_point t_start = std::chrono::high_resolution_clock::now();
    Graph res_test_gr(gr_test.get_shortest_paths());
    std::chrono::time_point t_end = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> elapsed_time = t_end - t_start;
    std::cout << "Elapsed perf time: " << elapsed_time.count() << std::endl;

    return 0;
}
