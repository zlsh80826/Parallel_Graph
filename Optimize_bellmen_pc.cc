#include <bits/stdc++.h>

uint32_t vertex_num, edge_num;
uint32_t user_threads;
uint32_t source;
uint32_t* parent;
char* input_file;
char* output_file;

bool done = false;

std::ofstream outfs;

std::vector<std::thread> threads_pool;

class Vertex {
  public:
    uint32_t distance_to_source;
    std::vector<uint32_t> neighbors;
    std::vector<uint32_t> weight;
    std::vector<uint32_t> neighbors_distance;
    std::mutex mutex_v;

    void add_neighbor(uint32_t to, uint32_t w) {
        std::lock_guard<std::mutex> lock(mutex_v);
        neighbors.emplace_back(to);
        weight.emplace_back(w);
    }

    void AllocateNeighbor() { neighbors_distance.resize(neighbors.size()); }

    Vertex() : distance_to_source(std::numeric_limits<uint32_t>::max()) {}
};

std::vector<Vertex> vertex;

void FindPath(uint32_t root) {
    std::vector<uint32_t> reverse_vec;
    uint32_t current = root;
    reverse_vec.emplace_back(current);

    while (current != source) {
        current = parent[current];
        reverse_vec.emplace_back(current);
    }

    std::reverse(reverse_vec.begin(), reverse_vec.end());

    for (size_t i = 0; i < reverse_vec.size() - 1; ++i) {
        outfs << reverse_vec[i] << " ";
    }
    outfs << reverse_vec[reverse_vec.size() - 1] << std::endl;
}

void GoToLine(std::ifstream& file, uint32_t num) {
    file.seekg(std::ios::beg);
    for (ssize_t i = 0; i < num - 1; ++i) {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

void ThreadGoGo(int32_t thread_id) {
    int32_t extra = vertex_num % user_threads;
    int32_t chunk_size = vertex_num / user_threads + 1;
    int32_t work_load;
    int32_t offset;

    if (thread_id < extra) {
        offset = thread_id * chunk_size + 1;
        work_load = chunk_size;
    } else {
        offset = chunk_size * extra + ((thread_id - extra) * (chunk_size - 1)) + 1;
        work_load = chunk_size - 1;
    }

    for (ssize_t i = offset; i < offset + work_load; ++i) {
        auto& v = vertex[i];
        for (size_t j = 0; j < v.neighbors_distance.size(); ++j) {
            if (v.neighbors_distance[j] == std::numeric_limits<uint32_t>::max())
                continue;
            const uint32_t you = v.neighbors_distance[j] + v.weight[j];
            if (you < v.distance_to_source) {
                v.distance_to_source = you;
                parent[i] = v.neighbors[j];
                done = false;
            }
        }
    }
}

void ThreadYoLo(int32_t thread_id) {
    int32_t extra = vertex_num % user_threads;
    int32_t chunk_size = vertex_num / user_threads + 1;
    int32_t work_load;
    int32_t offset;

    if (thread_id < extra) {
        offset = thread_id * chunk_size + 1;
        work_load = chunk_size;
    } else {
        offset = chunk_size * extra + ((thread_id - extra) * (chunk_size - 1)) + 1;
        work_load = chunk_size - 1;
    }

    for (ssize_t i = offset; i < offset + work_load; ++i) {
        auto& v = vertex[i];
        for (size_t j = 0; j < v.neighbors.size(); ++j) {
            v.neighbors_distance[j] = vertex[v.neighbors[j]].distance_to_source;
        }
    }
}

void ReadFile(uint32_t thread_id) {
    std::ifstream infs;
    infs.open(input_file);

    uint32_t total_chunk = user_threads * 10 + (user_threads - 1) * user_threads / 2;
    uint32_t extra = edge_num % total_chunk;
    uint32_t work_load;
    uint32_t start_line;

    if (thread_id == 0) {
        start_line = 2;
        work_load = (edge_num / total_chunk) * (9 + user_threads - thread_id) + extra;
    } else {
        start_line = 2 + extra +
                     (edge_num / total_chunk) * (thread_id * 10 + (2 * user_threads - thread_id - 1) * thread_id / 2);
        work_load = (edge_num / total_chunk) * (9 + user_threads - thread_id);
    }

    GoToLine(infs, start_line);

    for (ssize_t i = 0; i < work_load; ++i) {
        uint32_t start, to, weight;
        infs >> start >> to >> weight;
        vertex[start].add_neighbor(to, weight);
        vertex[to].add_neighbor(start, weight);
    }

    infs.close();
}

int main(int argc, char** argv) {

    user_threads = std::stoi(argv[1]);
    source = std::stoi(argv[4]);
    input_file = argv[2];
    output_file = argv[3];

    auto IstartTimer = std::chrono::high_resolution_clock::now();

    std::ifstream infs;
    infs.open(input_file);

    infs >> vertex_num >> edge_num;

    infs.close();
    parent = new uint32_t[vertex_num + 1];
    for (ssize_t i = 1; i <= vertex_num; ++i) {
        parent[i] = source;
    }

    vertex = std::vector<Vertex>(vertex_num + 1);

    for (ssize_t i = 0; i < user_threads; ++i) {
        threads_pool.emplace_back(ReadFile, i);
    }

    for (auto& t : threads_pool) {
        t.join();
    }

    threads_pool.clear();

    for (auto& v : vertex) {
        v.AllocateNeighbor();
    }

    auto IendTimer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> Idiff = IendTimer - IstartTimer;
    std::cout << "Input took " << Idiff.count() << " Second" << std::endl;

    auto startTimer = std::chrono::high_resolution_clock::now();
    vertex[source].distance_to_source = 0;

    for (size_t i = 0; i < vertex[source].neighbors.size(); ++i) {
        vertex[vertex[source].neighbors[i]].distance_to_source = vertex[source].weight[i];
    }

    while (not done) {
        done = true;
        for (ssize_t i = 0; i < user_threads; ++i) {
            threads_pool.emplace_back(ThreadYoLo, i);
        }

        for (auto& t : threads_pool) {
            t.join();
        }

        threads_pool.clear();

        for (ssize_t i = 0; i < user_threads; ++i) {
            threads_pool.emplace_back(ThreadGoGo, i);
        }

        for (auto& t : threads_pool) {
            t.join();
        }

        threads_pool.clear();
    }

    auto endTimer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> diff = endTimer - startTimer;
    std::cout << "Bellmen compute took " << diff.count() << " Second" << std::endl;

    auto OstartTimer = std::chrono::high_resolution_clock::now();
    outfs.open(output_file);

    for (uint32_t i = 1; i <= vertex_num; ++i) {
        if (i == source)
            outfs << source << " " << source << std::endl;
        else
            FindPath(i);
    }
    outfs.close();

    auto OendTimer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> Odiff = OendTimer - OstartTimer;
    std::cout << "Output took " << Odiff.count() << " Second" << std::endl;
    return 0;
}
