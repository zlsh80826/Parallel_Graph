#include <sys/syscall.h>
#include <unistd.h>
// #include <boost/iostreams/device/mapped_file.hpp>
// #include <boost/iostreams/stream.hpp>
#include <bits/stdc++.h>

uint32_t vertex_num, edge_num;
uint32_t user_threads;
uint32_t source;
uint32_t* parent;
char* input_file;
char* output_file;

std::ofstream outfs;

std::vector<std::thread> threads_pool;
std::mutex file_mutex;

class Edge {
  public:
    uint32_t to;
    uint32_t weight;

    Edge(uint32_t to, uint32_t weight) : to(to), weight(weight) {}

    bool operator<(const Edge& e) const { return this->weight < e.weight; }
};

class Vertex {
  public:
    uint32_t distance_to_source;
    std::vector<Edge> neighbors;
    std::mutex mutex_v;

    void add_neighbor(uint32_t to, uint32_t weight) {
        std::lock_guard<std::mutex> lock(mutex_v);
        neighbors.emplace_back(to, weight);
    }

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

/*void GoToLine(boost::iostreams::stream<boost::iostreams::mapped_file_source>& file, uint32_t num) {
    file.seekg(std::ios::beg);
    for (ssize_t i = 0; i < num - 1; ++ i) {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}*/

void ReadFile(uint32_t thread_id) {
    // std::cout << syscall(SYS_gettid) << std::endl;

    // boost::iostreams::mapped_file_source mmap(input_file);
    // boost::iostreams::stream<boost::iostreams::mapped_file_source> infs(mmap);
    std::ifstream infs;
    infs.open(input_file);
    // uint32_t extra = edge_num % user_threads;

    uint32_t total_chunk = user_threads * 10 + (user_threads - 1) * user_threads / 2;
    uint32_t extra = edge_num % total_chunk;

    // uint32_t chunk_size = edge_num / user_threads + 1;
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

    /*if (thread_id < extra) {
        start_line = chunk_size * thread_id + 2;
        work_load = chunk_size;
    } else {
        start_line = (chunk_size * extra) + (chunk_size - 1) * (thread_id - extra) + 2;
        work_load = chunk_size - 1;
    }*/

    GoToLine(infs, start_line);
    // std::cout << "Go to " << start_line << " done" << std::endl;

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
        threads_pool.push_back(std::thread(ReadFile, user_threads - i - 1));
    }

    for (auto& t : threads_pool) {
        t.join();
    }

    auto IendTimer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> Idiff = IendTimer - IstartTimer;
    std::cout << "Input took " << Idiff.count() << " Second" << std::endl;

    auto startTimer = std::chrono::high_resolution_clock::now();
    vertex[source].distance_to_source = 0;

    // uint32_t me;
    // uint32_t you;

    bool done = false;
    while (not done) {
        done = true;
        for (uint32_t i = source, j = 1; j <= vertex_num; i = i % vertex_num + 1, ++j) {
            const auto& v = vertex[i];
            if (v.distance_to_source == std::numeric_limits<uint32_t>::max())
                continue;
            for (const auto& u : v.neighbors) {
                uint32_t me = v.distance_to_source + u.weight;
                const uint32_t you = u.to;
                if (me < vertex[you].distance_to_source) {
                    vertex[you].distance_to_source = me;
                    parent[you] = i;
                    done = false;
                }
            }
        }
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
