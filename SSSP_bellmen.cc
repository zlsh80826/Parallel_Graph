#include <bits/stdc++.h>

int32_t vertex_num, edge_num;
int32_t user_threads;
int32_t source;
int32_t* parent;
char* input_file;
char* output_file;
std::ifstream infs;
std::ofstream outfs;

class Edge {
  public:
    int32_t to;
    int32_t weight;

    Edge(int32_t to, int32_t weight) : to(to), weight(weight) {}

    bool operator<(const Edge& e) const { return this->weight < e.weight; }
};

class Vertex {
  public:
    int32_t distance_to_source;
    std::vector<Edge> neighbors;

    Vertex() : distance_to_source(std::numeric_limits<int32_t>::max()) {}
};

void FindPath(int32_t root) {
    std::vector<int32_t> reverse_vec;
    int32_t current = root;
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

int main(int argc, char** argv) {

    user_threads = std::stoi(argv[1]);
    source = std::stoi(argv[4]);
    input_file = argv[2];
    output_file = argv[3];

    auto IstartTimer = std::chrono::high_resolution_clock::now();

    infs.open(input_file);

    infs >> vertex_num >> edge_num;

    parent = new int32_t[vertex_num + 1];
    for (int i = 1; i <= vertex_num; ++i) {
        parent[i] = source;
    }

    std::vector<Vertex> vertex(vertex_num + 1);

    int32_t start, to, weight;
    for (int32_t i = 0; i < edge_num; ++i) {
        infs >> start >> to >> weight;
        vertex[start].neighbors.emplace_back(to, weight);
        vertex[to].neighbors.emplace_back(start, weight);
    }

    infs.close();

    auto IendTimer = std::chrono::high_resolution_clock::now();
    std::chrono::duration<double> Idiff = IendTimer - IstartTimer;
    std::cout << "Input took " << Idiff.count() << " Second" << std::endl;

    auto startTimer = std::chrono::high_resolution_clock::now();
    vertex[source].distance_to_source = 0;

    int32_t me;
    int32_t you;

    bool done = false;
    while (not done) {
        done = true;
        for (int32_t i = source, j = 1; j <= vertex_num; i = i % vertex_num + 1, ++j) {
            const auto& v = vertex[i];
            if (v.distance_to_source == std::numeric_limits<int32_t>::max())
                continue;
            for (const auto& u : v.neighbors) {
                me = v.distance_to_source + u.weight;
                you = vertex[u.to].distance_to_source;
                if (me < you) {
                    vertex[u.to].distance_to_source = me;
                    parent[u.to] = i;
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

    for (int32_t i = 1; i <= vertex_num; ++i) {
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
