#include <bits/stdc++.h>

int32_t vertex_num, edge_num;
int32_t* parent;
bool* inSet;

int32_t user_threads;
int32_t source;
char* input_file;
char* output_file;
std::fstream infs;
std::fstream outfs;

class heapNode {
  public:
    int32_t distance_to_source;
    int32_t index;

    bool operator<(const heapNode& h) const { return this->distance_to_source < h.distance_to_source; }

    bool operator>(const heapNode& h) const { return this->distance_to_source > h.distance_to_source; }
};

class vertex {
  public:
    int32_t distance_to_source;
    int32_t index;
    std::vector<std::pair<int32_t, int32_t>> neighbors;

    vertex() : distance_to_source(std::numeric_limits<int32_t>::max()) {}

    void add_neighbor(int32_t index, int32_t weight) { neighbors.emplace_back(index, weight); }

    bool operator<(const vertex& v) const { return this->distance_to_source < v.distance_to_source; }

    bool operator>(const vertex& v) const { return this->distance_to_source > v.distance_to_source; }
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
    input_file = argv[2];
    output_file = argv[3];
    source = std::stoi(argv[4]);

    infs.open(input_file, std::fstream::in);

    infs >> vertex_num >> edge_num;

    std::vector<vertex> graph(vertex_num + 1);
    parent = new int32_t[vertex_num + 1];
    inSet = new bool[vertex_num + 1];

    for (int32_t i = 0; i <= vertex_num; ++i) {
        parent[i] = source;
        graph[i].index = i;
        inSet[i] = false;
    }

    for (int32_t i = 0; i < edge_num; ++i) {
        int32_t start, to, weight;
        infs >> start >> to >> weight;
        graph[start].add_neighbor(to, weight);
        graph[to].add_neighbor(start, weight);
    }

    infs.close();

    std::priority_queue<heapNode, std::vector<heapNode>, std::greater<heapNode>> Q;
    graph[source].distance_to_source = 0;
    parent[source] = source;
    Q.push((heapNode){0, source});

    while (not Q.empty()) {
        auto node = Q.top();
        Q.pop();
        if (inSet[node.index])
            continue;
        inSet[node.index] = true;
        const auto u = graph[node.index];
        for (const auto& v : u.neighbors) {
            if ((not inSet[graph[v.first].index]) &&
                (graph[v.first].distance_to_source > u.distance_to_source + v.second)) {
                graph[v.first].distance_to_source = u.distance_to_source + v.second;
                Q.push((heapNode){graph[v.first].distance_to_source, v.first});
                parent[v.first] = u.index;
            }
        }
    }

    outfs.open(output_file, std::fstream::out);

    for (int32_t i = 1; i <= vertex_num; ++i) {
        if (i == source)
            outfs << source << " " << source << std::endl;
        else
            FindPath(i);
    }

    outfs.close();
}
