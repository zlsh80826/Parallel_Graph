#include <bits/stdc++.h>

int64_t vertex_num, edge_num;
int64_t* parent;

int64_t user_threads;
int64_t source;
char* input_file;
char* output_file;
std::fstream infs;
std::fstream outfs;

class vertex {
  public:
    int64_t distance_to_source;
    int64_t index;
    bool inSet;
    std::vector<std::pair<int64_t, int64_t>> neighbors;

    vertex() : distance_to_source(std::numeric_limits<int64_t>::max()), inSet(false) {}

    void add_neighbor(int64_t index, int64_t weight) { neighbors.emplace_back(index, weight); }

    bool operator<(const vertex& v) const { return this->distance_to_source < v.distance_to_source; }
};

void FindPath(int64_t root) {
    std::vector<int64_t> reverse_vec;
    int64_t current = root;
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
    source = atoi(argv[4]);

    infs.open(input_file, std::fstream::in);

    infs >> vertex_num >> edge_num;

    std::vector<vertex> graph(vertex_num + 1);
    parent = new int64_t[vertex_num + 1];

    for (int64_t i = 0; i <= vertex_num; ++i) {
        parent[i] = source;
        graph[i].index = i;
    }

    for (int64_t i = 0; i < edge_num; ++i) {
        int64_t start, to, weight;
        infs >> start >> to >> weight;
        graph[start].add_neighbor(to, weight);
        graph[to].add_neighbor(start, weight);
    }

    infs.close();

    std::priority_queue<vertex> Q;
    graph[source].distance_to_source = 0;
    parent[source] = source;
    Q.push(graph[source]);

    while (not Q.empty()) {
        auto u = Q.top();
        Q.pop();
        u.inSet = true;
        for (auto v : u.neighbors) {
            if ((not graph[v.first].inSet) && (graph[v.first].distance_to_source > u.distance_to_source + v.second)) {
                graph[v.first].distance_to_source = u.distance_to_source + v.second;
                Q.push(graph[v.first]);
                parent[v.first] = u.index;
            }
        }
    }

    outfs.open(output_file, std::fstream::out);

    for (int64_t i = 1; i <= vertex_num; ++i) {
    	if ( i == source )
    		outfs << source << " " << source << std::endl;
    	else
        	FindPath(i);
    }

    outfs.close();
}
