#include <bits/stdc++.h>
#include <cstdio>
#include <mpi.h>


MPI_Comm graph_communicator;
int32_t vertex_num, edge_num;
int32_t user_threads;
int32_t source;
int32_t Rank, nprocs;
int32_t* parents;
int32_t activeConstructor;

char* input_file;
char* output_file;

std::ofstream outfs;

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

    Vertex() : distance_to_source(std::numeric_limits<uint32_t>::max()) {}
};

void FindPath(int32_t root) {
    std::vector<int32_t> reverse_vec;
    int32_t current = root;
    reverse_vec.emplace_back(current);
    while (current != source) {
        current = parents[current];
        reverse_vec.emplace_back(current);
    }

    std::reverse(reverse_vec.begin(), reverse_vec.end());

    for (size_t i = 0; i < reverse_vec.size() - 1; ++i) {
        outfs << reverse_vec[i] + 1 << " ";
    }
    outfs << reverse_vec[reverse_vec.size() - 1] + 1 << std::endl;
}

void GoToLine(std::ifstream& file, uint32_t num) {
    file.seekg(std::ios::beg);
    for (ssize_t i = 0; i < num - 1; ++ i) {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

void ParallelConstructGraph() {
    int32_t extra = edge_num % activeConstructor;
    int32_t chunk_size = edge_num / activeConstructor + 1;
    int32_t work_load;
    int32_t start_line;

    if (Rank < extra) {
        start_line = chunk_size * Rank + 2;
        work_load = chunk_size;
    } else {
        start_line = (chunk_size * extra) + (chunk_size - 1) * (Rank - extra) + 2;
        work_load = chunk_size - 1;
    }

    std::ifstream infs;
    infs.open(input_file);

    GoToLine(infs, start_line);

    int32_t* sources;
    int32_t* degrees;
    int32_t* destinations;
    int32_t* weights;

    std::vector<Vertex> vertex(vertex_num);

    for (ssize_t i = 0; i < work_load; ++ i) {
        uint32_t start, to, weight;
        infs >> start >> to >> weight;
        -- start; -- to;
        vertex[start].neighbors.emplace_back(to, weight);
        vertex[to].neighbors.emplace_back(start, weight);
    }

    const int32_t double_edge = work_load * 2;
    sources = new int32_t[vertex_num];
    degrees = new int32_t[vertex_num];
    destinations = new int32_t[double_edge];
    weights = new int32_t[double_edge];

    int32_t current = -1;
    for (ssize_t i = 0; i < vertex_num; ++ i) {
        sources[i] = i;
        degrees[i] = vertex[i].neighbors.size();
        for (auto const v : vertex[i].neighbors) {
            destinations[++ current] = v.to;
            weights[current] = v.weight;
        }
    }

    infs.close();

    MPI_Dist_graph_create(MPI_COMM_WORLD, vertex_num, sources, degrees, destinations, weights, MPI_INFO_NULL, false, &graph_communicator);
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    user_threads = std::stoi(argv[1]);
    source = std::stoi(argv[4]) - 1;
    input_file = argv[2];
    output_file = argv[3];
    activeConstructor = (nprocs < 16) ? nprocs : 16;

    if (Rank < activeConstructor) {
        std::ifstream infs;
        infs.open(input_file);
        infs >> vertex_num >> edge_num;
        infs.close();

        auto startTimer = std::chrono::high_resolution_clock::now();
        ParallelConstructGraph();
        auto endTimer = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = endTimer - startTimer;
        std::cout << "Construct Graph Took " << diff.count() << " Second" << std::endl;
    } else {
        int32_t* sources = nullptr;
        int32_t* degrees = nullptr;
        int32_t* destinations = nullptr;
        int32_t* weights = nullptr;
        MPI_Dist_graph_create(MPI_COMM_WORLD, 0, sources, degrees, destinations, weights, MPI_INFO_NULL, false, &graph_communicator);
    }

    MPI_Comm_rank(graph_communicator, &Rank);

    int32_t id, od, we;
    MPI_Dist_graph_neighbors_count(graph_communicator, &id, &od, &we);

    int32_t* clear_source = new int32_t[id];
    int32_t* clear_source_weight = new int32_t[id];
    int32_t* clear_dest = new int32_t[od];
    int32_t* clear_dest_weight = new int32_t[od];
    int32_t* neighbors_dest_to_source = new int32_t[id];

    MPI_Dist_graph_neighbors(graph_communicator, id, clear_source, clear_source_weight, od, clear_dest, clear_dest_weight);

    bool global_done = false;
    int32_t distance_to_source = std::numeric_limits<int32_t>::max();
    int32_t parent = source;

    if ( source == Rank ) {
        distance_to_source = 0;
    }

    while ( not global_done ) {
        MPI_Neighbor_allgather(&distance_to_source, 1, MPI_INT, neighbors_dest_to_source, 1, MPI_INT, graph_communicator);
        bool local_done = true;
        for ( ssize_t i = 0; i < id; ++i) {
            if (neighbors_dest_to_source[i] == std::numeric_limits<int32_t>::max())
                continue;
            const int32_t your_weight = neighbors_dest_to_source[i] + clear_source_weight[i];
            if ( your_weight < distance_to_source) {
                distance_to_source = your_weight;
                parent = clear_source[i];
                local_done = false;
            }
        }
        MPI_Allreduce(&local_done, &global_done, 1, MPI::BOOL, MPI_LAND, graph_communicator);
    }

    if (Rank == 0)
        parents = new int32_t[vertex_num];

    MPI_Gather(&parent, 1, MPI_INT, parents, 1, MPI_INT, 0, graph_communicator);

    if (Rank == 0) {
        outfs.open(output_file);

        for (ssize_t i = 0; i < vertex_num; ++i) {
            if (i == source)
                outfs << source + 1 << " " << source + 1 << std::endl;
            else
                FindPath(i);
        }
        outfs.close();
    }

    MPI_Finalize();
}
