#include <bits/stdc++.h>
#include <cstdio>
#include <mpi.h>

uint32_t vertex_num, edge_num;
uint32_t user_threads;
uint32_t source;
int32_t Rank, nprocs;
int32_t* parents;

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

void FindPath(uint32_t root) {
    std::vector<uint32_t> reverse_vec;
    uint32_t current = root;
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

int main(int argc, char** argv) {

    // auto TstartTimer = std::chrono::high_resolution_clock::now();
    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    user_threads = std::stoi(argv[1]);
    source = std::stoi(argv[4]) - 1;
    input_file = argv[2];
    output_file = argv[3];

    MPI_Comm graph_communicator;
    // MPI_Info info;
    int32_t* sources;
    int32_t* degrees;
    int32_t* destinations;
    int32_t* weights;

    // auto TendTimer = std::chrono::high_resolution_clock::now();
    // std::chrono::duration<double> Tdiff = TendTimer - TstartTimer;
    // std::cout << "Construct Graph Took " << Tdiff.count() << " Second" << std::endl;

    if (Rank == 0) {
        auto IstartTimer = std::chrono::high_resolution_clock::now();
        std::ifstream infs;
        infs.open(input_file);

        infs >> vertex_num >> edge_num;
        std::cout << vertex_num << " " << edge_num << std::endl;
        std::vector<Vertex> vertex(vertex_num);
        for (size_t i = 0; i < edge_num; ++ i) {
            uint32_t start, to, weight;
            infs >> start >> to >> weight;
            -- start;
            -- to;
            // std::cout << start << " " << to << " " << weight << std::endl;
            vertex.at(start).neighbors.emplace_back(to, weight);
            vertex.at(to).neighbors.emplace_back(start, weight);
        }

        infs.close();

        const int32_t double_edge = edge_num * 2;
        sources = new int32_t[vertex_num];
        // int32_t* sources_weight = new int32_t[vertex_num];
        degrees = new int32_t[vertex_num];
        destinations = new int32_t[double_edge];
        weights = new int32_t[double_edge];

        int32_t current = -1;
        for (size_t i = 0; i < vertex_num; ++ i) {
            sources[i] = i;
            // sources_weight[i] = std::numeric_limits<int32_t>::max();
            degrees[i] = vertex[i].neighbors.size();
            for (auto const v : vertex[i].neighbors) {
                destinations[++ current] = v.to;
                weights[current] = v.weight;
            }
        }
        auto IendTimer = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> Idiff = IendTimer - IstartTimer;
        std::cout << "Reading Input " << Idiff.count() << " Second" << std::endl;

        auto startTimer = std::chrono::high_resolution_clock::now();
        MPI_Dist_graph_create(MPI_COMM_WORLD, vertex_num, sources, degrees, destinations, weights, MPI_INFO_NULL, false, &graph_communicator);
        auto endTimer = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = endTimer - startTimer;
        std::cout << "Construct Graph Took " << diff.count() << " Second" << std::endl;
    } else {
        sources = nullptr;
        degrees = nullptr;
        destinations = nullptr;
        weights = nullptr;
        MPI_Dist_graph_create(MPI_COMM_WORLD, 0, sources, degrees, destinations, weights, MPI_INFO_NULL, false, &graph_communicator);
    }

    int32_t status;
    MPI_Topo_test(graph_communicator, &status);
    // std::cout << status << " " << MPI_DIST_GRAPH << std::endl;
    MPI_Comm_rank(graph_communicator, &Rank);
    // std::cout << old_rank << " " << Rank << std::endl;
    int32_t id, od, we;
    MPI_Dist_graph_neighbors_count(graph_communicator, &id, &od, &we);
    // printf("Rank: %d has %d indegree %d outdegree %d weights\n", Rank + 1, id, od, we);

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
            const int32_t you_weight = neighbors_dest_to_source[i] + clear_dest_weight[i];
            if ( you_weight < distance_to_source) {
                distance_to_source = you_weight;
                parent = clear_dest[i];
                local_done = false;
            }
        }
        MPI_Allreduce(&local_done, &global_done, 1, MPI::BOOL, MPI_LAND, MPI_COMM_WORLD);
    }

    if (Rank == 0)
        parents = new int32_t[vertex_num];

    MPI_Gather(&parent, 1, MPI_INT, parents, 1, MPI_INT, 0, MPI_COMM_WORLD);

    if (Rank == 0) {
        outfs.open(output_file);

        for (uint32_t i = 0; i < vertex_num; ++i) {
            if (i == source)
                outfs << source + 1 << " " << source + 1 << std::endl;
            else
                FindPath(i);
        }
        outfs.close();
    }

    MPI_Finalize();
}
