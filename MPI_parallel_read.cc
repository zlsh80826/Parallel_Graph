#include <bits/stdc++.h>
#include <mpi.h>

#define GRAPH_DATA 666

uint32_t vertex_num, edge_num;
uint32_t user_threads;
uint32_t source;
int32_t Rank, nprocs, sender_num;

MPI_Datatype MPI_Edge;

char* input_file;
char* output_file;

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

std::vector<Vertex> vertex;
Vertex me;

void GoToLine(std::ifstream& file, uint32_t num) {
    file.seekg(std::ios::beg);
    for (ssize_t i = 0; i < num - 1; ++ i) {
        file.ignore(std::numeric_limits<std::streamsize>::max(), '\n');
    }
}

void recv_neighbors() {
    for ( ssize_t i = 0; i < sender_num; ++ i) {
        std::vector<Edge> mirror;
        MPI_Request request;
        MPI_Irecv(mirror.data(), vertex_num, MPI_Edge, MPI_ANY_SOURCE, GRAPH_DATA, MPI_COMM_WORLD, &request);
    }
}

void send_neighbors() {

}

void parallelReadFile(int total_read_rank) {

    std::ifstream infs;
    infs.open(input_file);
    uint32_t extra = edge_num % total_read_rank;
    uint32_t chunk_size = edge_num / total_read_rank + 1;
    uint32_t work_load;
    uint32_t start_line;

    if ((uint32_t)Rank < extra) {
        start_line = chunk_size * Rank + 2;
        work_load = chunk_size;
    } else {
        start_line = (chunk_size * extra) + (chunk_size - 1) * (Rank - extra) + 2;
        work_load = chunk_size - 1;
    }

    GoToLine(infs, start_line);

    for (ssize_t i = 0; i < work_load; ++ i) {
        uint32_t start, to, weight;
        infs >> start >> to >> weight;
        // std::cout << start << " " << to << " " << weight << std::endl;
        vertex[start].neighbors.emplace_back(to, weight);
        vertex[to].neighbors.emplace_back(start, weight);
    }

    infs.close();
}

int main(int argc, char** argv) {

    MPI_Init(&argc, &argv);
    MPI_Comm_rank(MPI_COMM_WORLD, &Rank);
    MPI_Comm_size(MPI_COMM_WORLD, &nprocs);

    MPI_Datatype type[2] = {MPI_UNSIGNED, MPI_UNSIGNED};
    int blocklen[2] = {1, 1};
    MPI_Aint disp[2];
    disp[0] = 0;
    disp[1] = 4;
    MPI_Type_create_struct(2, blocklen, disp, type, &MPI_Edge);
    MPI_Type_commit(&MPI_Edge);

    double startTime = MPI_Wtime();

    user_threads = std::stoi(argv[1]);
    source = std::stoi(argv[4]);
    input_file = argv[2];
    output_file = argv[3];

    std::ifstream infs;
    infs.open(input_file);

    infs >> vertex_num >> edge_num;

    infs.close();

    sender_num = (nprocs < 16) ? nprocs : 16;

    if (Rank < sender_num) {
        vertex = std::vector<Vertex>(vertex_num + 1);
        parallelReadFile(16);
    }

    MPI_Barrier(MPI_COMM_WORLD);

    double endTime = MPI_Wtime();

    std::cout << "Rank: " << Rank << " took " << endTime - startTime << std::endl;
    MPI_Finalize();
}
