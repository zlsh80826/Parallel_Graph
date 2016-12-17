#include <bits/stdc++.h>
#include <cassert>
#include <cstdio>
#include <mpi.h>
#include <sys/types.h>
#include <unistd.h>

#define WHITE 0
#define BLACK 1

#define TOKEN 80826
#define DISTANCE 30904
#define TERMINTAE 77777

MPI_Comm graph_communicator;
MPI_Datatype MPI_TOKEN;
int32_t vertex_num, edge_num;
int32_t user_threads;
int32_t source;
int32_t Rank, nprocs;
int32_t* parents;
int32_t activeConstructor;

char* input_file;
char* output_file;

std::ofstream outfs;

/*struct Token {
public:
    int32_t color;
    int32_t terminate;

    void init() {
        this -> color = WHITE;
        this -> terminate = false;
    }

    void setTerminate() {
        this -> terminate = true;
    }

    void print() {
        fprintf(stderr,"Rank %d has Token: color = %d, terminate = %d\n", Rank + 1, this -> color, this -> terminate);
    }

    Token() : color(WHITE), terminate(false) {}
} token;*/

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
    for (ssize_t i = 0; i < num - 1; ++i) {
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

    for (ssize_t i = 0; i < work_load; ++i) {
        uint32_t start, to, weight;
        infs >> start >> to >> weight;
        --start;
        --to;
        vertex[start].neighbors.emplace_back(to, weight);
        vertex[to].neighbors.emplace_back(start, weight);
    }

    const int32_t double_edge = work_load * 2;
    sources = new int32_t[vertex_num];
    degrees = new int32_t[vertex_num];
    destinations = new int32_t[double_edge];
    weights = new int32_t[double_edge];

    int32_t current = -1;
    for (ssize_t i = 0; i < vertex_num; ++i) {
        sources[i] = i;
        degrees[i] = vertex[i].neighbors.size();
        for (auto const v : vertex[i].neighbors) {
            destinations[++current] = v.to;
            weights[current] = v.weight;
        }
    }

    infs.close();

    MPI_Dist_graph_create(MPI_COMM_WORLD, vertex_num, sources, degrees, destinations, weights, MPI_INFO_NULL, false,
                          &graph_communicator);
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
    vertex_num = 48;
    if (Rank < activeConstructor) {
        std::ifstream infs;
        infs.open(input_file);
        infs >> vertex_num >> edge_num;
        infs.close();

        auto startTimer = std::chrono::high_resolution_clock::now();
        ParallelConstructGraph();
        auto endTimer = std::chrono::high_resolution_clock::now();
        std::chrono::duration<double> diff = endTimer - startTimer;
        std::cerr << "Construct Graph Took " << diff.count() << " Second" << std::endl;
    } else {
        int32_t* sources = nullptr;
        int32_t* degrees = nullptr;
        int32_t* destinations = nullptr;
        int32_t* weights = nullptr;
        MPI_Dist_graph_create(MPI_COMM_WORLD, 0, sources, degrees, destinations, weights, MPI_INFO_NULL, false,
                              &graph_communicator);
    }

    MPI_Comm_rank(graph_communicator, &Rank);

    int32_t indegree, outdegree, weighted;
    MPI_Dist_graph_neighbors_count(graph_communicator, &indegree, &outdegree, &weighted);

    assert(indegree > 0);
    assert(outdegree > 0);
    int32_t* sources = new int32_t[indegree];
    int32_t* sources_weight = new int32_t[indegree];
    int32_t* destinations = new int32_t[outdegree];
    int32_t* destinations_weight = new int32_t[outdegree];
    // int32_t* neighbors_dest_to_source = new int32_t[indegree];
    MPI_Barrier(graph_communicator);
    MPI_Dist_graph_neighbors(graph_communicator, indegree, sources, sources_weight, outdegree, destinations,
                             destinations_weight);

    fprintf(stderr, "Rank %d has %d indegree %d outdegree\n", Rank + 1, indegree, outdegree);

    // int32_t* cache = new int32_t[vertex_num];
    std::vector<int32_t> cache(vertex_num);
    assert(vertex_num > 0);
    assert(indegree <= vertex_num);
    for (int i = 0; i < vertex_num; ++i) {
        cache.at(i) = -1;
    }
    assert(indegree == outdegree);
    for (int i = 0; i < indegree; ++i) {
        // std::cerr << i << " " << vertex_num << " " << destinations[i] << " " << destinations_weight[i] << std::endl;
        fprintf(stderr, "%d connected to %d weights: %d\n", Rank + 1, destinations[i] + 1, destinations_weight[i]);
        cache.at(destinations[i]) = destinations_weight[i];
    }

    cache.at(Rank) = 0;
    // MPI_Finalize();
    // return 0;

    // fprintf(stderr,"Rank %d done\n", Rank);

    int32_t distance_to_source = std::numeric_limits<int32_t>::max();
    int32_t new_distance_to_source = 0;
    int32_t zero = 0;
    int32_t parent = source;

    bool color = WHITE;
    bool terminate = false;
    int32_t prev_rank = (Rank == 0) ? nprocs - 1 : Rank - 1;
    int32_t next_rank = (Rank == nprocs - 1) ? 0 : Rank + 1;
    int32_t source_prev = (source == 0) ? nprocs - 1 : source - 1;
    MPI_Request requestArr[3];

    MPI_Request& tokenRequest = requestArr[0];
    MPI_Request& distanceRequest = requestArr[1];
    MPI_Request& terminateRequest = requestArr[2];
    MPI_Status distanceStatus;
    MPI_Status tokenStatus;
    MPI_Status terminateStatus;

    int32_t hasWorkload = false;
    int32_t hasToken = false;
    int32_t needTerminate = false;

    int32_t token_with_color = WHITE;
    int32_t token_with_terminate = false;

    MPI_Barrier(graph_communicator);
    std::cerr << "BARRIERRRRRRRRRRRRRRRRRRRRRRR\n";
    MPI_Barrier(graph_communicator);
    MPI_Comm_rank(graph_communicator, &Rank);

    MPI_Irecv(&new_distance_to_source, 1, MPI_INT, MPI_ANY_SOURCE, DISTANCE, graph_communicator, &distanceRequest);
    // MPI_Irecv(&token, 1, MPI_TOKEN, MPI_ANY_SOURCE, TOKEN, graph_communicator, &tokenRequest);
    MPI_Irecv(&token_with_color, 1, MPI_INT, MPI_ANY_SOURCE, TOKEN, graph_communicator, &tokenRequest);
    MPI_Irecv(&token_with_terminate, 1, MPI_INT, MPI_ANY_SOURCE, TERMINTAE, graph_communicator, &terminateRequest);

    MPI_Barrier(graph_communicator);
    fprintf(stderr, "%d %d\n", Rank, getpid());
    MPI_Barrier(graph_communicator);

    if (Rank == source) {
        // token.init();
        // token.color = BLACK;
        token_with_color = BLACK;
        assert(zero == 0);
        MPI_Send(&zero, 1, MPI_INT, source, DISTANCE, graph_communicator);
        std::cerr << "SEND ZERO DONE\n";
        MPI_Send(&token_with_color, 1, MPI_INT, source, TOKEN, graph_communicator);
        std::cerr << "SEND DONE\n";
    }

    MPI_Status SSS;
    std::cerr << "?" << std::endl;
    MPI_Barrier(graph_communicator);

    for (int i = 0; i < 1; ++i) {
        // usleep(500);
        fprintf(stderr, "Rank %d %d\n", Rank, i);
        // MPI_Test(&distanceRequest, &hasWorkload, &distanceStatus);
        int32_t index;
        MPI_Waitany(3, requestArr, &index, &SSS);
        if (index == 1) {
            fprintf(stderr, "Rank %d receive %d cache: %d distance_to_source: %d\n", Rank + 1, new_distance_to_source,
                    cache.at(SSS.MPI_SOURCE), distance_to_source);
            if (new_distance_to_source != std::numeric_limits<int32_t>::max()) {

                int tmp = distance_to_source;
                if (new_distance_to_source + cache.at(SSS.MPI_SOURCE) < distance_to_source) {
                    distance_to_source = new_distance_to_source + cache.at(SSS.MPI_SOURCE);
                    parent = SSS.MPI_SOURCE;
                    fprintf(stderr, "Rank %d was update by %d\nnew distance is %d, old distance is %d, weight is %d\n",
                            Rank + 1, parent + 1, distance_to_source, tmp, cache.at(SSS.MPI_SOURCE));
                    for (ssize_t i = 0; i < outdegree; ++i) {
                        MPI_Send(&distance_to_source, 1, MPI_INT, destinations[i], DISTANCE, graph_communicator);
                        fprintf(stderr, "%d Send to %d\n", Rank, destinations[i]);
                        if (color == WHITE) {
                            if (Rank > source && destinations[i] < Rank && destinations[i] >= source) {
                                color = BLACK;
                            } else if (Rank < source && (destinations[i] < Rank || destinations[i] >= source)) {
                                color = BLACK;
                            }
                        }
                    }
                }
            }
            MPI_Irecv(&new_distance_to_source, 1, MPI_INT, MPI_ANY_SOURCE, DISTANCE, graph_communicator,
                      &distanceRequest);
            // MPI_Test(&distanceRequest, &hasWorkload, &distanceStatus);

        } else if (index == 0) {
            MPI_Irecv(&token_with_color, 1, MPI_INT, MPI_ANY_SOURCE, TOKEN, graph_communicator, &tokenRequest);

        } else {
        }

        /*MPI_Test(&terminateRequest, &needTerminate, &terminateStatus);
        if (needTerminate) {
            if (Rank != source_prev)
                MPI_Send(&terminate, 1, MPI_INT, next_rank, TERMINTAE, graph_communicator);
        }

        MPI_Test(&tokenRequest, &hasToken, &tokenStatus);
        if (hasToken) {
            // token.print();
            fprintf(stderr,"Rank %d receive %d token\n", Rank + 1, token_with_color);
            if (Rank == source) {
                /*if (token.color == WHITE) {
                    token.setTerminate();
                    terminate = true;
                    MPI_Send(&token, 1, MPI_TOKEN, next_rank, TOKEN, graph_communicator);
                } else {
                    token.init();
                    MPI_Send(&token, 1, MPI_TOKEN, next_rank, TOKEN, graph_communicator);
                    MPI_Irecv(&token, 1, MPI_TOKEN, MPI_ANY_SOURCE, TOKEN, graph_communicator, &tokenRequest);
                }*/
        /*if (token_with_color == WHITE) {
            terminate = true;
            token_with_terminate = true;
            MPI_Send(&token_with_terminate, 1, MPI_INT, next_rank, TOKEN, graph_communicator);
        } else {
            token_with_color = WHITE;
            MPI_Send(&token_with_color, 1, MPI_INT, next_rank, TOKEN, graph_communicator);
            MPI_Irecv(&token_with_color, 1, MPI_INT, MPI_ANY_SOURCE, TOKEN, graph_communicator, &tokenRequest);
        }
    } else {
        /*if (token.terminate == true) {
            terminate = true;
            if (Rank != source_prev)
                MPI_Send(&token, 1, MPI_TOKEN, next_rank, TOKEN, graph_communicator);
        } else {
            if (color == BLACK) {
                token.color = BLACK;
                color = WHITE;
            }
            MPI_Send(&token, 1, MPI_TOKEN, next_rank, TOKEN, graph_communicator);
            MPI_Irecv(&token, 1, MPI_TOKEN, MPI_ANY_SOURCE, TOKEN, graph_communicator, &tokenRequest);
        }*/
        /*if (color == BLACK) {
            token_with_color = BLACK;
            color = WHITE;
        }
        MPI_Send(&token_with_color, 1, MPI_INT, next_rank, TOKEN, graph_communicator);
        MPI_Irecv(&token_with_color, 1, MPI_INT, MPI_ANY_SOURCE, TOKEN, graph_communicator, &tokenRequest);
    }
}*/
        // MPI_Barrier(graph_communicator);
    }
    MPI_Barrier(graph_communicator);
    if (Rank == 0) {
        std::cerr << "need to allocate\n";
        parents = new int32_t[vertex_num];
        std::cerr << "allocate done\n";
    }
    MPI_Barrier(graph_communicator);
    fprintf(stderr, "Start Gather\n");
    MPI_Gather(&parent, 1, MPI_INT, parents, 1, MPI_INT, 0, MPI_COMM_WORLD);
    fprintf(stderr, "End Gather\n");
    MPI_Barrier(graph_communicator);
    if (Rank == 0) {
        std::cerr << "Gather done\n";
        /*for (int i = 0; i < vertex_num; ++i)
        {
            std::cerr << parents[i] << std::endl;
        }*/
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
