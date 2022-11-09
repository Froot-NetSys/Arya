// graph.h
#ifndef GRAPH_H
#define GRAPH_H

#include <iostream>
#include <fstream>
#include <sstream>
#include <vector>
#include <map>
#include <random>
#include <math.h>
#include <algorithm>
#include <string>
#include <unordered_map>
#include <iterator>
#include <assert.h>
#include <sys/socket.h>
#include <arpa/inet.h>
#include <unistd.h>
#include <libmemcached/memcached.h>

using namespace std;

#define FLAG_print false

/* Vertex structure
* (degree, edge_index to edge_list of its first neighbor edge)
*/
struct Machine {
/*    Machine(string node_ip_, uint32_t start_v_id_, uint32_t end_v_id_, uint32_t start_e_id_, uint32_t end_e_id_)
            : node_ip(node_ip_), start_v_id(start_v_id_), end_v_id(end_v_id_), start_e_id(start_e_id_),
              end_e_id(end_e_id_) {
        memcached_server_st *memcached_servers_parse(char *server_strings);
        graph_memc = memcached_create(NULL);
        memcached_return rc;
        graph_servers = memcached_server_list_append(graph_servers, node_ip_.c_str(), 11211, &rc);
        rc = memcached_server_push(graph_memc, graph_servers);
        if (rc != MEMCACHED_SUCCESS)
            cout << "Couldn't add server: " << memcached_strerror(graph_memc, rc) << endl;
    }

    memcached_server_st *graph_servers = NULL;
    memcached_st *graph_memc;
*/

    Machine(string node_ip_, uint32_t start_v_id_, uint32_t end_v_id_, uint32_t start_e_id_, uint32_t end_e_id_)
            : node_ip(node_ip_), start_v_id(start_v_id_), end_v_id(end_v_id_), start_e_id(start_e_id_),
              end_e_id(end_e_id_) {
        struct sockaddr_in serv_addr;

        if ((sock = socket(AF_INET, SOCK_STREAM, 0)) < 0) {
            printf("\n Socket creation error \n");
        }
        serv_addr.sin_family = AF_INET;
        serv_addr.sin_port = htons(11211);
        if (inet_pton(AF_INET, node_ip_.c_str(), &serv_addr.sin_addr) <= 0)
            printf("\nInvalid address/ Address not supported \n");
        if (connect(sock, (struct sockaddr *) &serv_addr, sizeof(serv_addr)) < 0)
            printf("\nConnection Failed \n");
    }


    // useless
    memcached_server_st *graph_servers = NULL;
    memcached_st *graph_memc;

    string node_ip;
    uint32_t start_v_id;
    uint32_t end_v_id;
    uint32_t start_e_id;
    uint32_t end_e_id;

    // structures for communication
    int sock = 0;
    char msg[100000] = {0};
    uint32_t part_len = 0;
    char buffer[5000000] = {0};
};

enum NodeType {
    CycleRoot, CycleLeave, StarRoot, StarLeave
};
string NodeTypeStr[] = {
        "CycleRoot",
        "CycleLeave",
        "StarRoot",
        "StarLeave"
};

/* Edge structure
* (v_start, v_end)
*/
struct Edge {
    Edge(uint32_t v_start_, uint32_t v_end_) : v_start(v_start_), v_end(v_end_) {
    }

    uint32_t v_start;
    uint32_t v_end;
};

/* Vertex structure
* (degree, edge_index to edge_list of its first neighbor edge)
*/
struct Vertex {
    Vertex(uint32_t degree_, uint32_t edge_index_) : degree(degree_), edge_index(edge_index_) {
    }

    uint32_t degree;
    uint32_t edge_index;
};

/* Pattern class
*/

class Pattern {
    /* This class is responsible of
    * 1. storing cycles
    * 2. storing stars
    * 3. storing other edges, missing from cycles and stars
    */

public:

    /* Initialize function
    * 1. store graph file path
    */
    Pattern(string pattern_file_path) : pattern_file_path_(pattern_file_path) {
    }

    /* dump
    * dump graph
    */
    void dump() {
        cout << "Pattern has " << cycle_num << " cycles, and " << star_num << " stars" << endl;

        cout << "cycles: " << endl;
        for (uint32_t i = 0; i < cycles.size(); i++) {
            for (uint32_t j = 0; j < cycles[i].size(); j++)
                cout << cycles[i][j] << " ";
            cout << endl;
        }

        cout << "stars: " << endl;
        for (uint32_t i = 0; i < stars.size(); i++) {
            for (uint32_t j = 0; j < stars[i].size(); j++)
                cout << stars[i][j] << " ";
            cout << endl;
        }

        cout << "other edges: " << endl;
        for (uint32_t i = 0; i < other_edges_.size(); i++) {
            cout << other_edges_[i].v_start << " " << other_edges_[i].v_end << endl;
        }
    }

    /* generate_candidate_maps
    * generate candidate vertex maps for cycles and stars
    */
    void generate_candidate_maps() {
        // cout << "Pattern has " << cycle_num << " cycles, and " << star_num << " stars" << endl;

        for (uint32_t i = 0; i < cycles.size(); i++) {
            vector<vector<uint32_t>> cur_cycle_candidate_map;
            for (uint32_t j = 0; j < cycles[i].size(); j++) {
                uint32_t cycle_len = cycles[i].size();
                // clock
                vector<uint32_t> cur_cycle_clock;
                cur_cycle_clock.push_back(cycles[i][j]);
                for (uint32_t k = (j + 1) % cycle_len; k != j;) {
                    cur_cycle_clock.push_back(cycles[i][k]);
                    k = (k + 1) % cycle_len;
                }
                cur_cycle_candidate_map.push_back(cur_cycle_clock);

                // reverse clock
                /*
                vector<uint32_t> cur_cycle_reverse;
                cur_cycle_reverse.push_back(cycles[i][j]);
                for (uint32_t k = (j - 1 + cycle_len)%cycle_len; k != j; ) {
                    cur_cycle_reverse.push_back(cycles[i][k]);
                    k = (k - 1 + cycle_len) % cycle_len;
                }
                cur_cycle_candidate_map.push_back(cur_cycle_reverse);
                */
            }
            candidate_cycle_v_maps.push_back(cur_cycle_candidate_map);
        }

        /*
        for (uint32_t i = 0; i < stars.size(); i++) {
            vector<vector<uint32_t>> cur_star_candidate_map;

            std::vector<uint32_t> petals;
            for (uint32_t j = 1; j < stars[i].size(); j++) {
                petals.push_back(stars[i][j]);
            }

            do {
                vector<uint32_t> cur_star;
                cur_star.push_back(stars[i][0]);
                for (uint32_t j = 0; j < petals.size(); j++) {
                    cur_star.push_back(petals[j]);
                }
                cur_star_candidate_map.push_back(cur_star);
            } while (next_permutation(petals.begin(), petals.end()));

            candidate_star_v_maps.push_back(cur_star_candidate_map);
        }
         */

        for (uint32_t i = 0; i < stars.size(); i++) {
            vector<vector<uint32_t>> cur_star_candidate_map;
            cur_star_candidate_map.push_back(stars[i]);
            candidate_star_v_maps.push_back(cur_star_candidate_map);
        }


        return;
        // dump all candidates

        for (uint32_t i = 0; i < candidate_cycle_v_maps.size(); i++) {
            cout << "cycle " << i << ": ";
            for (uint32_t j = 0; j < candidate_cycle_v_maps[i].size(); j++) {
                for (uint32_t k = 0; k < candidate_cycle_v_maps[i][j].size(); k++)
                    cout << candidate_cycle_v_maps[i][j][k] << " ";
                cout << ";";
            }
            cout << endl;
        }

        for (uint32_t i = 0; i < candidate_star_v_maps.size(); i++) {
            cout << "star " << i << ": ";
            for (uint32_t j = 0; j < candidate_star_v_maps[i].size(); j++) {
                for (uint32_t k = 0; k < candidate_star_v_maps[i][j].size(); k++)
                    cout << candidate_star_v_maps[i][j][k] << " ";
                cout << ";";
            }
            cout << endl;
        }

    }

    /* load
     * load graph data from graph_file_path to store in data structures:
     * line 1: odd-cycle number star number
     * following lines are odd-cycle and star
     * after that: are other edge
    */
    void load() {
        // cout << "* To load the pattern file!" << endl;

        // Read in num cycles and num of stars
        string line;
        ifstream infile(pattern_file_path_);
        getline(infile, line);
        stringstream statstream(line);
        statstream >> cycle_num >> star_num;

        // Read in cycles (u1, u2, ..., w1)
        for (uint32_t i = 0; i < cycle_num; i++) {
            getline(infile, line);
            stringstream linestream(line);

            vector<uint32_t> cur_cycle;
            uint32_t cur_vertex;
            while (linestream >> cur_vertex) {
                cur_cycle.push_back(cur_vertex);
                vertex_count++;
            }

            cycles.push_back(cur_cycle);
        }

        // Read in stars (center, petal1, petal2, ...)
        for (uint32_t i = 0; i < star_num; i++) {
            getline(infile, line);
            stringstream linestream(line);

            vector<uint32_t> cur_star;
            uint32_t cur_vertex;
            while (linestream >> cur_vertex) {
                cur_star.push_back(cur_vertex);
                vertex_count++;
            }

            stars.push_back(cur_star);
        }


        // read in other edges
        while (getline(infile, line)) {
            stringstream linestream(line);
            uint32_t v_start, v_end;
            linestream >> v_start >> v_end;
            other_edges_.push_back(Edge(v_start, v_end));
        }

        // dump();

        // to test automorphism
        generate_candidate_maps();
    }

    string pattern_file_path_;        // the file that stores pattern data
    uint32_t vertex_count = 0;
    uint32_t cycle_num = 0;
    uint32_t star_num = 0;
    vector<vector<uint32_t>> cycles;
    vector<vector<uint32_t>> stars;

    vector<vector<vector<uint32_t>>>
            candidate_cycle_v_maps;
    vector<vector<vector<uint32_t>>>
            candidate_star_v_maps;

    vector<Edge> other_edges_;      // other edges not covered by stars and cycles
};

/* Graph class
 */
class Graph {
    /* This class is responsible of
    * 1. storing the graph data
    * 2. serve graph queries (degree, neighbor, pair, edge-sampling)
    */

public:

    /* Initialize function
    * 1. store graph file path
    */
    Graph() {
    }

    // helper functions
    /* insert key-value pair into memcached
    */
    void insert_memcached(memcached_st *memc, string key, string value) {
        // cout << "insert_memcached: " << key << " " << value << endl;
        memcached_return rc = memcached_set(memc, key.c_str(), key.length(), value.c_str(), value.length(),
                                            (time_t) 1000000, (uint32_t) 0);
        if (rc != MEMCACHED_SUCCESS)
            cout << "!!!!!! Couldn't store key: " << key << " " << memcached_strerror(memc, rc) << endl;
    }

    /* get neighbors of vertex, via query memcached
    */
    vector<uint32_t> get_neighbors_memcached(uint32_t thread_id, memcached_st *memc, string vertex_str) {
        size_t value_len;
        uint32_t flags;
        memcached_return_t error;
        char *value_str;

        value_str = memcached_get(memc, vertex_str.c_str(), vertex_str.length(), &value_len, &flags, &error);

        if (value_str == NULL) {
            cout << "!!!!!! get neighbors fail: " << thread_id << " " << vertex_str << " "
                 << memcached_strerror(memc, error) << endl;
            throw -1;
            return {};
        } else {
            // cout << "Get value: " << vertex_str << " " << value_str << endl;
            std::stringstream iss(value_str);
            uint32_t next_v;
            std::vector<uint32_t> neighbors;
            while (iss >> next_v)
                neighbors.push_back(next_v);
            free(value_str);
            return neighbors;
        }

    }

    /* get degree of vertex, via query memcached
    */
    uint32_t get_degree_memcached(uint32_t thread_id, memcached_st *memc, string vertex_str) {
        size_t value_len;
        uint32_t flags;
        memcached_return_t error;
        char *value_str;
        value_str = memcached_get(memc, vertex_str.c_str(), vertex_str.length(), &value_len, &flags, &error);

        if (value_str == NULL) {
            cout << "!!!!!! get degree fail: " << thread_id << " " << vertex_str << " "
                 << memcached_strerror(memc, error) << endl;
            throw -1;
            return -1;
        } else {
            // cout << "Get value: " << vertex_str << " " << value_str << endl;
            std::stringstream iss(value_str);
            uint32_t next_v;
            std::vector<uint32_t> neighbors;
            while (iss >> next_v)
                neighbors.push_back(next_v);

            free(value_str);
            return neighbors.size();
        }
    }

    /* get edge according to edge_id
    */
    Edge get_edge_memcached(uint32_t thread_id, memcached_st *memc, string edge_str) {
        size_t value_len;
        uint32_t flags;
        memcached_return_t error;
        char *value_str;

        value_str = memcached_get(memc, edge_str.c_str(), edge_str.length(), &value_len, &flags, &error);

        if (value_str == NULL) {
            cout << "get edge fail: " << thread_id << " " << edge_str << " " << memcached_strerror(memc, error) << endl;
            throw -1;
            return Edge(-1, -1);
        } else {
            // cout << "Get value: " << edge_str << " " << value_str << endl;
            uint32_t v_start, v_end;
            std::istringstream is(value_str);
            is >> v_start >> v_end;
            free(value_str);
            return Edge(v_start, v_end);
        }
    }

    /* load_machine_info
    * load vertices and edges index range from files
    */
    void load_machine_info(string machine_info_file_path, uint32_t thread_id) {
        ifstream infile(machine_info_file_path);
        string line;
        getline(infile, line);
        stringstream linestream(line);
        linestream >> partition_number >> total_edges_;

        // Step 1: set sqrt_m and edge_distribution
        if (thread_id == 0) {
            sqrt_m = sqrt(total_edges_);
            rand_generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
            uniform_int_distribution<uint32_t> new_edge_distribution(0, total_edges_ - 1);
            edge_distribution = new_edge_distribution;
        }

        // Step 2: load each machines (storing graphs) info
        string node_ip;
        uint32_t start_v_id;
        uint32_t end_v_id;
        uint32_t start_e_id;
        uint32_t end_e_id;

        for (uint32_t machine_id = 0; machine_id < partition_number; machine_id++) {
            getline(infile, line);
            stringstream linestream(line);
            linestream >> node_ip >> start_v_id >> end_v_id >> start_e_id >> end_e_id;
            machine_map_.emplace(thread_id * partition_number + machine_id,
                                 Machine(node_ip, start_v_id, end_v_id, start_e_id, end_e_id));
            // if (thread_id == 0)
            //    cout << node_ip << " " << start_v_id << " " << end_v_id << " " << start_e_id << " " << end_e_id << endl;
        }
    }

    // Graph queries

    /* Degree query
    * lookup the degree of Vertex vertex_id
    */
    uint32_t degree(uint32_t thread_id, uint32_t vertex_id) {
        memcached_st *memc = NULL;
        // retrieve the node saving corresponding subgraph
        for (auto it = machine_map_.begin(); it != machine_map_.end(); it++) {
            if (it->first >= thread_id * partition_number && it->first < (thread_id + 1) * partition_number)
                if (it->second.start_v_id <= vertex_id && vertex_id <= it->second.end_v_id) {
                    memc = it->second.graph_memc;
                    break;
                }
        }
        return get_degree_memcached(thread_id, memc, "v-" + to_string(vertex_id));
    }

    /* Neighborhood lookup query
    * retrieve the idx'th neighbor of vertex_id
    */
    uint32_t retrieve_neighbor(uint32_t thread_id, uint32_t vertex_id, uint32_t idx) {
        // return edge_list_[v_table_.at(vertex_id).edge_index + idx].v_end;

        memcached_st *memc = NULL;
        // retrieve the node saving corresponding subgraph
        for (auto it = machine_map_.begin(); it != machine_map_.end(); it++) {
            if (it->first >= thread_id * partition_number && it->first < (thread_id + 1) * partition_number)
                if (it->second.start_v_id <= vertex_id && vertex_id <= it->second.end_v_id) {
                    memc = it->second.graph_memc;
                    break;
                }
        }
        return get_neighbors_memcached(thread_id, memc, "v-" + to_string(vertex_id))[idx];
    }

    /* Adjency test query
    * test whether vertex_1 and vertex_2 is adjacent
    * TODO: this implementation is pretty expensive
    */
    bool neighbor_test(uint32_t thread_id, uint32_t vertex_id_1, uint32_t vertex_id_2) {

        memcached_st *memc = NULL;
        // retrieve the node saving corresponding subgraph
        for (auto it = machine_map_.begin(); it != machine_map_.end(); it++) {
            if (it->first >= thread_id * partition_number && it->first < (thread_id + 1) * partition_number)
                if (it->second.start_v_id <= vertex_id_1 && vertex_id_1 <= it->second.end_v_id) {
                    memc = it->second.graph_memc;
                    break;
                }
        }
        vector<uint32_t> neighbors = get_neighbors_memcached(thread_id, memc, "v-" + to_string(vertex_id_1));

        // iterate edges of vertex_1
        for (auto next_v: neighbors) {
            if (vertex_id_2 == next_v)
                return true;
        }

        return false;
    }

    /* Edge sampling query
    * uniformly random pick a edge from the edge list
    */
    Edge edge_sampling(uint32_t thread_id, default_random_engine & rand_generator) {
        // memcached: get edge_id -> edge
        memcached_st *memc = NULL;
        uint32_t edge_index = edge_distribution(rand_generator);

        for (auto it = machine_map_.begin(); it != machine_map_.end(); it++) {
            if (it->first >= thread_id * partition_number && it->first < (thread_id + 1) * partition_number)
                if (it->second.start_e_id <= edge_index && edge_index <= it->second.end_e_id) {
                    memc = it->second.graph_memc;
                    break;
                }
        }
        return get_edge_memcached(thread_id, memc, "edge-" + to_string(edge_index));
    }

    /* Vertex ordering
    * totally order vertex according to degree
    */
    bool vertex_before(uint32_t thread_id, uint32_t v_1, uint32_t v_2) {
        uint32_t v_1_degree = degree(thread_id, v_1);
        uint32_t v_2_degree = degree(thread_id, v_2);

        if (v_1_degree < v_2_degree) {
            return true;
        } else if (v_1_degree == v_2_degree && v_1 < v_2) {
            return true;
        }
        return false;
    }

    // mapping from machine_id to the partition
    unordered_map<uint32_t, Machine> machine_map_;

    /* Graph Data Structures to store the graph, and support graph queries
    * 1. Degree table: array/hash table for D(v) ; small
    * 2. N(v, i):
    * 3. P(u, v): ????? read in all list of edges, then check is there a v
    * 4. Edge array: array of all edges ; large, support edge-sampling
    * Vertex table: map vertex id -> (degree, start position of edge list)
    * edge array: a edge will be stored in both places (start and end point); not considering directed
    */

    // unordered_map <uint32_t, Vertex> v_table_;
    // vector <Edge> edge_list_;

    /* basic stats
    */
    uint64_t total_edges_ = 0;   // usually denoted as m
    double sqrt_m;

    /* for random sampling
    */
    default_random_engine rand_generator;
    uniform_int_distribution<uint32_t> edge_distribution;


    // Data structures for memcached connection
    uint32_t edge_idx = 0;
    /*
     * In memcached:
     * 1. vertex_id -> (adjacency list)
     * 2. edge_id -> edge
     */

    // for partition
    uint32_t partition_number;


    // structures for communications
    /*
    uint32_t **edge_ids_per_partition = (uint32_t **) malloc(sizeof(uint32_t *) * graph.partition_number);
    uint32_t *idx_per_partition = (uint32_t *) malloc(sizeof(uint32_t) * graph.partition_number);
    uint32_t **center_v_per_partition = (uint32_t **) malloc(sizeof(uint32_t *) * graph.partition_number);
    uint32_t *idx_v_per_partition = (uint32_t *) malloc(sizeof(uint32_t) * graph.partition_number);
    */
};


#endif
