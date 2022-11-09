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

using namespace std;

#define FLAG_print false
enum NodeType {
    CycleRoot, CycleLeave, StarRoot, StarLeave, CycleRootFail, StarRootFail, CycleLeaveFail, StarLeaveFail
};
string NodeTypeStr[] = {
        "CycleRoot",
        "CycleLeave",
        "StarRoot",
        "StarLeave",
        "CycleRootFail",
        "StarRootFail",
        "CycleLeaveFail",
        "StarLeaveFail"
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
                for (uint32_t k = (j + 1)%cycle_len; k != j; ) {
                    cur_cycle_clock.push_back(cycles[i][k]);
                    k = (k+1) % cycle_len;
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



        // dump all candidates

        /*
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
	*/

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

            vector <uint32_t> cur_cycle;
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

            vector <uint32_t> cur_star;
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
    vector <vector<uint32_t>> cycles;
    vector <vector<uint32_t>> stars;

    vector<vector<vector<uint32_t>>> candidate_cycle_v_maps;
    vector<vector<vector<uint32_t>>> candidate_star_v_maps;

    vector <Edge> other_edges_;      // other edges not covered by stars and cycles
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
    Graph(string graph_file_path)
            : graph_file_path_(graph_file_path) {
    }

    // helper functions
    /* print_info
    * print basic info for this graph class
    */
    void print_info() {
        cout << "Graph is at " << graph_file_path_ << endl;
    }

    /* dump vertex table
    * print basic info for this graph class
    */
    void dump_vertex_table() {
        cout << "Dump vertex table " << endl;
        for (auto it = v_table_.begin(); it != v_table_.end(); it++) {
            cout << it->first << " (" << it->second.degree << " " << it->second.edge_index << "): ";
            for (uint32_t idx = 0; idx < it->second.degree; idx++) {
                cout << " " << edge_list_[it->second.edge_index + idx].v_end;
            }
            cout << endl;
        }
    }

    /* enumerate triangles
    * enumerate triangles in the graph
    */
    uint64_t enumerate_triangles() {
        cout << "====== Start to enumerate triangles " << endl;
        uint64_t triangle_count = 0;

        uint32_t first_vertex_id, second_vertex_id, third_vertex_id;

        // enumerate first vertex
        for (auto it = v_table_.begin(); it != v_table_.end(); it++) {
            first_vertex_id = it->first;
            Vertex &first_vertex = it->second;

            // enumerate second vertex, from first_vertex's neighborhoods, and second_vertex > first_vertex
            for (uint32_t idx_1 = 0; idx_1 < first_vertex.degree; idx_1++) {
                second_vertex_id = edge_list_[first_vertex.edge_index + idx_1].v_end;
                if (second_vertex_id <= first_vertex_id)
                    continue;

                Vertex &second_vertex = v_table_.at(second_vertex_id);

                // enumerate third vertex, from second's neighborhoods, and thrid_vertex > second_vertex
                for (uint32_t idx_2 = 0; idx_2 < second_vertex.degree; idx_2++) {
                    third_vertex_id = edge_list_[second_vertex.edge_index + idx_2].v_end;
                    if (third_vertex_id <= second_vertex_id)
                        continue;

                    // test triangle
                    if (neighbor_test(first_vertex_id, third_vertex_id)) {
                        //cout << "!" << first_vertex_id << " " << second_vertex_id << " " << third_vertex_id << endl;
                        triangle_count += 1;
                    } else {
                        //cout << "x" << first_vertex_id << " " << second_vertex_id << " " << third_vertex_id << endl;
                    }
                }

            }

        }


        return triangle_count;
    }


    /* DFS helper function
    * for finding k odd-cycles
    */
    vector<bool> visited;
    uint32_t target_cycle_length;
    uint32_t total_cycles;

    void DFS_cycle(vector <uint32_t> &cycle_vertex_ids) {
        /*
        for (auto it = cycle_vertex_ids.begin(); it != cycle_vertex_ids.end(); it++)
            cout << " " << *it;
        cout << endl;
        */

        // set end cases
        if (cycle_vertex_ids.size() == target_cycle_length) {
            // check cycle
            if (neighbor_test(cycle_vertex_ids.front(), cycle_vertex_ids.back())) {
                // print out
                ///*
                if (FLAG_print) {
                    cout << "find a cycle:";
                    for (auto it = cycle_vertex_ids.begin(); it != cycle_vertex_ids.end(); it++)
                        cout << " " << *it;
                    cout << endl;
                }
                //*/
                // count odd cycle
                total_cycles++;

            }

            return;
        }

        // recursively iterate a new vertex
        uint32_t cur_vertex_id = cycle_vertex_ids.back();
        Vertex &cur_vertex = v_table_.at(cur_vertex_id);
        uint32_t cur_vertex_edge_index = cur_vertex.edge_index;
        for (uint32_t idx = 0; idx < cur_vertex.degree; idx++) {
            uint32_t next_vertex_id = edge_list_[cur_vertex_edge_index + idx].v_end;

            if (!visited[next_vertex_id]) {
                // update states
                cycle_vertex_ids.push_back(next_vertex_id);
                visited[next_vertex_id] = true;

                // iterate
                DFS_cycle(cycle_vertex_ids);

                // reset states
                cycle_vertex_ids.pop_back();
                visited[next_vertex_id] = false;
            }

        }

        return;

    }

    /* enumerate k odd-cycles
    * enumerate k odd-cycles in the graph with DFS
    */
    uint64_t enumerate_odd_cycles(uint32_t cycle_length) {
        cout << "====== Start to enumerate " << cycle_length << "-cycles " << endl;
        visited = vector<bool>(total_vertex_, false);
        target_cycle_length = cycle_length;
        total_cycles = 0;


        // enumerate first vertex
        uint32_t progress = 0;
        for (auto it = v_table_.begin(); it != v_table_.end(); it++) {
            vector <uint32_t> cycle_vertex_ids;
            cycle_vertex_ids.push_back(it->first);
            visited[it->first] = true;

            DFS_cycle(cycle_vertex_ids);
            visited[it->first] = false;

            progress++;
            //if (progress % 10 == 0)
            //cout << "finished check " << progress << " out of "<< total_vertex_ << " vertexes" << endl;
        }

        // remove same cycles; e.g. 1,2,3,4,5 = 3,4,5,1,2 = 1,5,4,3,2
        //return total_cycles / (cycle_length * 2);
        return total_cycles;
    }

    /* enumerate stars
    * enumerate stars in the graph
    */
    /*
    vector<vector<uint32_t>> star_grow_vertex(vector<uint32_t> candidate_ids, uint32_t petal_num) {
    // end cases
    if (petal_num == 0 || candidate_ids.size() < petal_num) {
        return {};
    }

    // common cases
    for (int i = 0; i < candidate_ids.size(); i++) {
        candidate_ids[i] + star_grow_vertex()
    }

        return {};
    }
     */


    uint64_t enumerate_stars(uint32_t petal_num) {
        cout << "====== Start to enumerate " << petal_num << "-stars " << endl;
        uint64_t star_count = 0;
        // iterate center node
        for (auto it = v_table_.begin(); it != v_table_.end(); it++) {
            /*
            // get neighbors
            vector<uint32_t> neighbor_ids;
            Vertex & cur_vertex = v_table_.at(it->first);
            uint32_t cur_vertex_edge_index = cur_vertex.edge_index;
            for (uint32_t idx = 0; idx < cur_vertex.degree; idx++) {
                neighbor_ids.push_back(edge_list_[cur_vertex_edge_index + idx].v_end);
            }

            // pick petal_num neighbors
            vector<vector<uint32_t>> petals = star_grow_vertex(neighbor_ids, petal_num);
            for (int i = 0; i < petals.size(); i++) {
                cout << it->first << " ";
                for (int j = 0; j < petals[i].size(); j++)
                    cout << petals[i][j] << " ";
                cout << endl;
            }

            */

            // get degree
            uint32_t cur_vertex_degree = v_table_.at(it->first).degree;

            // calculate CkN   n!/ (k!* (n-k)!)
            if (petal_num > cur_vertex_degree)
                continue;
            if (petal_num == cur_vertex_degree) {
                star_count += 1;
            } else if (petal_num < cur_vertex_degree) {
                uint64_t this_count = 1;

                for (uint32_t i = cur_vertex_degree; i > cur_vertex_degree - petal_num; i--)
                    this_count = this_count * i;

                for (uint32_t i = 2; i <= petal_num; i++)
                    this_count = this_count / i;

                star_count += this_count;
            }
        }
        return star_count;
    }

    /* load
    * load graph data from graph_file_path to store in data structures:
    *
    */
    void load() {
        // cout << "* To load the graph file!" << endl;

        // Read the file line by line, assuming the edges are already sorted
        std::ifstream infile(graph_file_path_);
        uint32_t v_start, v_end;
        uint32_t cur_v = numeric_limits<uint32_t>::max();
        uint32_t cur_v_degree = 0;
        uint32_t edge_idx = 0;

        while (infile >> v_start >> v_end) {
            // cout << v_start << " " << v_end << endl;
            // if first v_start, update vertex table
            if (v_start != cur_v) {
                v_table_.emplace(v_start, Vertex((uint32_t) 0, edge_idx));

                if (cur_v != numeric_limits<uint32_t>::max()) {
                    // update cur_v degree
                    v_table_.at(cur_v).degree = cur_v_degree;
                }

                cur_v = v_start;
                cur_v_degree = 0;
            }

            // update edge list
            edge_list_.push_back(Edge(v_start, v_end));
            cur_v_degree += 1;
            edge_idx += 1;
        }

        // update last vertex
        v_table_.at(cur_v).degree = cur_v_degree;

        // wrap up some stats
        total_edges_ = edge_list_.size();
        sqrt_m = sqrt(total_edges_);
        total_vertex_ = v_table_.size();
        // cout << "load finished: get total " << total_vertex_ << " vertexes, and " << total_edges_ << " edges" << endl;

        rand_generator.seed(std::chrono::system_clock::now().time_since_epoch().count());
        uniform_int_distribution <uint32_t> new_edge_distribution(0, total_edges_ - 1);
        edge_distribution = new_edge_distribution;
    }

    // Graph queries

    /* Degree query
    * lookup the degree of Vertex vertex_id
    */
    uint32_t degree(uint32_t vertex_id) {
	return v_table_.at(vertex_id).degree;
    }

    /* Neighborhood lookup query
    * retrieve the idx'th neighbor of vertex_id
    */
    uint32_t retrieve_neighbor(uint32_t vertex_id, uint32_t idx) {
        //cout << "to retrieve " << vertex_id << " " << idx;
        return edge_list_[v_table_.at(vertex_id).edge_index + idx].v_end;
    }

    /* Adjency test query
    * test whether vertex_1 and vertex_2 is adjacent
    * TODO: this implementation is pretty expensive
    */
    bool neighbor_test(uint32_t vertex_id_1, uint32_t vertex_id_2) {
        Vertex &vertex_1 = v_table_.at(vertex_id_1);
        uint32_t vertex_1_edge_index = vertex_1.edge_index;

        // iterate edges of vertex_1
        for (uint32_t idx = 0; idx < vertex_1.degree; idx++) {
            if (vertex_id_2 == edge_list_[vertex_1_edge_index + idx].v_end)
                return true;
        }

        return false;
    }

    /* Edge sampling query
    * uniformly random pick a edge from the edge list
    */
    Edge edge_sampling() {
        return edge_list_[edge_distribution(rand_generator)];
    }

    /* Vertex ordering
    * totally order vertex according to degree
    */
    bool vertex_before(uint32_t v_1, uint32_t v_2) {
        uint32_t v_1_degree = degree(v_1);
        uint32_t v_2_degree = degree(v_2);

        if (v_1_degree < v_2_degree) {
            return true;
        } else if (v_1_degree == v_2_degree && v_1 < v_2) {
            return true;
        }
        return false;
    }

//  private:
    string graph_file_path_;        // the file that stores graph data

    /* Graph Data Structures to store the graph, and support graph queries
    * 1. Degree table: array/hash table for D(v) ; small
    * 2. N(v, i):
    * 3. P(u, v): ????? read in all list of edges, then check is there a v
    * 4. Edge array: array of all edges ; large, support edge-sampling
    * Vertex table: map vertex id -> (degree, start position of edge list)
    * edge array: a edge will be stored in both places (start and end point); not considering directed
    */

    unordered_map <uint32_t, Vertex> v_table_;
    vector <Edge> edge_list_;

    /* basic stats
    */
    uint64_t total_vertex_;
    uint64_t total_edges_;   // usually denoted as m
    double sqrt_m;

    /* for random sampling
    */
    default_random_engine rand_generator;
    uniform_int_distribution <uint32_t> edge_distribution;
};


#endif
