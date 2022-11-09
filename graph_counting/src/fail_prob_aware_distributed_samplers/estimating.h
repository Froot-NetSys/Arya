// estimating.h
#ifndef ESTIMATING_H
#define ESTIMATING_H

#include <iostream>
#include "graph.h"
#include <cmath>
#include <set>
#include <queue>
#include <map>
#include <mutex>
#include <algorithm>
#include "assert.h"

using namespace std;

// related to batch sampling and optimizations
uint32_t sampling_batch = 100;
bool early_pruning = true;


/* SamplingTree strcuture
 * tree of nodes (vertex, probability)
*/
struct SamplingTreeNode {
    SamplingTreeNode(NodeType type_, double probability_) : node_type(type_), probability(probability_) {
        fast_fail = false;
    }

    SamplingTreeNode(NodeType type_, vector<uint32_t> &vertexes_, double probability_) : node_type(type_),
                                                                                         probability(probability_) {
        for (uint32_t i = 0; i < vertexes_.size(); i++) {
            vertexes.push_back(vertexes_[i]);
        }
        fast_fail = false;
    }

    void print() {
        cout << "Node(" << NodeTypeStr[node_type] << ", " << children.size() << ", " << probability << ", [ ";
        for (auto it = vertexes.begin(); it != vertexes.end(); it++)
            cout << *it << " ";
        cout << "]) ";
    }

    NodeType node_type;   // CYCLE_ROOT, CYCLE_LEAF, STAR_ROOT, STAR_LEAF
    vector<struct SamplingTreeNode *> children;
    vector<uint32_t> vertexes;
    double probability;
    string waiting_data;
    bool fast_fail;
    uint32_t d_u0 = 0;
    uint32_t d_v0 = 0;
};

/* PerThread structure
 * storing buffers etc
*/
class SamplingThread {
public:
    unordered_map<string, vector<uint32_t>> retrieved;
    uint32_t msg_vec[1000] = {0};
    uint32_t msg_vec_backup[1000] = {0};
};

SamplingThread thread_struct[100];     // Up to 100 threads

class succeed_conditions {
public:
    double probability;
    vector<pair<uint32_t, uint32_t>> neighbors_need;
};

std::mutex print_mutex;                                    // to cleanly print something
void batch_get_data(Graph &graph, SamplingThread &this_thread, char type, uint32_t msg_len, uint32_t *ids) {
    // for (uint32_t i = 0; i < msg_len; i++)
    //     cout << ids[i] << " ";
    // cout << endl;

    if (msg_len == 0)
        return;
    //step 1: create messages for each partition

    char *msg;   //  a pointer to any possible msg
    for (uint32_t i = 0; i < graph.partition_number; i++) {
        msg = graph.machine_map_.at(i).msg;
        msg[0] = 'g';
        msg[1] = 'e';
        msg[2] = 't';
        graph.machine_map_.at(i).part_len = 3;
    }

    uint32_t tmp_id;
    if (type == 'v') {
        for (uint32_t i = 0; i < msg_len; i++) {
            tmp_id = ids[i];
            // update vertex logic (mod)
            uint32_t tmp_p_id = tmp_id % graph.partition_number;
            msg = graph.machine_map_.at(tmp_p_id).msg;   // why [] cause error?
            uint32_t &part_len = graph.machine_map_.at(tmp_p_id).part_len;

            msg[part_len++] = ' ';
            msg[part_len++] = 'v';
            msg[part_len++] = '-';
            uint32_t digits = 1;
            while (tmp_id / digits > 9) digits *= 10;
            while (true) {
                msg[part_len++] = (tmp_id / digits) % 10 + '0';
                if (digits == 1) break;
                digits /= 10;
            }

/*
            for (auto it = graph.machine_map_.begin(); it != graph.machine_map_.end(); it++) {
                msg = it->second.msg;
                uint32_t &part_len = it->second.part_len;
                if (it->second.start_v_id <= tmp_id && tmp_id <= it->second.end_v_id) {
                    msg[part_len++] = ' ';
                    msg[part_len++] = 'v';
                    msg[part_len++] = '-';
                    uint32_t digits = 1;

                    while (tmp_id / digits > 9) digits *= 10;
                    while (true) {
                        msg[part_len++] = (tmp_id / digits) % 10 + '0';
                        if (digits == 1) break;
                        digits /= 10;
                    }
                    break;
                }
            }
*/
        }
    } else {
        for (uint32_t i = 0; i < msg_len; i++) {
            tmp_id = ids[i];
            for (auto it = graph.machine_map_.begin(); it != graph.machine_map_.end(); it++) {
                msg = it->second.msg;
                uint32_t &part_len = it->second.part_len;
                if (it->second.start_e_id <= tmp_id && tmp_id <= it->second.end_e_id) {
                    msg[part_len++] = ' ';
                    msg[part_len++] = 'e';
                    msg[part_len++] = 'd';
                    msg[part_len++] = 'g';
                    msg[part_len++] = 'e';
                    msg[part_len++] = '-';
                    uint32_t digits = 1;

                    while (tmp_id / digits > 9) digits *= 10;
                    while (true) {
                        msg[part_len++] = (tmp_id / digits) % 10 + '0';
                        if (digits == 1) break;
                        digits /= 10;
                    }
                    break;
                }
            }
        }
    }

    // step 2: send out messages for each partition
    for (uint32_t i = 0; i < graph.partition_number; i++) {
        msg = graph.machine_map_.at(i).msg;
        uint32_t &part_len = graph.machine_map_.at(i).part_len;
        if (part_len <= 3)
            continue;
        msg[part_len++] = '\r';
        msg[part_len++] = '\n';
        msg[part_len++] = '\0';
        //{
        //    std::lock_guard<std::mutex> iolock(print_mutex);
        //    cout << "print out the message: " << msg;
        //}

        int this_sock = graph.machine_map_.at(i).sock;
        send(this_sock, msg, strlen(msg), 0);
    }

    // step 3: read back message ; parse each key-value ; insert into hash map
    int inserted = 0;
    for (uint32_t i = 0; i < graph.partition_number; i++) {
        if (graph.machine_map_.at(i).part_len <= 3)
            continue;
        // read the message
        char *buffer = graph.machine_map_.at(i).buffer;
        int this_sock = graph.machine_map_.at(i).sock;

        uint32_t buffer_idx = 0;
        for (;;) {
            int valread = read(this_sock, buffer + buffer_idx, 1024 * 64);
            buffer_idx += valread;

            if (valread <= 0) {
                cout << valread << endl;
                cout << "error!!!!!!\n";
                exit(1);
            }

            if (buffer_idx > 5000000) {
                cout << "overflow!!!!!\n";
                exit(1);
            }
            if (buffer[buffer_idx - 5] == 'E' && buffer[buffer_idx - 4] == 'N' && buffer[buffer_idx - 3] == 'D') {
                buffer[buffer_idx] = '\0';
                break;
            }
        }

        // parse the return values; insert to (key, vector of values)
        stringstream ss(buffer);
        string to_key, to_val, item;
        int idx = 0;
        //{
        //    std::lock_guard<std::mutex> iolock(print_mutex);
        //    cout << "received message: " << buffer;
        //}
        while (getline(ss, to_key, '\n')) {
            if (to_key[0] == 'V') {
                // parse key
                item = "";
                idx = 6;
                while (to_key[idx] != ' ') {
                    item += to_key[idx];
                    idx += 1;
                }

                // parse value
                getline(ss, to_val, '\n');
                // cout << "get a pair: " << item << " " << to_val << endl;

                uint32_t parse_idx = 0, number = 0;
                uint32_t str_len = to_val.length() - 1;
                vector<uint32_t> neighbors;
                while (parse_idx < str_len) {
                    number = 0;
                    while (to_val[parse_idx] != ' ' && parse_idx < str_len) {
                        number = number * 10 + (to_val[parse_idx] - '0');
                        parse_idx++;
                    }
                    parse_idx++;
                    neighbors.push_back(number);
                }

                // insert into hash table
                this_thread.retrieved[item] = neighbors;
                inserted += 1;
            }
        }
    }

    if (inserted < msg_len) {
        cout << "got is less than queried " << inserted << " " << msg_len << endl;

        for (uint32_t i = 0; i < graph.partition_number; i++) {
            if (graph.machine_map_.at(i).part_len <= 3)
                continue;
            // read the message
            char *buffer = graph.machine_map_.at(i).buffer;
            char * msg = graph.machine_map_.at(i).msg;

            printf("partition %d: %s \n %s\n", i, msg, buffer);
        }
    }

    return;
}

/* Sample a odd-cycle
* from graph, <u1, u2, ..., w1>
* node([(v1, v2), (v3, v4), ...]) -> k x w1
*/
SamplingTreeNode *odd_cycle_sampler(uint32_t thread_id, Graph &graph, uint32_t cycle_length, default_random_engine & rand_generator) {
    //cout << "====== Start to sample a " << cycle_length << "-cycles" << endl;
    // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());
    uint32_t ue_id, d_ue;

    // step 1: sample cycle_length/2 edges (u1, v1), (x, x), (x, x), ... u1 < v1
    vector<uint32_t> sample_edge_vector;
    Edge e1 = graph.edge_sampling(thread_id, rand_generator);
    if (graph.vertex_before(thread_id, e1.v_start, e1.v_end)) {
        sample_edge_vector.push_back(e1.v_start);
        sample_edge_vector.push_back(e1.v_end);
    } else {
        sample_edge_vector.push_back(e1.v_end);
        sample_edge_vector.push_back(e1.v_start);
    }

    for (uint32_t j = 0 + 1; j < cycle_length / 2; j++) {
        Edge sampled_edge = graph.edge_sampling(thread_id, rand_generator);
        sample_edge_vector.push_back(sampled_edge.v_start);
        sample_edge_vector.push_back(sampled_edge.v_end);
    }

    SamplingTreeNode *root = new SamplingTreeNode(CycleRoot,
                                                  sample_edge_vector,
                                                  pow(graph.total_edges_, (cycle_length / 2)) / 2.0);
    assert(root != nullptr);

    // step 2: sample w1, a neighbor of u1
    ue_id = sample_edge_vector[0];
    d_ue = graph.degree(thread_id, ue_id);

    // cout << "cycle u1 degree: " << d_ue << endl;

    uniform_int_distribution<uint32_t> next_edge_distribution(0, d_ue - 1);

    for (uint32_t j = 0; j < ceil(d_ue / graph.sqrt_m); j++) {
        // step 2: sample an vertex from u1's neighbor
        vector<uint32_t> sampled_vertex_id;
        sampled_vertex_id.push_back(graph.retrieve_neighbor(thread_id, ue_id, next_edge_distribution(rand_generator)));
        root->children.push_back(new SamplingTreeNode(CycleLeave, sampled_vertex_id, d_ue));
    }

    return root;
}

/* Sample a star
* from graph, <center, w1, w2, ...>
*/
SamplingTreeNode *star_sampler(uint32_t thread_id, Graph &graph, uint32_t petal_num, default_random_engine & rand_generator) {
    //cout << "====== Start to sample an " << petal_num << "-star" << endl;
    // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());

    // Step 1: choose a vertex, proportional to its degree, choose an edge, its first vertex
    Edge e1 = graph.edge_sampling(thread_id, rand_generator);
    uint32_t center_vertex_id = e1.v_start;
    uint32_t center_vertex_degree = graph.degree(thread_id, center_vertex_id);
    vector<uint32_t> center_vertex;
    center_vertex.push_back(center_vertex_id);

    SamplingTreeNode *root = new SamplingTreeNode(StarRoot,
                                                  center_vertex,
                                                  1.0 * graph.total_edges_ / center_vertex_degree);

    // Step 2: sample petal_num vertexes from N(X), without replacement
    set<uint32_t> vertex_set;
    uniform_int_distribution<uint32_t> next_edge_distribution(0, center_vertex_degree - 1);

    vector<uint32_t> sampled_petals;
    if (center_vertex_degree >= petal_num) {
        // indicate can find enough petals
        for (uint32_t i = 0; i < petal_num; i++) {
            uint32_t rand_index = next_edge_distribution(rand_generator);
            while (vertex_set.find(rand_index) != vertex_set.end()) {
                rand_index = next_edge_distribution(rand_generator);
            }
            vertex_set.insert(rand_index);
            uint32_t sampled_vertex_id = graph.retrieve_neighbor(thread_id, center_vertex_id, rand_index);
            sampled_petals.push_back(sampled_vertex_id);
        }
    }

    // cout << "star vertex degree: " << center_vertex_degree << " , petal num: " << petal_num << endl;
    uint64_t this_count = 1;
    if (petal_num < center_vertex_degree) {
        for (uint32_t i = center_vertex_degree; i > center_vertex_degree - petal_num; i--)
            this_count = this_count * i;

        for (uint32_t i = 2; i <= petal_num; i++)
            this_count = this_count / i;
    }

    root->children.push_back(new SamplingTreeNode(StarLeave, sampled_petals, 1.0 * this_count));

    // cout << center_vertex_id << " " << center_vertex_degree
    //        << " " << root->probability << " " <<  root->children[0]->probability
    //        << " " << root->probability * root->children[0]->probability << endl;

    return root;
}

/* BFS Iterate sampling trees to calculate a value
*/
void printTree(SamplingTreeNode *root) {
    // data structure
    queue<SamplingTreeNode *> to_visit_list;
    int depth = 0;

    // insert root
    to_visit_list.push(root);

    // BFS to print
    while (to_visit_list.size() > 0) {
        cout << "level " << depth << ": ";
        int sz = to_visit_list.size();
        for (int i = 0; i < sz; i++) {
            // pop the front
            SamplingTreeNode *next_print = to_visit_list.front();
            if (next_print == nullptr) {
                cout << "!!!!!! get a nullptr" << endl;
            }
            next_print->print();

            // expand next level
            for (uint32_t j = 0; j < next_print->children.size(); j++) {
                to_visit_list.push(next_print->children[j]);
            }

            to_visit_list.pop();
        }
        cout << endl;
        depth++;
    }
    return;
}

/* Create sampling trees
*/
SamplingTreeNode *construct_sampling_tree(uint32_t thread_id, Graph &graph, Pattern &pattern, default_random_engine & rand_generator) {
    SamplingTreeNode *root = nullptr;
    vector<SamplingTreeNode *> leaves;
    vector<SamplingTreeNode *> new_leaves;

    // Step 1: iterate odd-cycles to add 2-levels
    // cout << "grow tree with cycles" << endl;
    for (uint32_t i = 0; i < pattern.cycle_num; i++) {
        if (i == 0) {
            root = odd_cycle_sampler(thread_id, graph, pattern.cycles[i].size(), rand_generator);
            for (uint32_t j = 0; j < root->children.size(); j++) {
                leaves.push_back(root->children[j]);
            }
        } else {
            // extend every leaf by a cycle tree
            for (uint32_t j = 0; j < leaves.size(); j++) {
                leaves[j]->children.push_back(odd_cycle_sampler(thread_id, graph, pattern.cycles[i].size(), rand_generator));
            }

            // update leaves and new_leaves
            new_leaves.clear();
            for (uint32_t j = 0; j < leaves.size(); j++) {
                SamplingTreeNode *last_level_root = leaves[j]->children[0];
                for (uint32_t k = 0; k < last_level_root->children.size(); k++) {
                    new_leaves.push_back(last_level_root->children[k]);
                }
            }
            leaves = new_leaves;
        }
    }

    // Step 2: iterate stars to add 2-levels
    // cout << "grow tree with stars" << endl;
    for (uint32_t i = 0; i < pattern.star_num; i++) {
        if (root == nullptr) {
            root = star_sampler(thread_id, graph, pattern.stars[i].size() - 1, rand_generator);
            for (uint32_t j = 0; j < root->children.size(); j++) {
                leaves.push_back(root->children[j]);
            }
        } else {
            // extend every leaf by a star tree
            for (uint32_t j = 0; j < leaves.size(); j++) {
                leaves[j]->children.push_back(star_sampler(thread_id, graph, pattern.stars[i].size() - 1, rand_generator));
            }

            // update leaves and new_leaves
            new_leaves.clear();
            for (uint32_t j = 0; j < leaves.size(); j++) {
                SamplingTreeNode *last_level_root = leaves[j]->children[0];
                for (uint32_t k = 0; k < last_level_root->children.size(); k++) {
                    new_leaves.push_back(last_level_root->children[k]);
                }
            }
            leaves = new_leaves;
        }
    }

    return root;
}

/* Testing whether profiles -> actual copy of patterns
* count probabilities as well
*/
bool test_circle_copy(uint32_t thread_id, Graph &graph, vector<uint32_t> &e_vertexes, uint32_t w) {
    // test u1 is the smallest vertex
    uint32_t u1 = e_vertexes[0];
    uint32_t v1 = e_vertexes[1];
    for (uint32_t i = 1; i < e_vertexes.size(); i++) {
        if (!graph.vertex_before(thread_id, u1, e_vertexes[i]))
            return false;
    }

    if (!graph.vertex_before(thread_id, u1, w))
        return false;

    // test v1 < w
    if (!graph.vertex_before(thread_id, v1, w))
        return false;

    // test edges of (v1,u2), ... (vk, w)
    if (!graph.neighbor_test(thread_id, e_vertexes.back(), w))
        return false;

    for (uint32_t i = 1; i < e_vertexes.size() - 1; i += 2) {
        if (!graph.neighbor_test(thread_id, e_vertexes[i], e_vertexes[i + 1]))
            return false;
    }

    return true;
}

/* Major helper function to test a sampled path
*/
void DFS_judge_remaining_edges(uint32_t thread_id, uint32_t idx, Graph &graph, Pattern &pattern,
                               vector<vector<uint32_t>> &cycle_vertexes, vector<vector<uint32_t>> &star_vertexes,
                               map<uint32_t, uint32_t> &v_map, uint32_t &found_instances) {
    // corner cases
    if (idx == cycle_vertexes.size() + star_vertexes.size()) {

        // judge remaining edges
        bool test_succeed = true;
        for (uint32_t i = 0; i < pattern.other_edges_.size(); i++) {
            uint32_t test_v_start = v_map.at(pattern.other_edges_[i].v_start);
            uint32_t test_v_end = v_map.at(pattern.other_edges_[i].v_end);

            if (!graph.neighbor_test(thread_id, test_v_start, test_v_end)) {
                test_succeed = false;
            }
        }

        if (test_succeed) {
            found_instances += 1;
        }
        return;
    }

    // normal grow mapping
    if (idx < cycle_vertexes.size()) {
        for (uint32_t i = 0; i < pattern.candidate_cycle_v_maps[idx].size(); i++) {
            // map candidate_cycle_v_maps[idx][i] to cycle_vertexes[idx]
            for (uint32_t j = 0; j < pattern.candidate_cycle_v_maps[idx][i].size(); j++) {
                v_map.emplace(pattern.candidate_cycle_v_maps[idx][i][j], cycle_vertexes[idx][j]);
            }

            // DFS_judge_further
            DFS_judge_remaining_edges(thread_id, idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map,
                                      found_instances);

            // restore
            for (uint32_t j = 0; j < pattern.candidate_cycle_v_maps[idx][i].size(); j++) {
                v_map.erase(pattern.candidate_cycle_v_maps[idx][i][j]);
            }
        }
    } else {
        uint32_t star_idx = idx - cycle_vertexes.size();
        for (uint32_t i = 0; i < pattern.candidate_star_v_maps[star_idx].size(); i++) {
            // map candidate_star_v_maps[idx][i] to star_vertexes[idx]
            for (uint32_t j = 0; j < pattern.candidate_star_v_maps[star_idx][i].size(); j++) {
                v_map.emplace(pattern.candidate_star_v_maps[star_idx][i][j], star_vertexes[star_idx][j]);
            }

            // DFS_judge_further
            DFS_judge_remaining_edges(thread_id, idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map,
                                      found_instances);

            // restore
            for (uint32_t j = 0; j < pattern.candidate_star_v_maps[star_idx][i].size(); j++) {
                v_map.erase(pattern.candidate_star_v_maps[star_idx][i][j]);
            }
        }
    }
}


uint32_t judge_remaining_edges(uint32_t thread_id, Graph &graph, Pattern &pattern,
                               vector<vector<uint32_t>> &cycle_vertexes, vector<vector<uint32_t>> &star_vertexes) {

    if (pattern.other_edges_.size() == 0) {
        return true;
    }

    uint32_t found_instances = 0;
    map<uint32_t, uint32_t> v_map;
    DFS_judge_remaining_edges(thread_id, 0, graph, pattern, cycle_vertexes, star_vertexes, v_map, found_instances);

    return found_instances;
}

/* Major helper function to iterate a sampling tree
*/
void tree_DFS_helper(uint32_t thread_id, SamplingTreeNode *root, Graph &graph, Pattern &pattern,
                     vector<SamplingTreeNode *> &path, map<uint32_t, uint32_t> &v_map, double &this_cross_prob_inverse,
                     uint64_t &num_cross_partition, double &this_intra_prob_inverse, uint64_t &num_intra_partition) {
    // corner cases
    if (root->children.size() == 0) {
        path.push_back(root);
        v_map.clear();

        // task 1: print path
        assert(path.size() == (pattern.cycle_num + pattern.star_num) * 2);

        if (FLAG_print) {
            cout << "found a path: ";
        }

        set<uint32_t> vertex_set;
        bool test_succeed = true;
        double prob_inverse = 1.0;
        // total_sampling_times += 1;

        // task 2: judge cycles
        vector<vector<uint32_t>> cycle_vertexes;
        for (uint32_t i = 0; i < pattern.cycle_num; i++) {
            set<uint32_t> cycle_vertex_set;
            SamplingTreeNode *cur_root = path[i * 2];
            SamplingTreeNode *cur_leave = path[i * 2 + 1];

            if (FLAG_print) {
                cur_root->print();
                cout << ", ";
                cur_leave->print();
                cout << ", ";
            }

            // handle cycle root
            prob_inverse *= cur_root->probability;

            // map this node's vertexes and cycles[i]'s first several nodes
            vector<uint32_t> cur_cycle;
            for (uint32_t j = 0; j < cur_root->vertexes.size(); j++) {
                v_map.emplace(pattern.cycles[i][j], cur_root->vertexes[j]);
                vertex_set.insert(cur_root->vertexes[j]);
                cur_cycle.push_back(cur_root->vertexes[j]);
                cycle_vertex_set.insert(cur_root->vertexes[j]);
            }

            // handle cycle leave
            prob_inverse *= cur_leave->probability;

            v_map.emplace(pattern.cycles[i].back(), cur_leave->vertexes[0]);
            vertex_set.insert(cur_leave->vertexes[0]);
            cur_cycle.push_back(cur_leave->vertexes[0]);
            cycle_vertexes.push_back(cur_cycle);
            cycle_vertex_set.insert(cur_leave->vertexes[0]);

            // judge cycle construct
            // test overlapped vertex
            if (cycle_vertex_set.size() != pattern.cycles[i].size()) {
                test_succeed = false;
                if (FLAG_print)
                    cout << "cycle " << i << " has overlapped vertexes" << endl;
            }

            // test copy
            if (!test_circle_copy(thread_id, graph, cur_root->vertexes, cur_leave->vertexes[0])) {
                test_succeed = false;
                if (FLAG_print)
                    cout << "cycle " << i << " does not pass copy test" << endl;
            }
        }

        // task 3: judge stars
        vector<vector<uint32_t>> star_vertexes;
        if (test_succeed != false) {
            for (uint32_t i = 0; i < pattern.star_num; i++) {
                SamplingTreeNode *cur_root = path[(pattern.cycle_num + i) * 2];
                SamplingTreeNode *cur_leave = path[(pattern.cycle_num + i) * 2 + 1];

                if (FLAG_print) {
                    cur_root->print();
                    cout << ", ";
                    cur_leave->print();
                    cout << ", ";
                }

                // handle star root
                prob_inverse *= cur_root->probability;
                v_map.emplace(pattern.stars[i].front(), cur_root->vertexes[0]);
                vector<uint32_t> cur_star;
                cur_star.push_back(cur_root->vertexes[0]);
                vertex_set.insert(cur_root->vertexes[0]);

                // handle leaves and and test copy
                prob_inverse *= cur_leave->probability;
                uint32_t petal_num = cur_leave->vertexes.size();
                if (petal_num == 0) {
                    test_succeed = false;
                    if (FLAG_print)
                        cout << "star " << i << " does not have enough petals" << endl;
                } /* else if (petal_num == 1 &&
                           !graph.vertex_before(thread_id, cur_root->vertexes[0], cur_leave->vertexes[0])) {
                    test_succeed = false;
                    if (FLAG_print)
                        cout << "star " << i << " does not form a copy" << endl;
                } */
                else {
                    for (uint32_t j = 0; j < cur_leave->vertexes.size(); j++) {
                        v_map.emplace(pattern.stars[i][j + 1], cur_leave->vertexes[j]);
                        cur_star.push_back(cur_leave->vertexes[j]);
                        vertex_set.insert(cur_leave->vertexes[j]);
                    }
                    star_vertexes.push_back(cur_star);
                }
            }
            if (FLAG_print)
                cout << endl;
        }

        // print v_map
        // assert(v_map.size() == pattern.vertex_count);

        // task 4: judge other edges
        if (test_succeed == false) {
            if (FLAG_print)
                cout << "previous cycle/star tests already failed" << endl;
        } else if (vertex_set.size() != pattern.vertex_count) {
            if (FLAG_print)
                cout << "this path has overlapped sampled vertexes" << endl;
            test_succeed = false;
        } else {
            // judge other edges
            //if (!judge_remaining_edges(graph, pattern, cycle_vertexes, star_vertexes)) {
            //    test_succeed = false;
            //}

            uint32_t num_matches = judge_remaining_edges(thread_id, graph, pattern, cycle_vertexes, star_vertexes);
            if (num_matches == 0) {
                test_succeed = false;
            } else {
                prob_inverse *= num_matches;
                // cout << num_matches << " " << prob_inverse << endl;
            }

        }


        // for distributed version, judge whether the pattern is inter/intra partition
        bool is_cross_partition = false;
        /*
        set<uint32_t> partition_set;
        uint32_t partition_id;

        for (auto vertex = vertex_set.begin(); vertex != vertex_set.end(); vertex++) {
            partition_id = -1;
            for (auto it = graph.machine_map_.begin(); it != graph.machine_map_.end(); it++) {
                if (it->first >= thread_id * graph.partition_number &&
                    it->first < (thread_id + 1) * graph.partition_number)
                    if (it->second.start_v_id <= *vertex && *vertex <= it->second.end_v_id) {
                        partition_id = it->first;
                        break;
                    }
            }
            if (!partition_set.empty() && partition_set.find(partition_id) == partition_set.end()) {
                is_cross_partition = true;
                break;
            }
            partition_set.insert(partition_id);
        }
         */

        if (is_cross_partition)
            num_cross_partition += 1;
        else
            num_intra_partition += 1;

        // task 4: if pass test, count
        if (test_succeed) {
            if (FLAG_print) {
                cout << "the path passed the test!" << endl;
            }
            // total_prob_inverse += prob_inverse;
            if (is_cross_partition)
                this_cross_prob_inverse += prob_inverse;
            else
                this_intra_prob_inverse += prob_inverse;
        }

        path.pop_back();
        return;
    }

    // DFS to iterate children
    path.push_back(root);
    for (uint32_t i = 0; i < root->children.size(); i++) {
        tree_DFS_helper(thread_id, root->children[i], graph, pattern, path, v_map, this_cross_prob_inverse,
                        num_cross_partition, this_intra_prob_inverse, num_intra_partition);
    }
    path.pop_back();

    return;
}

/* Count a sampling tree
*/
void sampling_tree_count(uint32_t thread_id, SamplingTreeNode *root, Graph &graph, Pattern &pattern,
                         double &this_cross_prob_inverse, uint64_t &num_cross_partition,
                         double &this_intra_prob_inverse, uint64_t &num_intra_partition) {
    if (FLAG_print)
        cout << "====== start to count this tree" << endl;

    // DFS iterate from root to a leave
    vector<SamplingTreeNode *> path;
    map<uint32_t, uint32_t> v_map;

    tree_DFS_helper(thread_id, root, graph, pattern, path, v_map, this_cross_prob_inverse, num_cross_partition,
                    this_intra_prob_inverse, num_intra_partition);

    return;
}

void free_sampling_tree(SamplingTreeNode *root) {
    if (root->children.size() == 0) {
        delete root;
        root = nullptr;
        return;
    }

    for (uint32_t i = 0; i < root->children.size(); i++) {
        free_sampling_tree(root->children[i]);
    }

    delete root;
    root = nullptr;
}

/* Estimating arbitrary patterns
* from graph, with sampling_times
*/
tuple<double, uint64_t, double, uint64_t>
estimate_pattern(uint32_t thread_id, Graph &graph, Pattern &pattern, uint32_t sampling_times, default_random_engine & rand_generator) {
    //cout << "= Start to estimate patterns, " << sampling_times << " times" << endl;
    double this_cross_prob_inverse = 0.0;
    double this_intra_prob_inverse = 0.0;
    uint64_t num_cross_partition = 0;
    uint64_t num_intra_partition = 0;

    /*
    SamplingTreeNode *shared_root = new SamplingTreeNode;
    assert(shared_root != nullptr);
    */

    for (uint32_t i = 0; i < sampling_times; i++) {
        try {
            // Step 1: construct the trees
            SamplingTreeNode *root = construct_sampling_tree(thread_id, graph, pattern, rand_generator);

            // Step 2: iterate every path from root to leave
            sampling_tree_count(thread_id, root, graph, pattern, this_cross_prob_inverse, num_cross_partition,
                                this_intra_prob_inverse, num_intra_partition);

            // TODO: need to free tree
            free_sampling_tree(root);
        }
        catch (...) {
            i--;
        }
    }
    return make_tuple(this_cross_prob_inverse, num_cross_partition, this_intra_prob_inverse, num_intra_partition);
}


/* Estimating arbitrary patterns
* from graph, with sampling_times
*/
void parse(char *str_to_parse) {
    stringstream ss(str_to_parse);
    string to_key, to_val, item;
    int idx = 0;

    while (getline(ss, to_key, '\n')) {
        if (to_key[0] == 'E')
            break;

        if (to_key[0] == 'V') {
            item = "";
            idx = 6;
            while (to_key[idx] != ' ') {
                item += to_key[idx];
                idx += 1;
            }

            getline(ss, to_val, '\n');
        }
    }
}

tuple<double, uint64_t, double, uint64_t>
estimate_stars(uint32_t thread_id, Graph &graph, Pattern &pattern, uint32_t sampling_times, default_random_engine & rand_generator) {
    //cout << "= Start to estimate patterns, " << sampling_times << " times" << endl;
    double this_cross_prob_inverse = 0.0;
    double this_intra_prob_inverse = 0.0;
    uint64_t num_cross_partition = 0;
    uint64_t num_intra_partition = 0;

    // option 1: sample one by one
    bool batch_sampling = true;
    if (!batch_sampling) {
        for (uint32_t i = 0; i < sampling_times; i++) {
            try {
                // Step 1: construct the trees
                SamplingTreeNode *root = construct_sampling_tree(thread_id, graph, pattern, rand_generator);

                // Step 2: iterate every path from root to leave
                sampling_tree_count(thread_id, root, graph, pattern, this_cross_prob_inverse, num_cross_partition,
                                    this_intra_prob_inverse, num_intra_partition);

                // TODO: need to free tree
                free_sampling_tree(root);
            }
            catch (...) {
                i--;
            }
        }
        return make_tuple(this_cross_prob_inverse, num_cross_partition, this_intra_prob_inverse, num_intra_partition);
    } else {
        // assumes it's a N-star
        // Step 1: decide random edge ids
        uint32_t **edge_ids_per_partition = (uint32_t **) malloc(sizeof(uint32_t *) * graph.partition_number);
        uint32_t *idx_per_partition = (uint32_t *) malloc(sizeof(uint32_t) * graph.partition_number);
        uint32_t **center_v_per_partition = (uint32_t **) malloc(sizeof(uint32_t *) * graph.partition_number);
        uint32_t *idx_v_per_partition = (uint32_t *) malloc(sizeof(uint32_t) * graph.partition_number);
        double total_prob_verse = 0;
        // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());

        for (uint32_t i = 0; i < graph.partition_number; i++) {
            idx_per_partition[i] = 0; // just a 100 x int
            edge_ids_per_partition[i] = (uint32_t *) malloc(
                    sizeof(uint32_t) * sampling_times);  // just a sampling_times x 100 array
            idx_v_per_partition[i] = 0;  // 100 x int
            center_v_per_partition[i] = (uint32_t *) malloc(
                    sizeof(uint32_t) * sampling_times);  // sampling_times x 100 array
        }

        uint32_t tmp_id;
        for (uint32_t i = 0; i < sampling_times; i++) {
            // tmp_id = graph.edge_distribution(graph.rand_generator);
            tmp_id = graph.edge_distribution(rand_generator);
            for (auto it = graph.machine_map_.begin(); it != graph.machine_map_.end(); it++) {
                if (it->second.start_e_id <= tmp_id && tmp_id <= it->second.end_e_id) {
                    edge_ids_per_partition[it->first][idx_per_partition[it->first]] = tmp_id;
                    idx_per_partition[it->first]++;
                    break;
                }
            }
        }

        // Step 2: batch and send to get random edges
        for (uint32_t i = 0; i < graph.partition_number; i++) {
            string query_str = "";
            char *msg = graph.machine_map_.at(i).msg;
            // char msg[100000] = {0};

            query_str = "get";
            for (uint32_t j = 0; j < idx_per_partition[i]; j++)
                query_str += " edge-" + to_string(edge_ids_per_partition[i][j]);

            query_str += "\r\n\0";
            strcpy(msg, query_str.c_str());
            // cout << "message string: " << msg << endl;

            // Step 2: send
            int this_sock = graph.machine_map_.at(i).sock;
            send(this_sock, msg, strlen(msg), 0);
        }

        // read the return message and creat the new message
        for (uint32_t i = 0; i < graph.partition_number; i++) {
            char *buffer = graph.machine_map_.at(i).buffer;
            // char buffer[1000000] = {0};
            // memset(buffer ,0 , 1000000);
            int this_sock = graph.machine_map_.at(i).sock;

            uint32_t buffer_idx = 0;
            for (;;) {
                int valread = read(this_sock, buffer + buffer_idx, 1024 * 64);
                buffer_idx += valread;

                // cout << buffer_idx << " : " << buffer[buffer_idx-5] << buffer[buffer_idx-4] << buffer[buffer_idx-3] << endl;
                // cout << buffer << endl;

                if (buffer_idx > 1000000) {
                    cout << "overflow!!!!!\n";
                    return make_tuple(0, 0, 0, 0);
                }
                if (buffer[buffer_idx - 5] == 'E') {
                    buffer[buffer_idx] = '\0';
                    break;
                }
            }


            // parse the return edges
            stringstream ss(buffer);
            string to_key, to_val, item;
            int idx = 0;

            while (getline(ss, to_key, '\n')) {
                if (to_key[0] == 'V') {
                    item = "";
                    idx = 6;
                    while (to_key[idx] != ' ') {
                        item += to_key[idx];
                        idx += 1;
                    }

                    getline(ss, to_val, '\n');


                    // parse the edge
                    // stringstream is(to_val);   // this was the bottleneck!!!!!
                    // is >> start_v >> end_v;

                    uint32_t start_v = 0, end_v = 0;
                    uint32_t parse_idx = 0;
                    uint32_t str_len = to_val.length() - 1;
                    while (to_val[parse_idx] != ' ') {
                        start_v = start_v * 10 + (to_val[parse_idx] - '0');
                        parse_idx++;
                    }
                    parse_idx++;
                    while (parse_idx < str_len) {
                        end_v = end_v * 10 + (to_val[parse_idx] - '0');
                        parse_idx++;
                    }
                    // cout << item << " : " << start_v << " " << end_v << endl;
                    // update vertex logic (mod)
                    int tmp_p_id = start_v % graph.partition_number;
                    center_v_per_partition[tmp_p_id][idx_v_per_partition[tmp_p_id]] = start_v;
                    idx_v_per_partition[tmp_p_id]++;

/*
                    for (auto it = graph.machine_map_.begin(); it != graph.machine_map_.end(); it++) {
                        if (it->second.start_v_id <= start_v && start_v <= it->second.end_v_id) {
                            center_v_per_partition[it->first][idx_v_per_partition[it->first]] = start_v;
                            idx_v_per_partition[it->first]++;
                            break;
                        }
                    }
*/

                }
            }
        }


        for (uint32_t i = 0; i < graph.partition_number; i++) {
            // send the v message
            string query_str = "get";
            char *msg = graph.machine_map_.at(i).msg;
            // char msg[100000] = {0};

            for (uint32_t j = 0; j < idx_v_per_partition[i]; j++)
                query_str += " v-" + to_string(center_v_per_partition[i][j]);

            query_str += "\r\n";
            strcpy(msg, query_str.c_str());

            // cout << "message string: " << msg << endl;

            // Step 2: send
            int this_sock = graph.machine_map_.at(i).sock;
            send(this_sock, msg, strlen(msg), 0);
        }


        for (uint32_t i = 0; i < graph.partition_number; i++) {
            char *buffer = graph.machine_map_.at(i).buffer;
            // char buffer[1000000] = {0};
            // memset(buffer ,0 , 1000000);
            int this_sock = graph.machine_map_.at(i).sock;

            uint32_t buffer_idx = 0;
            for (;;) {
                int valread = read(this_sock, buffer + buffer_idx, 1024 * 64);
                buffer_idx += valread;

                // cout << "get here! with " << valread << endl;
                // cout << buffer_idx << " : " << buffer[buffer_idx-5] << buffer[buffer_idx-4] << buffer[buffer_idx-3] << endl;
                // cout << buffer << endl;

                if (buffer_idx > 1000000) {
                    cout << "overflow!!!!!\n";
                    return make_tuple(0, 0, 0, 0);
                }
                if (buffer[buffer_idx - 5] == 'E') {
                    buffer[buffer_idx] = '\0';
                    break;
                }
            }



            // cout << "valread: " << valread << " " << buffer << "\n=========== " << endl;
            // cout << buffer << "\n============= " << endl;
            // parse the return edges
            stringstream ss(buffer);
            string to_key, to_val, item;
            int idx = 0;
            while (getline(ss, to_key, '\n')) {
                // if (to_key[0] == 'E') {
                //    cout << "really get to the vertex end!!!!!!!!!!!\n\n\n\n\n\n" << endl;
                // }

                if (to_key[0] == 'V') {
                    item = "";
                    idx = 6;
                    while (to_key[idx] != ' ') {
                        item += to_key[idx];
                        idx += 1;
                    }

                    getline(ss, to_val, '\n');
                    // continue;

                    // parse the neighbors
                    uint32_t parse_idx = 0, number = 0;
                    uint32_t str_len = to_val.length() - 1;
                    vector<uint32_t> neighbors;
                    while (parse_idx < str_len) {
                        number = 0;
                        while (to_val[parse_idx] != ' ' && parse_idx < str_len) {
                            number = number * 10 + (to_val[parse_idx] - '0');
                            parse_idx++;
                        }
                        parse_idx++;
                        neighbors.push_back(number);
                    }
                    /*
                     std::istringstream is(to_val);
                     vector<uint32_t> neighbors;
                     uint32_t number;
                     while (is >> number)
                         neighbors.push_back(number);
                     */
                    uint32_t center_vertex_degree = neighbors.size();

                    /*
                    cout << item << " " << neighbors.size() << " : ";
                    for (auto it : neighbors)
                         cout << it << " ";
                    cout << endl;
                    */


                    uint32_t petal_num = 2;
                    double root_prob = 1.0 * graph.total_edges_ / center_vertex_degree;
                    // Step 2: sample petal_num vertexes from N(X), without replacement
                    /*
                    set<uint32_t> vertex_set;
                    uniform_int_distribution<uint32_t> next_edge_distribution(0, center_vertex_degree - 1);
                    vector<uint32_t> sampled_petals;
                    if (center_vertex_degree >= petal_num) {
                        // indicate can find enough petals
                        for (uint32_t i = 0; i < petal_num; i++) {
                            uint32_t rand_index = next_edge_distribution(rand_generator);
                            while (vertex_set.find(rand_index) != vertex_set.end()) {
                                rand_index = next_edge_distribution(rand_generator);
                            }
                            vertex_set.insert(rand_index);
                            uint32_t sampled_vertex_id = neighbors[rand_index];
                            sampled_petals.push_back(sampled_vertex_id);
                        }
                    }
                    */
                    uint64_t this_count = 1;
                    if (petal_num < center_vertex_degree) {
                        for (uint32_t i = center_vertex_degree; i > center_vertex_degree - petal_num; i--)
                            this_count = this_count * i;

                        for (uint32_t i = 2; i <= petal_num; i++)
                            this_count = this_count / i;
                    }

                    double child_prob = 1.0 * this_count;
                    total_prob_verse += root_prob * child_prob;
                }
            }

        }
        return make_tuple(total_prob_verse, sampling_times, 0, 0);
    }


}


tuple<double, uint64_t, double, uint64_t>
estimate_odd_cycles(uint32_t thread_id, Graph &graph, Pattern &pattern, uint32_t sampling_times, default_random_engine & rand_generator) {
    cout << "\n\n========== Start to estimate " << sampling_times << " times  ==========\n";
    double total_prob_verse = 0;
    uint32_t num_estimators = 0;

    uint32_t N_odd_cycle = 5;
    // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());


    SamplingThread &this_thread = thread_struct[thread_id];
    SamplingTreeNode *estimators[sampling_times];
    for (uint32_t i = 0; i < sampling_times; i++)
        estimators[i] = new SamplingTreeNode(CycleRoot, pow(graph.total_edges_, (N_odd_cycle / 2)) / 2.0);

    // TODO: possibly move these to the thread structure
    uint32_t msg_vec[1000] = {0};
    uint32_t msg_vec_backup[1000] = {0};
    uint32_t msg_idx = 0;

    // Step 1: Handle u1,v1 and related testing ; get u1, v1; decide u1, v1 ; choose w; get w; test u1 < w and v1 < w
    // Step 1.1 get u1, v1 for every estimator
    cout << "========== get u0, v0 ==========\n";
    uint32_t tmp_id;
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        // tmp_id = graph.edge_distribution(graph.rand_generator);
        tmp_id = graph.edge_distribution(rand_generator);
        estimators[estimator_idx]->waiting_data = "edge-" + to_string(tmp_id);
        msg_vec[estimator_idx] = tmp_id;
    }

    batch_get_data(graph, this_thread, 'e', sampling_times,
                   msg_vec);  // send out the list of edges, and put them in the hash table


    // Step 1.2 get u1, v1 neighbors -> create two arrays store u1, v1 neighbors
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        try {
            vector<uint32_t> &edge = this_thread.retrieved.at(this_estimator->waiting_data);

            this_estimator->vertexes.push_back(edge[0]);
            this_estimator->vertexes.push_back(edge[1]);
            msg_vec[estimator_idx] = edge[0];
            msg_vec_backup[estimator_idx] = edge[1];
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            cout << "cannot find!!!! 1\n";
            continue;
        }
    }
    batch_get_data(graph, this_thread, 'v', sampling_times,
                   msg_vec);  // send out the list of edges, and put them in the hash table
    batch_get_data(graph, this_thread, 'v', sampling_times,
                   msg_vec_backup);  // send out the list of edges, and put them in the hash table

    // Step 1.3 compare u1, v1, decide u1, v1, w -> get w neighbors, and create leave SamplingNode
    cout << "========== compare u0, v0, decide w ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        if (this_estimator->fast_fail)
            continue;

        uint32_t v0_id = this_estimator->vertexes[0];
        uint32_t v1_id = this_estimator->vertexes[1];

        try {
            vector<uint32_t> &neighbors_0 = this_thread.retrieved.at("v-" + to_string(v0_id));
            vector<uint32_t> &neighbors_1 = this_thread.retrieved.at("v-" + to_string(v1_id));

            uint32_t d0 = neighbors_0.size();
            uint32_t d1 = neighbors_1.size();

            // decide u1, v1
            if ((d0 > d1) || (d0 == d1 && v0_id > v1_id)) {
                // cout << "want to swap ";
                //<< v0_id << " : " << v1_id << " : ";
                this_estimator->d_u0 = d1;
                this_estimator->d_v0 = d0;
                this_estimator->vertexes[0] = v1_id;
                this_estimator->vertexes[1] = v0_id;

                // choose w   // TODO create as many leave SamplingNodes as in the algorithm
                cout << this_estimator->vertexes[0] << " vs. " << this_estimator->vertexes[1] << " : "
                     << this_estimator->d_u0 << " " << this_estimator->d_v0;
                uniform_int_distribution<uint32_t> next_edge_distribution(0, d1 - 1);
                // for (uint32_t j = 0; j < ceil(d0 / graph.sqrt_m); j++) {
                vector<uint32_t> sampled_vertex_id;
                sampled_vertex_id.push_back(neighbors_1[next_edge_distribution(rand_generator)]);

                cout << " sampled w: " << sampled_vertex_id[0] << endl;

                this_estimator->children.push_back(new SamplingTreeNode(CycleLeave, sampled_vertex_id, d1));
                msg_vec[estimator_idx] = sampled_vertex_id[0]; // TODO note this should be changed if multiple sampled w
                // }
            } else {
                this_estimator->d_u0 = d0;
                this_estimator->d_v0 = d1;

                // choose w   // TODO create as many leave SamplingNodes as in the algorithm
                cout << this_estimator->vertexes[0] << " vs. " << this_estimator->vertexes[1] << " : "
                     << this_estimator->d_u0 << " " << this_estimator->d_v0;
                uniform_int_distribution<uint32_t> next_edge_distribution(0, d0 - 1);
                // for (uint32_t j = 0; j < ceil(d0 / graph.sqrt_m); j++) {
                vector<uint32_t> sampled_vertex_id;
                sampled_vertex_id.push_back(neighbors_0[next_edge_distribution(rand_generator)]);

                cout << " sampled w: " << sampled_vertex_id[0] << endl;

                this_estimator->children.push_back(new SamplingTreeNode(CycleLeave, sampled_vertex_id, d0));
                msg_vec[estimator_idx] = sampled_vertex_id[0]; // TODO note this should be changed if multiple sampled w
                // }
            }
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            cout << "cannot find!!!! 2\n";
            continue;
        }

    }
    batch_get_data(graph, this_thread, 'v', sampling_times, msg_vec);  // TODO note this should be changed if multiple sampled w

    // Step 1.4 test u1 < w and v1 < w
    cout << "========== test ordering of w and u0,v0 ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        if (this_estimator->fast_fail)
            continue;

        // test ordering
        bool fail_root = true;
        for (uint32_t j = 0; j < this_estimator->children.size(); j++) {
            SamplingTreeNode *this_child = this_estimator->children[j];
            try {
                uint32_t w_id = this_child->vertexes[0];
                uint32_t w_degree = this_thread.retrieved.at("v-" + to_string(w_id)).size();
                // test v0 < w only (no need to test u0 < w because u0 < v0 already)
                cout << "w: " << w_id << ", " << w_degree << " vs v0: " << this_estimator->vertexes[1] << ", "
                     << this_estimator->d_v0 << " ";
                if ((w_degree < this_estimator->d_v0) ||
                    (w_degree == this_estimator->d_v0 && w_id <= this_estimator->vertexes[1])) {
                    this_child->fast_fail = true;
                    cout << "ordering fail!";
                } else {
                    fail_root = false;  // one child is not failed, then cannot fail the root
                }
                cout << endl;
            } catch (const std::out_of_range &oor) {
                this_child->fast_fail = true;
                cout << "cannot find!!!! 3\n";
                continue;
            }
        }

        if (fail_root)
            this_estimator->fast_fail = true;
    }


    // Step 2 get remaining edges, and test ordering and neighboring
    for (uint32_t edge_idx = 1; edge_idx < N_odd_cycle / 2; edge_idx++) {
        cout << "========== get another edge " << edge_idx << " ==========\n";
        // Step 2.1  get a random edge
        uint32_t tmp_id;
        msg_idx = 0;
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // tmp_id = graph.edge_distribution(graph.rand_generator);
            tmp_id = graph.edge_distribution(rand_generator);
            this_estimator->waiting_data = "edge-" + to_string(tmp_id);
            msg_vec[msg_idx++] = tmp_id;
        }
        batch_get_data(graph, this_thread, 'e', msg_idx,
                       msg_vec);  // send out the list of edges, and put them in the hash table

        // Step 2.2 store ui, vi; test vi-1 and ui neighboring; if pass, get ui neighbors
        cout << "========== test vi-1, ui neighboring ==========\n";
        msg_idx = 0;
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // parse the edge
            try {
                vector<uint32_t> &edge = this_thread.retrieved.at(this_estimator->waiting_data);
                this_estimator->vertexes.push_back(edge[0]);
                this_estimator->vertexes.push_back(edge[1]);

                // test vi-1 and ui neighboring
                cout << "neighboring test: vi-1: " << this_estimator->vertexes[edge_idx * 2 - 1] << ", ui: " << edge[0];
                vector<uint32_t> &last_v_neighbors = this_thread.retrieved.at(
                        "v-" + to_string(this_estimator->vertexes[edge_idx * 2 - 1]));
                if (find(last_v_neighbors.begin(), last_v_neighbors.end(), edge[0]) == last_v_neighbors.end()) {
                    this_estimator->fast_fail = true;
                    cout << "     fail!";
                } else {
                    // if pass test, then get ui neighbors
                    msg_vec[msg_idx++] = edge[0];
                }
                cout << endl;
            } catch (const std::out_of_range &oor) {
                this_estimator->fast_fail = true;
                cout << "cannot find!!!! 4\n";
                continue;
            }
        }
        batch_get_data(graph, this_thread, 'v', msg_idx,
                       msg_vec);  // send out the list of edges, and put them in the hash table

        // Step 2.3 test u0 < ui; if pass, get vi neighbors
        cout << "========== test ordering of ui and u0 ==========\n";
        msg_idx = 0;
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // test u0 < ui
            try {
                uint32_t ui_id = this_estimator->vertexes[2 * edge_idx];
                uint32_t ui_degree = this_thread.retrieved.at("v-" + to_string(ui_id)).size();
                cout << "ui: " << ui_id << ", " << ui_degree << " vs u0: " << this_estimator->vertexes[0] << ", "
                     << this_estimator->d_u0 << " ";
                if ((ui_degree < this_estimator->d_u0) ||
                    (ui_degree == this_estimator->d_u0 && ui_id <= this_estimator->vertexes[0])) {
                    this_estimator->fast_fail = true;
                    cout << "ordering fail!";
                } else {
                    // if pass, get vi neighbors
                    msg_vec[msg_idx++] = this_estimator->vertexes[2 * edge_idx + 1];
                }
                cout << endl;
            } catch (const std::out_of_range &oor) {
                this_estimator->fast_fail = true;
                cout << "cannot find!!!! 5\n";
                continue;
            }
        }
        batch_get_data(graph, this_thread, 'v', msg_idx, msg_vec);

        // Step 2.4 test u0 < vi
        cout << "========== test ordering of vi and u0 ==========\n";
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // test u0 < ui
            try {
                uint32_t vi_id = this_estimator->vertexes[2 * edge_idx + 1];
                uint32_t vi_degree = this_thread.retrieved.at("v-" + to_string(vi_id)).size();
                cout << "vi: " << vi_id << ", " << vi_degree << " vs u0: " << this_estimator->vertexes[0] << ", "
                     << this_estimator->d_u0 << " ";
                if ((vi_degree < this_estimator->d_u0) ||
                    (vi_degree == this_estimator->d_u0 && vi_id <= this_estimator->vertexes[0])) {
                    this_estimator->fast_fail = true;
                    cout << "ordering fail!";
                }
                cout << endl;
            } catch (const std::out_of_range &oor) {
                this_estimator->fast_fail = true;
                cout << "cannot find!!!! 6\n";
                continue;
            }
        }
    }

    // Step 3 test w and vk neighboring
    // cout << "========== test w and vk neighboring ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        if (this_estimator->fast_fail)
            continue;

        // get last_v_id and neighbors
        try {
            uint32_t last_v_id = this_estimator->vertexes[N_odd_cycle - 2];
            vector<uint32_t> &last_v_neighbors = this_thread.retrieved.at("v-" + to_string(last_v_id));
            // iterate w
            bool fail_root = true;
            for (uint32_t j = 0; j < this_estimator->children.size(); j++) {
                SamplingTreeNode *this_child = this_estimator->children[j];

                if (this_child->fast_fail)
                    continue;

                uint32_t w_id = this_child->vertexes[0];
                cout << "test neighboring (v_last, w): " << last_v_id << " " << w_id;
                if (find(last_v_neighbors.begin(), last_v_neighbors.end(), w_id) == last_v_neighbors.end()) {
                    this_child->fast_fail = true;
                    cout << "     fail!";
                } else {
                    fail_root = false;
                }
                cout << endl;
            }
            if (fail_root)
                this_estimator->fast_fail = true;
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            cout << "cannot find!!!! 7\n";
            continue;
        }
    }

    // Step 4(not needed in final version): calculate probability
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        for (uint32_t j = 0; j < this_estimator->children.size(); j++) {
            SamplingTreeNode *this_child = this_estimator->children[j];
            num_estimators++;

            if (this_estimator->fast_fail || this_child->fast_fail)
                continue;
            // calculate probablity here!
            total_prob_verse += this_estimator->probability * this_child->probability;
        }
    }
    //TODO Note: we haven't tested the v_set, there could be repeated vertexes

    return make_tuple(total_prob_verse, num_estimators, 0, 0);

    /*
    cout << "After batch get: " << endl;
    for (auto it: this_thread.retrieved) {
        cout << it.first << " : " << it.second.size() << " " ;
        // for (auto it_neighbors: it.second) {
        //     cout << it_neighbors << " ";
        // }
        cout << endl;
    }
     */



    /*
!    // get edge u1, v1 -> store edge
!    // get u1 neighbors -> store u1 neighbors
!    // get v1 neighbors -> store v1 neighbors, decide who is u, who is w
!    // get w neighbors -> store w neighbors, test u1 < w and v1  < w
! Loop through all edges{
    // get edge u2, v2 -> store edge, test v1, u2 neighboring
    // get u2 neighbors -> test u1 < u2
    // get v2 neighbors -> test u1 < v2
}
    // get edge u3, v3 -> store edge, test v2, u3 neighboring
    // get u3 neighbors -> test u1 < u3
    // get v3 neighbors -> test u1 < v3
    // ...
!    // final step, no need to get anything -> test w, vk neighboring
    // if get here, continue = true, else = false    // TODO write into paper, this kind of batching can avoid some edges
!    // if get here, calculate the probability; note: we already tested that this cycle is good + get all neighbors of everybody

    // high level principle, make sure that remainging edge test is separated and use the hash map for everything
*/
    return make_tuple(total_prob_verse, sampling_times, 0, 0);
}

// N-star sampler - batch
vector<SamplingTreeNode *>
star_sampler_batch(uint32_t thread_id, Graph &graph, uint32_t petal_num, uint32_t sampling_times, default_random_engine & rand_generator) {
    vector<SamplingTreeNode *> estimators;
    SamplingThread &this_thread = thread_struct[thread_id];
    // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());
    uint32_t *msg_vec = this_thread.msg_vec;
    uint32_t *msg_vec_backup = this_thread.msg_vec_backup;
    uint32_t msg_idx = 0;
    uint32_t tmp_id;
    for (uint32_t i = 0; i < sampling_times; i++)
        estimators.push_back(new SamplingTreeNode(StarRoot, 0.0)); // probability not finalized


    // Step 1: randomly get edges -> first vertex is the center
    // cout << "========== get an edge ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        // tmp_id = graph.edge_distribution(graph.rand_generator);
        tmp_id = graph.edge_distribution(rand_generator);
        estimators[estimator_idx]->waiting_data = "edge-" + to_string(tmp_id);
        msg_vec[estimator_idx] = tmp_id;
    }
    batch_get_data(graph, this_thread, 'e', sampling_times, msg_vec);

    // Step 2: get neighbors of all centers : parse edge -> batch get center vertex neighbors
    // cout << "========== get center neighbors ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        try {
            vector<uint32_t> &edge = this_thread.retrieved.at(this_estimator->waiting_data);
            this_estimator->vertexes.push_back(edge[0]);
            msg_vec[estimator_idx] = edge[0];
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            cout << "cannot find!!!! 8\n";
            continue;
        }
    }
    batch_get_data(graph, this_thread, 'v', sampling_times, msg_vec);

    // Step 3: get degree, update root probability, pick neighbors, create leave; if petal_num == 1, need to get the chosen neighbor degree, compare
    // cout << "========== calculate probability ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        try {
            // get center neighbors and degree
            uint32_t center_id = this_estimator->vertexes[0];
            vector<uint32_t> &center_neighbors = this_thread.retrieved.at("v-" + to_string(center_id));
            uint32_t center_vertex_degree = center_neighbors.size();
            this_estimator->d_u0 = center_vertex_degree;
            // update root probability
            this_estimator->probability = 1.0 * graph.total_edges_ / center_vertex_degree;
            // pick neighbors, create leave; if petal_num == 1, need to judge to fast fail
            uniform_int_distribution<uint32_t> next_edge_distribution(0, center_vertex_degree - 1);
            set<uint32_t> vertex_set;
            vector<uint32_t> sampled_petals;
            if (center_vertex_degree >= petal_num) {
                // indicate can find enough petals
                for (uint32_t i = 0; i < petal_num; i++) {
                    uint32_t rand_index = next_edge_distribution(rand_generator);
                    while (vertex_set.find(rand_index) != vertex_set.end())
                        rand_index = next_edge_distribution(rand_generator);
                    vertex_set.insert(rand_index);
                    sampled_petals.push_back(center_neighbors[rand_index]);
                }
            } else {
                this_estimator->fast_fail = true;
            }

            uint64_t this_count = 1;
            if (petal_num < center_vertex_degree) {
                for (uint32_t i = center_vertex_degree; i > center_vertex_degree - petal_num; i--)
                    this_count = this_count * i;

                for (uint32_t i = 2; i <= petal_num; i++)
                    this_count = this_count / i;
            }

            this_estimator->children.push_back(new SamplingTreeNode(StarLeave, sampled_petals, 1.0 * this_count));
            //cout << this_estimator->probability << ", " << this_estimator->children[0]->probability << " : ";
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            cout << "cannot find!!!! 9\n";
            continue;
        }
    }
    // cout << endl;

    // finally, if petal_num = 1, need to compare ordering for uniqueness
    /*
    if (petal_num == 1) {
        // get neighbors of the petal
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_child = estimators[estimator_idx]->children[0];
            msg_vec[estimator_idx] = this_child->vertexes[0];
        }
        batch_get_data(graph, this_thread, 'v', sampling_times, msg_vec);

        // compare degree
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            SamplingTreeNode *this_child = estimators[estimator_idx]->children[0];
            try {
                vector<uint32_t> &petal_neighbors = this_thread.retrieved.at("v-" + to_string(this_child->vertexes[0]));
                if (petal_neighbors.size() > this_estimator->d_u0) {
                    this_estimator->fast_fail = true;
                    this_child->fast_fail = true;
                }
            } catch (const std::out_of_range &oor) {
                this_estimator->fast_fail = true;
                cout << "cannot find!!!! 10\n";
                continue;
            }
        }
    }
    */
    
    return estimators;
}

// Odd cycle sampler - batch
vector<SamplingTreeNode *>
odd_cycle_sampler_batch(uint32_t thread_id, Graph &graph, uint32_t N_odd_cycle, uint32_t sampling_times, default_random_engine & rand_generator) {
    vector<SamplingTreeNode *> estimators;
    SamplingThread &this_thread = thread_struct[thread_id];

    // default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());
    for (uint32_t i = 0; i < sampling_times; i++)
        estimators.push_back(new SamplingTreeNode(CycleRoot, pow(graph.total_edges_, (N_odd_cycle / 2)) / 2.0));
    uint32_t *msg_vec = this_thread.msg_vec;
    uint32_t *msg_vec_backup = this_thread.msg_vec_backup;
    uint32_t msg_idx = 0;
    uint32_t tmp_id;

    // Step 1: Handle u1,v1 and related testing ; get u1, v1; decide u1, v1 ; choose w; get w; test u1 < w and v1 < w
    // Step 1.1 get u1, v1 for every estimator
    // cout << "========== get u0, v0 ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        // tmp_id = graph.edge_distribution(graph.rand_generator);
        tmp_id = graph.edge_distribution(rand_generator);
        estimators[estimator_idx]->waiting_data = "edge-" + to_string(tmp_id);
        msg_vec[estimator_idx] = tmp_id;
    }
    batch_get_data(graph, this_thread, 'e', sampling_times, msg_vec);

    // Step 1.2 get u1, v1 neighbors -> create two arrays store u1, v1 neighbors
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        try {
            vector<uint32_t> &edge = this_thread.retrieved.at(this_estimator->waiting_data);
            this_estimator->vertexes.push_back(edge[0]);
            this_estimator->vertexes.push_back(edge[1]);
            msg_vec[estimator_idx] = edge[0];
            msg_vec_backup[estimator_idx] = edge[1];
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            cout << "cannot find!!!! 11\n";
            continue;
        }
    }
    batch_get_data(graph, this_thread, 'v', sampling_times, msg_vec);
    batch_get_data(graph, this_thread, 'v', sampling_times, msg_vec_backup);

    // Step 1.3 compare u0, v0, decide u0, v0, w -> get w neighbors, and create leave SamplingNode
    // cout << "========== compare u0, v0, decide w ==========\n";
    msg_idx = 0;
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        if (this_estimator->fast_fail) continue;

        try {
            uint32_t v0_id = this_estimator->vertexes[0];
            uint32_t v1_id = this_estimator->vertexes[1];
            vector<uint32_t> &neighbors_0 = this_thread.retrieved.at("v-" + to_string(v0_id));
            vector<uint32_t> &neighbors_1 = this_thread.retrieved.at("v-" + to_string(v1_id));
            uint32_t d0 = neighbors_0.size();
            uint32_t d1 = neighbors_1.size();

            // decide u0, v0
            if ((d0 > d1) || (d0 == d1 && v0_id > v1_id)) {
                this_estimator->d_u0 = d1;
                this_estimator->d_v0 = d0;
                this_estimator->vertexes[0] = v1_id;
                this_estimator->vertexes[1] = v0_id;

                // choose w s
                // cout << this_estimator->vertexes[0] << " vs. " << this_estimator->vertexes[1] << " : "
                //      << this_estimator->d_u0 << " " << this_estimator->d_v0;
                uniform_int_distribution<uint32_t> next_edge_distribution(0, d1 - 1);
                for (uint32_t j = 0; j < ceil(this_estimator->d_u0 / graph.sqrt_m); j++) {
                    vector<uint32_t> sampled_vertex_id;
                    sampled_vertex_id.push_back(neighbors_1[next_edge_distribution(rand_generator)]);

                    // cout << " sampled w: " << sampled_vertex_id[0] << endl;

                    this_estimator->children.push_back(new SamplingTreeNode(CycleLeave, sampled_vertex_id, d1));
                    msg_vec[msg_idx++] = sampled_vertex_id[0];
                }
            } else {
                this_estimator->d_u0 = d0;
                this_estimator->d_v0 = d1;

                // choose w s
                // cout << this_estimator->vertexes[0] << " vs. " << this_estimator->vertexes[1] << " : "
                //     << this_estimator->d_u0 << " " << this_estimator->d_v0;
                uniform_int_distribution<uint32_t> next_edge_distribution(0, d0 - 1);
                for (uint32_t j = 0; j < ceil(this_estimator->d_u0 / graph.sqrt_m); j++) {
                    vector<uint32_t> sampled_vertex_id;
                    sampled_vertex_id.push_back(neighbors_0[next_edge_distribution(rand_generator)]);

                    // cout << " sampled w: " << sampled_vertex_id[0] << endl;

                    this_estimator->children.push_back(new SamplingTreeNode(CycleLeave, sampled_vertex_id, d0));
                    msg_vec[msg_idx++] = sampled_vertex_id[0];
                }
            }
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            cout << "cannot find!!!! 12\n";
            continue;
        }

    }
    batch_get_data(graph, this_thread, 'v', msg_idx, msg_vec);

    // Step 1.4 test u1 < w and v1 < w
    // cout << "========== test ordering of w and u0,v0 ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        if (this_estimator->fast_fail) continue;

        // test ordering
        bool fail_root = true;
        for (uint32_t j = 0; j < this_estimator->children.size(); j++) {
            SamplingTreeNode *this_child = this_estimator->children[j];
            try {
                uint32_t w_id = this_child->vertexes[0];
                uint32_t w_degree = this_thread.retrieved.at("v-" + to_string(w_id)).size();
                // test v0 < w only (no need to test u0 < w because u0 < v0 already)
                // cout << "w: " << w_id << ", " << w_degree << " vs v0: " << this_estimator->vertexes[1] << ", "
                //     << this_estimator->d_v0 << " ";
                if ((w_degree < this_estimator->d_v0) ||
                    (w_degree == this_estimator->d_v0 && w_id <= this_estimator->vertexes[1])) {
                    this_child->fast_fail = true;
                    // cout << "ordering fail!";
                } else {
                    fail_root = false;  // one child is not failed, then cannot fail the root
                }
                // cout << endl;
            } catch (const std::out_of_range &oor) {
                this_child->fast_fail = true;
                cout << "cannot find!!!! 13\n";
                continue;
            }
        }

        if (fail_root)
            this_estimator->fast_fail = true;
    }

    // Step 2 get remaining edges, and test ordering and neighboring
    for (uint32_t edge_idx = 1; edge_idx < N_odd_cycle / 2; edge_idx++) {
        // cout << "========== get another edge " << edge_idx << " ==========\n";
        // Step 2.1  get a random edge
        uint32_t tmp_id;
        msg_idx = 0;
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // tmp_id = graph.edge_distribution(graph.rand_generator);
            tmp_id = graph.edge_distribution(rand_generator);
            this_estimator->waiting_data = "edge-" + to_string(tmp_id);
            msg_vec[msg_idx++] = tmp_id;
        }
        batch_get_data(graph, this_thread, 'e', msg_idx, msg_vec);

        // Step 2.2 store ui, vi; test vi-1 and ui neighboring; if pass, get ui neighbors
        // cout << "========== test vi-1, ui neighboring ==========\n";
        msg_idx = 0;
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // parse the edge
            try {
                vector<uint32_t> &edge = this_thread.retrieved.at(this_estimator->waiting_data);
                this_estimator->vertexes.push_back(edge[0]);
                this_estimator->vertexes.push_back(edge[1]);

                // test vi-1 and ui neighboring
                // cout << "neighboring test: vi-1: " << this_estimator->vertexes[edge_idx * 2 - 1] << ", ui: " << edge[0];
                vector<uint32_t> &last_v_neighbors = this_thread.retrieved.at(
                        "v-" + to_string(this_estimator->vertexes[edge_idx * 2 - 1]));
                if (find(last_v_neighbors.begin(), last_v_neighbors.end(), edge[0]) == last_v_neighbors.end()) {
                    this_estimator->fast_fail = true;
                    // cout << "     fail!";
                } else {
                    // if pass test, then get ui neighbors
                    msg_vec[msg_idx++] = edge[0];
                }
                // cout << endl;
            } catch (const std::out_of_range &oor) {
                this_estimator->fast_fail = true;
                cout << "cannot find!!!! 14\n";
                continue;
            }
        }
        batch_get_data(graph, this_thread, 'v', msg_idx, msg_vec);

        // Step 2.3 test u0 < ui; if pass, get vi neighbors
        // cout << "========== test ordering of ui and u0 ==========\n";
        msg_idx = 0;
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // test u0 < ui
            try {
                uint32_t ui_id = this_estimator->vertexes[2 * edge_idx];
                uint32_t ui_degree = this_thread.retrieved.at("v-" + to_string(ui_id)).size();
                // cout << "ui: " << ui_id << ", " << ui_degree << " vs u0: " << this_estimator->vertexes[0] << ", "
                //     << this_estimator->d_u0 << " ";
                if ((ui_degree < this_estimator->d_u0) ||
                    (ui_degree == this_estimator->d_u0 && ui_id <= this_estimator->vertexes[0])) {
                    this_estimator->fast_fail = true;
                    // cout << "ordering fail!";
                } else {
                    // if pass, get vi neighbors
                    msg_vec[msg_idx++] = this_estimator->vertexes[2 * edge_idx + 1];
                }
                // cout << endl;
            } catch (const std::out_of_range &oor) {
                this_estimator->fast_fail = true;
                cout << "cannot find!!!! 15\n";
                continue;
            }
        }
        batch_get_data(graph, this_thread, 'v', msg_idx, msg_vec);

        // Step 2.4 test u0 < vi
        // cout << "========== test ordering of vi and u0 ==========\n";
        for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
            SamplingTreeNode *this_estimator = estimators[estimator_idx];
            if (this_estimator->fast_fail)
                continue;

            // test u0 < ui
            try {
                uint32_t vi_id = this_estimator->vertexes[2 * edge_idx + 1];
                uint32_t vi_degree = this_thread.retrieved.at("v-" + to_string(vi_id)).size();
                // cout << "vi: " << vi_id << ", " << vi_degree << " vs u0: " << this_estimator->vertexes[0] << ", "
                //     << this_estimator->d_u0 << " ";
                if ((vi_degree < this_estimator->d_u0) ||
                    (vi_degree == this_estimator->d_u0 && vi_id <= this_estimator->vertexes[0])) {
                    this_estimator->fast_fail = true;
                //    cout << "ordering fail!";
                }
                // cout << endl;
            } catch (const std::out_of_range &oor) {
                this_estimator->fast_fail = true;
                cout << "cannot find!!!! 16\n";
                continue;
            }
        }
    }

    // Step 3 test w and vk neighboring
    // cout << "========== test w and vk neighboring ==========\n";
    for (uint32_t estimator_idx = 0; estimator_idx < sampling_times; estimator_idx++) {
        SamplingTreeNode *this_estimator = estimators[estimator_idx];
        if (this_estimator->fast_fail)
            continue;

        // get last_v_id and neighbors
        try {
            uint32_t last_v_id = this_estimator->vertexes[N_odd_cycle - 2];
            vector<uint32_t> &last_v_neighbors = this_thread.retrieved.at("v-" + to_string(last_v_id));
            // iterate w
            bool fail_root = true;
            for (uint32_t j = 0; j < this_estimator->children.size(); j++) {
                SamplingTreeNode *this_child = this_estimator->children[j];

                if (this_child->fast_fail)
                    continue;

                uint32_t w_id = this_child->vertexes[0];
                // cout << "test neighboring (v_last, w): " << last_v_id << " " << w_id;
                if (find(last_v_neighbors.begin(), last_v_neighbors.end(), w_id) == last_v_neighbors.end()) {
                    this_child->fast_fail = true;
                    // cout << "     fail!";
                } else {
                    fail_root = false;
                }
                // cout << endl;
            }
            if (fail_root)
                this_estimator->fast_fail = true;
        } catch (const std::out_of_range &oor) {
            this_estimator->fast_fail = true;
            // cout << "cannot find!!!! 17\n";
            continue;
        }
    }

    return estimators;
}

/* Create sampling trees
 * in batch of num_estimators
*/
vector<SamplingTreeNode *>
construct_sampling_tree_batch(uint32_t thread_id, Graph &graph, Pattern &pattern, uint32_t num_estimators, default_random_engine & rand_generator) {
    SamplingThread &this_thread = thread_struct[thread_id];
    vector<SamplingTreeNode *> roots;
    vector<SamplingTreeNode *> leaves;
    vector<SamplingTreeNode *> new_leaves;
    bool tree_empty = true;

    for (uint32_t i = 0; i < pattern.blocks.size(); i++)
    {
        if (pattern.blocks[i].type == Cycle)
        {
            if (tree_empty) {
                // batch sample some odd-cycle samplers -> update roots, and leaves
                roots = odd_cycle_sampler_batch(thread_id, graph, pattern.cycles[pattern.blocks[i].idx].size(), num_estimators, rand_generator);

                for (SamplingTreeNode *it: roots) {
                    for (uint32_t j = 0; j < it->children.size(); j++) {
                        if (early_pruning == true) {
                            if(it->fast_fail == false && it->children[j]->fast_fail == false)
                                leaves.push_back(it->children[j]);
                        } else {
                            leaves.push_back(it->children[j]);
                        }
                    }
                }
                tree_empty = false;
            } else {
                // batch sample leaves.size() odd-cycle samplers -> extend it to every leaf and update new leaves
                vector<SamplingTreeNode *> new_odd_cycle_roots = odd_cycle_sampler_batch(thread_id, graph,
                                                                                        pattern.cycles[pattern.blocks[i].idx].size(),
                                                                                        leaves.size(), rand_generator);
                for (uint32_t j = 0; j < leaves.size(); j++)
                    leaves[j]->children.push_back(new_odd_cycle_roots[j]);

                // update leaves and new_leaves
                new_leaves.clear();
                for (uint32_t j = 0; j < leaves.size(); j++) {
                    SamplingTreeNode *last_level_root = leaves[j]->children[0];
                    for (uint32_t k = 0; k < last_level_root->children.size(); k++) {
                        if (early_pruning == true) {
                            if (last_level_root->fast_fail == false && last_level_root->children[k]->fast_fail == false)
                                new_leaves.push_back(last_level_root->children[k]);
                        } else {
                            new_leaves.push_back(last_level_root->children[k]);
                        }
                    }
                }
                leaves = new_leaves;
            }
        }
        else 
        {
            if (tree_empty) {
                roots = star_sampler_batch(thread_id, graph, pattern.stars[pattern.blocks[i].idx].size() - 1, num_estimators, rand_generator);

                for (SamplingTreeNode *it: roots) {
                    if (early_pruning == true) {
                        if(it->fast_fail == false && it->children[0]->fast_fail == false)
                            leaves.push_back(it->children[0]);
                    } else {
                        leaves.push_back(it->children[0]);
                    }
                }
                tree_empty = false;
            } else {
                vector<SamplingTreeNode *> new_star_roots = star_sampler_batch(thread_id, graph,
                                                                            pattern.stars[pattern.blocks[i].idx].size() - 1, leaves.size(), rand_generator);

                for (uint32_t j = 0; j < leaves.size(); j++)
                    leaves[j]->children.push_back(new_star_roots[j]);

                // update leaves and new_leaves
                new_leaves.clear();
                for (uint32_t j = 0; j < leaves.size(); j++) {
                    if (early_pruning == true) {
                        if (leaves[j]->children[0]->fast_fail == false && leaves[j]->children[0]->children[0]->fast_fail == false)
                            new_leaves.push_back(leaves[j]->children[0]->children[0]);
                    } else {
                        new_leaves.push_back(leaves[j]->children[0]->children[0]);
                    }
                }
                leaves = new_leaves;
            }
        }
    }
    return roots;
}


void DFS_judge_remaining_edges_batch(uint32_t thread_id, uint32_t idx, Graph &graph, Pattern &pattern,
                                     vector<vector<uint32_t>> &cycle_vertexes, vector<vector<uint32_t>> &star_vertexes,
                                     map<uint32_t, uint32_t> &v_map, vector<succeed_conditions> &conditions,
                                     double prob) {
    // corner cases
    if (idx == cycle_vertexes.size() + star_vertexes.size()) {
        // judge remaining edges
        succeed_conditions this_condition;
        this_condition.probability = prob;
        for (uint32_t i = 0; i < pattern.other_edges_.size(); i++) {
            uint32_t test_v_start = v_map.at(pattern.other_edges_[i].v_start);
            uint32_t test_v_end = v_map.at(pattern.other_edges_[i].v_end);
            this_condition.neighbors_need.push_back(make_pair(test_v_start, test_v_end));
        }
        conditions.push_back(this_condition);
        return;
    }

    // normal grow mapping
    if (idx < cycle_vertexes.size()) {
        for (uint32_t i = 0; i < pattern.candidate_cycle_v_maps[idx].size(); i++) {
            // map candidate_cycle_v_maps[idx][i] to cycle_vertexes[idx]
            for (uint32_t j = 0; j < pattern.candidate_cycle_v_maps[idx][i].size(); j++) {
                v_map.emplace(pattern.candidate_cycle_v_maps[idx][i][j], cycle_vertexes[idx][j]);
            }

            // DFS_judge_further
            DFS_judge_remaining_edges_batch(thread_id, idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map,
                                            conditions, prob);

            // restore
            for (uint32_t j = 0; j < pattern.candidate_cycle_v_maps[idx][i].size(); j++) {
                v_map.erase(pattern.candidate_cycle_v_maps[idx][i][j]);
            }
        }
    } else {
        uint32_t star_idx = idx - cycle_vertexes.size();
        for (uint32_t i = 0; i < pattern.candidate_star_v_maps[star_idx].size(); i++) {
            // map candidate_star_v_maps[idx][i] to star_vertexes[idx]
            for (uint32_t j = 0; j < pattern.candidate_star_v_maps[star_idx][i].size(); j++) {
                v_map.emplace(pattern.candidate_star_v_maps[star_idx][i][j], star_vertexes[star_idx][j]);
            }

            // DFS_judge_further
            DFS_judge_remaining_edges_batch(thread_id, idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map,
                                            conditions, prob);

            // restore
            for (uint32_t j = 0; j < pattern.candidate_star_v_maps[star_idx][i].size(); j++) {
                v_map.erase(pattern.candidate_star_v_maps[star_idx][i][j]);
            }
        }
    }
}


void judge_remaining_edges_batch(uint32_t thread_id, Graph &graph, Pattern &pattern,
                                 vector<vector<uint32_t>> &cycle_vertexes, vector<vector<uint32_t>> &star_vertexes,
                                 vector<succeed_conditions> &conditions, double prob) {
    if (pattern.other_edges_.size() == 0) {
        succeed_conditions this_condition;
        this_condition.probability = prob;
        conditions.push_back(this_condition);
        return;
    }

    map<uint32_t, uint32_t> v_map;
    DFS_judge_remaining_edges_batch(thread_id, 0, graph, pattern, cycle_vertexes, star_vertexes, v_map, conditions,
                                    prob);

    return;
}

/* Major helper function to iterate a sampling tree
*/
void tree_DFS_helper_batch(uint32_t thread_id, SamplingTreeNode *root, Graph &graph, Pattern &pattern,
                           vector<SamplingTreeNode *> &path, map<uint32_t, uint32_t> &v_map,
                           uint64_t &num_cross_partition, vector<succeed_conditions> &conditions) {
    // corner cases
    if (root->children.size() == 0) {
        path.push_back(root);
        v_map.clear();

        // task 1: print path
        if (path.size() < (pattern.cycle_num + pattern.star_num) * 2) {
            num_cross_partition += 1;
            return;
        }
        // assert(path.size() == (pattern.cycle_num + pattern.star_num) * 2);

        if (FLAG_print) {
            cout << "found a path: ";
        }

        set<uint32_t> vertex_set;
        bool test_succeed = true;
        double prob_inverse = 1.0;
        num_cross_partition += 1;


        // task 2, 3: judge cycles and stars
        vector<vector<uint32_t>> cycle_vertexes;
        vector<vector<uint32_t>> star_vertexes;
        for (uint32_t i = 0; i < pattern.block_num; i++)
        {
            if (test_succeed == false) 
                break;
            if (pattern.blocks[i].type == Cycle)
            {
                pattern.blocks[i].total_sampled_times += 1;
                set<uint32_t> cycle_vertex_set;
                SamplingTreeNode *cur_root = path[i * 2];
                SamplingTreeNode *cur_leave = path[i * 2 + 1];

                if (FLAG_print) {
                    cur_root->print();
                    cout << ", ";
                    cur_leave->print();
                    cout << ", ";
                }

                // handle cycle root
                prob_inverse *= cur_root->probability;

                // map this node's vertexes and cycles[i]'s first several nodes
                vector<uint32_t> cur_cycle;
                for (uint32_t j = 0; j < cur_root->vertexes.size(); j++) {
                    v_map.emplace(pattern.cycles[pattern.blocks[i].idx][j], cur_root->vertexes[j]);
                    vertex_set.insert(cur_root->vertexes[j]);
                    cur_cycle.push_back(cur_root->vertexes[j]);
                    cycle_vertex_set.insert(cur_root->vertexes[j]);
                }

                // handle cycle leave
                prob_inverse *= cur_leave->probability;

                v_map.emplace(pattern.cycles[pattern.blocks[i].idx].back(), cur_leave->vertexes[0]);
                vertex_set.insert(cur_leave->vertexes[0]);
                cur_cycle.push_back(cur_leave->vertexes[0]);
                cycle_vertexes.push_back(cur_cycle);
                cycle_vertex_set.insert(cur_leave->vertexes[0]);

                // judge cycle construct
                // test overlapped vertex
                if (cycle_vertex_set.size() != pattern.cycles[pattern.blocks[i].idx].size() || cur_root->fast_fail == true || cur_leave->fast_fail == true) {
                    test_succeed = false;
                    pattern.blocks[i].fail_probability += 1;
                    if (FLAG_print)
                        cout << "cycle " << i << " has overlapped vertexes" << endl;
                }
            }
            else 
            {
                SamplingTreeNode *cur_root = path[i * 2];
                SamplingTreeNode *cur_leave = path[i * 2 + 1];
                pattern.blocks[i].total_sampled_times += 1;

                if (FLAG_print) {
                    cur_root->print();
                    cout << ", ";
                    cur_leave->print();
                    cout << ", ";
                }

                // handle star root
                prob_inverse *= cur_root->probability;
                v_map.emplace(pattern.stars[pattern.blocks[i].idx].front(), cur_root->vertexes[0]);
                vector<uint32_t> cur_star;
                cur_star.push_back(cur_root->vertexes[0]);
                vertex_set.insert(cur_root->vertexes[0]);

                // handle leaves and and test copy
                prob_inverse *= cur_leave->probability;
                uint32_t petal_num = cur_leave->vertexes.size();
                if (petal_num == 0 || cur_root->fast_fail == true || cur_leave->fast_fail == true) {
                    test_succeed = false;
                    pattern.blocks[i].fail_probability += 1;
                    if (FLAG_print)
                        cout << "star " << i << " does not have enough petals" << endl;
                } else {
                    for (uint32_t j = 0; j < cur_leave->vertexes.size(); j++) {
                        v_map.emplace(pattern.stars[pattern.blocks[i].idx][j + 1], cur_leave->vertexes[j]);
                        cur_star.push_back(cur_leave->vertexes[j]);
                        vertex_set.insert(cur_leave->vertexes[j]);
                    }
                    star_vertexes.push_back(cur_star);
                }
            }
        }

        // print v_map
        // assert(v_map.size() == pattern.vertex_count);

        // task 4: judge other edges
        if (test_succeed == false) {
            if (FLAG_print)
                cout << "previous cycle/star tests already failed" << endl;
        } else if (vertex_set.size() != pattern.vertex_count) {
            if (FLAG_print)
                cout << "this path has overlapped sampled vertexes" << endl;
            test_succeed = false;
        } else {
            // judge other edges
            judge_remaining_edges_batch(thread_id, graph, pattern, cycle_vertexes, star_vertexes, conditions, prob_inverse);
        }
        path.pop_back();
        return;
    }

    // DFS to iterate children
    path.push_back(root);
    for (uint32_t i = 0; i < root->children.size(); i++)
        tree_DFS_helper_batch(thread_id, root->children[i], graph, pattern, path, v_map, num_cross_partition, conditions);
    path.pop_back();

    return;
}

/* Count a sampling tree
*/
void sampling_tree_count_batch(uint32_t thread_id, vector<SamplingTreeNode *> &roots, Graph &graph, Pattern &pattern,
                               double &this_cross_prob_inverse, uint64_t &num_cross_partition,
                               double &num_succeed, uint64_t &num_intra_partition) {
    // cout << "====== start to count batch of trees ======" << endl;

    SamplingThread &this_thread = thread_struct[thread_id];
    uint32_t *msg_vec = this_thread.msg_vec;
    uint32_t msg_idx = 0;


    // Step 2.1: iterate every tree and get a bunch of probability and conditions (neighboring of some edges);
    vector<succeed_conditions> conditions;
    for (SamplingTreeNode *root: roots) {
        vector<SamplingTreeNode *> path;
        map<uint32_t, uint32_t> v_map;
        // create conditions for succeed, and update num_cross_partition for num_estimators
        tree_DFS_helper_batch(thread_id, root, graph, pattern, path, v_map, num_cross_partition, conditions);
    }

    // Step 2.2: iterate the conditions, if need some information, create the message
    msg_idx = 0;
    for (succeed_conditions &condition: conditions) {
        for (pair<uint32_t, uint32_t> &need_neighbors: condition.neighbors_need) {
            uint32_t v1 = need_neighbors.first;
            uint32_t v2 = need_neighbors.second;
            if (this_thread.retrieved.find("v-" + to_string(v1)) == this_thread.retrieved.end() &&
                this_thread.retrieved.find("v-" + to_string(v2)) == this_thread.retrieved.end()) {
                msg_vec[msg_idx++] = v1;
            }
        }
    }
    batch_get_data(graph, this_thread, 'v', msg_idx, msg_vec);

    // Step 2.3: iterate again
    for (succeed_conditions &condition: conditions) {
        bool succeed = true;
        // cout << "get a condition " << condition.probability << ";\n";
        for (pair<uint32_t, uint32_t> &need_neighbors: condition.neighbors_need) {
            uint32_t v1 = need_neighbors.first;
            uint32_t v2 = need_neighbors.second;
            if (this_thread.retrieved.find("v-" + to_string(v1)) != this_thread.retrieved.end()) {
                vector<uint32_t> &v1_neighbors = this_thread.retrieved.at("v-" + to_string(v1));
                if (find(v1_neighbors.begin(), v1_neighbors.end(), v2) == v1_neighbors.end())
                    succeed = false;
            } else if (this_thread.retrieved.find("v-" + to_string(v2)) != this_thread.retrieved.end()) {
                vector<uint32_t> &v2_neighbors = this_thread.retrieved.at("v-" + to_string(v2));
                if (find(v2_neighbors.begin(), v2_neighbors.end(), v1) == v2_neighbors.end())
                    succeed = false;
            } else {
                cout << "Cannot find anybody!\n";
                succeed = false;
            }
        }

        if (succeed) {
            this_cross_prob_inverse += condition.probability;
            num_succeed += 1;
        }
    }
    return;
}

/* Estimating arbitrary patterns
* from graph, with sampling_times, and all samplers are executed vectorized, to batch network requests
*/
tuple<double, uint64_t, double, uint64_t>
estimate_pattern_batch(uint32_t thread_id, Graph &graph, Pattern &pattern, uint32_t sampling_times, default_random_engine & rand_generator) {
    double this_cross_prob_inverse = 0.0;
    double num_succeed = 0.0;
    uint64_t num_cross_partition = 0;
    uint64_t num_intra_partition = 0;

    // Step 1: construct the trees
    vector<SamplingTreeNode *> roots = construct_sampling_tree_batch(thread_id, graph, pattern, sampling_times, rand_generator);

    // Step 2: iterate every path from root to leave
    sampling_tree_count_batch(thread_id, roots, graph, pattern, this_cross_prob_inverse, num_cross_partition,
                              num_succeed, num_intra_partition);

    return make_tuple(this_cross_prob_inverse, num_cross_partition, num_succeed, num_intra_partition);
}

void adjust_sampling_block_order(Pattern & pattern)
{
    sort(pattern.blocks.begin(), pattern.blocks.end());
}

#endif
