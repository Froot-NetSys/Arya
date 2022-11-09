// estimating.h
#ifndef ESTIMATING_H
#define ESTIMATING_H

#include <iostream>
#include "graph.h"
#include <cmath>
#include <set>
#include <queue>
#include <map>
#include "assert.h"

using namespace std;

/* SamplingTree strcuture
 * tree of nodes (vertex, probability)
*/
struct SamplingTreeNode {
    SamplingTreeNode(NodeType type_, vector <uint32_t> &vertexes_, double probability_) : node_type(type_),
                                                                                          probability(probability_) {
        for (uint32_t i = 0; i < vertexes_.size(); i++) {
            vertexes.push_back(vertexes_[i]);
        }
    }

    void print() {
        cout << "Node(" << NodeTypeStr[node_type] << ", " << children.size() << ", " << probability << ", [ ";
        for (auto it = vertexes.begin(); it != vertexes.end(); it++)
            cout << *it << " ";
        cout << "]) ";
    }

    NodeType node_type;   // CYCLE_ROOT, CYCLE_LEAF, STAR_ROOT, STAR_LEAF
    vector <uint32_t> vertexes;
    vector<struct SamplingTreeNode *> children;
    double probability;
};


/* Sample a odd-cycle
* from graph, <u1, u2, ..., w1>
* node([(v1, v2), (v3, v4), ...]) -> k x w1
*/
SamplingTreeNode *odd_cycle_sampler(Graph &graph, uint32_t cycle_length) {
    //cout << "====== Start to sample a " << cycle_length << "-cycles" << endl;
    default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());
    uint32_t ue_id, d_ue;

    // step 1: sample cycle_length/2 edges (u1, v1), (x, x), (x, x), ... u1 < v1
    vector <uint32_t> sample_edge_vector;
    Edge e1 = graph.edge_sampling();
    if (graph.vertex_before(e1.v_start, e1.v_end)) {
        sample_edge_vector.push_back(e1.v_start);
        sample_edge_vector.push_back(e1.v_end);
    } else {
        sample_edge_vector.push_back(e1.v_end);
        sample_edge_vector.push_back(e1.v_start);
    }

    for (uint32_t j = 0 + 1; j < cycle_length / 2; j++) {
        Edge sampled_edge = graph.edge_sampling();
        sample_edge_vector.push_back(sampled_edge.v_start);
        sample_edge_vector.push_back(sampled_edge.v_end);
    }

    SamplingTreeNode *root = new SamplingTreeNode(CycleRoot,
                                                  sample_edge_vector,
                                                  pow(graph.total_edges_, (cycle_length / 2)) / 2.0);
    assert(root != nullptr);


    // step 2: sample w1, a neighbor of u1
    ue_id = sample_edge_vector[0];
    d_ue = graph.degree(ue_id);

    // cout << "cycle u1 degree: " << d_ue << endl;

    uniform_int_distribution <uint32_t> next_edge_distribution(0, d_ue - 1);

    for (uint32_t j = 0; j < ceil(d_ue / graph.sqrt_m); j++) {
        // step 2: sample an vertex from u1's neighbor
        vector <uint32_t> sampled_vertex_id;
        sampled_vertex_id.push_back(graph.retrieve_neighbor(ue_id, next_edge_distribution(rand_generator)));
        root->children.push_back(new SamplingTreeNode(CycleLeave, sampled_vertex_id, d_ue));
    }

    return root;
}

/* Sample a star
* from graph, <center, w1, w2, ...>
*/
SamplingTreeNode *star_sampler(Graph &graph, uint32_t petal_num) {
    //cout << "====== Start to sample an " << petal_num << "-star" << endl;
    default_random_engine rand_generator(std::chrono::system_clock::now().time_since_epoch().count());

    // Step 1: choose a vertex, proportional to its degree, choose an edge, its first vertex
    Edge e1 = graph.edge_sampling();
    uint32_t center_vertex_id = e1.v_start;
    uint32_t center_vertex_degree = graph.degree(center_vertex_id);
    vector <uint32_t> center_vertex;
    center_vertex.push_back(center_vertex_id);
    SamplingTreeNode *root = new SamplingTreeNode(StarRoot,
                                                  center_vertex,
                                                  1.0 * graph.total_edges_ / center_vertex_degree);

    // Step 2: sample petal_num vertexes from N(X), without replacement
    set <uint32_t> vertex_set;
    uniform_int_distribution <uint32_t> next_edge_distribution(0, center_vertex_degree - 1);

    vector <uint32_t> sampled_petals;
    if (center_vertex_degree >= petal_num) {
        // indicate can find enough petals
        for (uint32_t i = 0; i < petal_num; i++) {
            uint32_t rand_index = next_edge_distribution(rand_generator);
            while (vertex_set.find(rand_index) != vertex_set.end()) {
                rand_index = next_edge_distribution(rand_generator);
            }
            vertex_set.insert(rand_index);
            uint32_t sampled_vertex_id = graph.retrieve_neighbor(center_vertex_id, rand_index);
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
    queue < SamplingTreeNode * > to_visit_list;
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
SamplingTreeNode *construct_sampling_tree(Graph &graph, Pattern &pattern) {
    SamplingTreeNode *root = nullptr;
    vector < SamplingTreeNode * > leaves;
    vector < SamplingTreeNode * > new_leaves;

    // Step 1: iterate odd-cycles to add 2-levels
    // cout << "grow tree with cycles" << endl;
    for (uint32_t i = 0; i < pattern.cycle_num; i++) {
        if (i == 0) {
            root = odd_cycle_sampler(graph, pattern.cycles[i].size());
            for (uint32_t j = 0; j < root->children.size(); j++) {
                leaves.push_back(root->children[j]);
            }
        } else {
            // extend every leaf by a cycle tree
            for (uint32_t j = 0; j < leaves.size(); j++) {
                leaves[j]->children.push_back(odd_cycle_sampler(graph, pattern.cycles[i].size()));
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
            root = star_sampler(graph, pattern.stars[i].size() - 1);
            for (uint32_t j = 0; j < root->children.size(); j++) {
                leaves.push_back(root->children[j]);
            }
        } else {
            // extend every leaf by a cycle tree
            for (uint32_t j = 0; j < leaves.size(); j++) {
                leaves[j]->children.push_back(star_sampler(graph, pattern.stars[i].size() - 1));
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
bool test_circle_copy(Graph &graph, vector <uint32_t> &e_vertexes, uint32_t w) {
    // test u1 is the smallest vertex
    uint32_t u1 = e_vertexes[0];
    uint32_t v1 = e_vertexes[1];
    for (uint32_t i = 1; i < e_vertexes.size(); i++) {
        if (!graph.vertex_before(u1, e_vertexes[i]))
            return false;
    }

    if (!graph.vertex_before(u1, w))
        return false;

    // test v1 < w
    if (!graph.vertex_before(v1, w))
        return false;

    // test edges of (v1,u2), ... (vk, w)
    if (!graph.neighbor_test(e_vertexes.back(), w))
        return false;

    for (uint32_t i = 1; i < e_vertexes.size() - 1; i += 2) {
        if (!graph.neighbor_test(e_vertexes[i], e_vertexes[i + 1]))
            return false;
    }

    return true;
}

/* Major helper function to test a sampled path
*/
void DFS_judge_remaining_edges(uint32_t idx, Graph &graph, Pattern &pattern,
                               vector <vector<uint32_t>> &cycle_vertexes, vector <vector<uint32_t>> &star_vertexes,
                               map <uint32_t, uint32_t> &v_map, uint32_t &found_instances, uint32_t *lost_edge_cases) {
    // corner cases
    if (idx == cycle_vertexes.size() + star_vertexes.size()) {

        // judge remaining edges
        bool test_succeed = true;
        // uint32_t lost_edges = 0;
        for (uint32_t i = 0; i < pattern.other_edges_.size(); i++) {
            uint32_t test_v_start = v_map.at(pattern.other_edges_[i].v_start);
            uint32_t test_v_end = v_map.at(pattern.other_edges_[i].v_end);

            if (!graph.neighbor_test(test_v_start, test_v_end)) {
                test_succeed = false;
                // lost_edges += 1;
                break; // early return 
            }
        }

        if (test_succeed) {
            found_instances += 1;
        }
        // lost_edge_cases[lost_edges] += 1;
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
            DFS_judge_remaining_edges(idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map, found_instances, lost_edge_cases);

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
            DFS_judge_remaining_edges(idx + 1, graph, pattern, cycle_vertexes, star_vertexes, v_map, found_instances, lost_edge_cases);

            // restore
            for (uint32_t j = 0; j < pattern.candidate_star_v_maps[star_idx][i].size(); j++) {
                v_map.erase(pattern.candidate_star_v_maps[star_idx][i][j]);
            }
        }
    }
}


uint32_t judge_remaining_edges(Graph &graph, Pattern &pattern,
                               vector <vector<uint32_t>> &cycle_vertexes, vector <vector<uint32_t>> &star_vertexes, uint32_t *lost_edge_cases) {
    /*
    for (uint32_t i = 0; i < cycle_vertexes.size(); i++) {
        for (uint32_t j = 0; j < cycle_vertexes[i].size(); j++)
            cout << cycle_vertexes[i][j] << " ";
        cout << endl;
    }

    for (uint32_t i = 0; i < star_vertexes.size(); i++) {
        for (uint32_t j = 0; j < star_vertexes[i].size(); j++)
            cout << star_vertexes[i][j] << " ";
        cout << endl;
    } */

    if (pattern.other_edges_.size() == 0) {
        return true;
    }

    uint32_t found_instances = 0;
    map <uint32_t, uint32_t> v_map;
    DFS_judge_remaining_edges(0, graph, pattern, cycle_vertexes, star_vertexes, v_map, found_instances, lost_edge_cases);

    return found_instances;
}

/* Major helper function to iterate a sampling tree
*/
void tree_DFS_helper(SamplingTreeNode *root, Graph &graph, Pattern &pattern,
                     vector<SamplingTreeNode *> &path, map <uint32_t, uint32_t> &v_map, double &total_prob_inverse, uint64_t & total_sampling_times,
                     uint32_t *lost_edge_cases) {
    // corner cases
    if (root->children.size() == 0) {
        path.push_back(root);
        v_map.clear();

        // task 1: print path
        assert(path.size() == (pattern.cycle_num + pattern.star_num) * 2);

        if (FLAG_print) {
            cout << "found a path: ";
        }

        set <uint32_t> vertex_set;
        bool test_succeed = true;
        double prob_inverse = 1.0;
        total_sampling_times += 1;

        // task 2: judge cycles
        vector <vector<uint32_t>> cycle_vertexes;
        for (uint32_t i = 0; i < pattern.cycle_num; i++) {
            set <uint32_t> cycle_vertex_set;
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
            vector <uint32_t> cur_cycle;
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
            if (!test_circle_copy(graph, cur_root->vertexes, cur_leave->vertexes[0])) {
                test_succeed = false;
                if (FLAG_print)
                    cout << "cycle " << i << " does not pass copy test" << endl;
            }
        }

        // task 3: judge stars
        vector <vector<uint32_t>> star_vertexes;
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
                vector <uint32_t> cur_star;
                cur_star.push_back(cur_root->vertexes[0]);
                vertex_set.insert(cur_root->vertexes[0]);

                // handle leaves and and test copy
                prob_inverse *= cur_leave->probability;
                uint32_t petal_num = cur_leave->vertexes.size();
                if (petal_num == 0) {
                    test_succeed = false;
                    if (FLAG_print)
                        cout << "star " << i << " does not have enough petals" << endl;
                } else if (petal_num == 1 && !graph.vertex_before(cur_root->vertexes[0], cur_leave->vertexes[0])) {
                    test_succeed = false;
                    if (FLAG_print)
                        cout << "star " << i << " does not form a copy" << endl;
                } else {
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

            uint32_t num_matches = judge_remaining_edges(graph, pattern, cycle_vertexes, star_vertexes, lost_edge_cases);
            if (num_matches == 0) {
                test_succeed = false;
            } else {
                prob_inverse *= num_matches;
                // cout << num_matches << " " << prob_inverse << endl;
            }

        }

        // task 4: if pass test, count
        if (test_succeed) {
            if (FLAG_print) {
                cout << "the path passed the test!" << endl;
            }
            total_prob_inverse += prob_inverse;
        }

        path.pop_back();
        return;
    }

    // DFS to iterate children
    path.push_back(root);
    for (uint32_t i = 0; i < root->children.size(); i++) {
        tree_DFS_helper(root->children[i], graph, pattern, path, v_map, total_prob_inverse, total_sampling_times, lost_edge_cases);
    }
    path.pop_back();

    return;
}

/* Count a sampling tree
*/
void sampling_tree_count(SamplingTreeNode *root, Graph &graph, Pattern &pattern, double &total_prob_inverse, uint64_t & total_sampling_times, uint32_t *lost_edge_cases) {
    if (FLAG_print)
        cout << "====== start to count this tree" << endl;

    // DFS iterate from root to a leave
    vector < SamplingTreeNode * > path;
    map <uint32_t, uint32_t> v_map;

    tree_DFS_helper(root, graph, pattern, path, v_map, total_prob_inverse, total_sampling_times, lost_edge_cases);

    return;
}

void free_sampling_tree(SamplingTreeNode *root)
{
    if (root->children.size() == 0)
    {
        delete root;
        root = nullptr;
        return ;
    }

    for (uint32_t i = 0; i < root->children.size(); i++)
    {
        free_sampling_tree(root->children[i]);
    }
    
    delete root;
    root = nullptr;
}

/* Estimating arbitrary patterns
* from graph, with sampling_times
*/
tuple<double, uint64_t> estimate_pattern(Graph &graph, Pattern &pattern, uint64_t sampling_times, uint32_t *lost_edge_cases) {
    //cout << "= Start to estimate patterns, " << sampling_times << " times" << endl;
    double total_prob_inverse = 0.0;
    double this_prob_inverse = 0.0;
    uint64_t this_num_leaves = 0;
    double this_estimated_pattern_count = 0.0;

    for (uint64_t i = 0; i < sampling_times; i++) {
        // Step 1: construct the trees
        SamplingTreeNode *root = construct_sampling_tree(graph, pattern);

        // Step 2: iterate every path from root to leave
        sampling_tree_count(root, graph, pattern, this_prob_inverse, this_num_leaves, lost_edge_cases);
        
        // Step 3: free sampling tree
        free_sampling_tree(root);

        this_estimated_pattern_count = this_prob_inverse / this_num_leaves;
        total_prob_inverse += this_estimated_pattern_count;

    }

    // cout << total_prob_inverse << " " << total_sampling_times << endl;

    return make_tuple(total_prob_inverse, sampling_times);
}

/* ASAP 
   chain estimator
*/
bool ASAP_ConditionalClose(Graph & graph, Subgraph & pattern, Subgraph & subgraph)
{
    if (pattern.total_vertices != subgraph.total_vertices)
        return false;

    return true;
}

double ASAP_chain_neighbor_sampler(Graph &graph, uint32_t chain_length, default_random_engine & rand_generator)
{   
    
    
    double prob_inverse = 0;
    Subgraph subgraph;
    map<uint32_t, Vertex> open_vertex_set;
    
    uint32_t last_edge_idx = graph.ASAP_edge_sampling(rand_generator);
    Edge e1 = graph.edge_list_[last_edge_idx];
    subgraph.insert_edge(e1);
    open_vertex_set.emplace(e1.v_start, graph.v_table_.at(e1.v_start));
    open_vertex_set.emplace(e1.v_end, graph.v_table_.at(e1.v_end));

    prob_inverse = graph.total_edges_;

    uint32_t total_d_extending_edge = 0;

    for (uint32_t j = 0; j < chain_length - 2; j++)
    {
        total_d_extending_edge = 0;
        vector<uint32_t> open_set_neighbor;
        for (auto it = open_vertex_set.begin(); it != open_vertex_set.end(); it++)
        {
            // uint32_t begin_idx = max(last_edge_idx + 1, it->second.edge_index);
            uint32_t begin_idx = it->second.edge_index;
            uint32_t end_idx = it->second.edge_index + it->second.degree;
            for (uint32_t idx = begin_idx; idx < end_idx; idx++) 
            {
                struct Edge edge = graph.edge_list_[idx];
                if (subgraph.vertex_set.find(edge.v_end) != subgraph.vertex_set.end()) 
                    continue;
                total_d_extending_edge++;
                open_set_neighbor.push_back(idx);
            }
        }
        prob_inverse *= total_d_extending_edge;

        if (total_d_extending_edge == 0) 
            return 0;

        uniform_int_distribution <uint32_t> next_edge_distribution(0, total_d_extending_edge - 1);
        last_edge_idx = open_set_neighbor[next_edge_distribution(rand_generator)];
        struct Edge extending_edge = graph.edge_list_[last_edge_idx];
        subgraph.insert_edge(extending_edge);

        // update open vertex set for chain
        open_vertex_set.erase(extending_edge.v_start);
        open_vertex_set.emplace(extending_edge.v_end, graph.v_table_.at(extending_edge.v_end));
    }

    // No test ConditionalComplete for chain

    return prob_inverse;
}

double ASAP_estimate_chain(Graph &graph, uint32_t chain_length, uint64_t sampling_times, default_random_engine & rand_generator)
{
    double total_prob_inverse = 0;
    for (uint64_t i = 0; i < sampling_times; i++) 
    {
        total_prob_inverse += ASAP_chain_neighbor_sampler(graph, chain_length, rand_generator);
    }

    return total_prob_inverse;
}

double ASAP_4_clique_neighbor_sampler(Graph &graph, default_random_engine & rand_generator)
{   
    
    double prob_inverse = 0;
    Subgraph subgraph;
    map<uint32_t, Vertex> open_vertex_set;
    
    uint32_t last_edge_idx = graph.ASAP_edge_sampling(rand_generator);
    Edge e1 = graph.edge_list_[last_edge_idx];
    subgraph.insert_edge(e1);
    open_vertex_set.emplace(e1.v_start, graph.v_table_.at(e1.v_start));
    open_vertex_set.emplace(e1.v_end, graph.v_table_.at(e1.v_end));

    prob_inverse = graph.total_edges_;

    uint32_t total_d_extending_edge = 0;
    uint32_t node1, node2, node3, node4, node3_;
    node1 = node2 = node3 = node3_ = node4 = 0;

    for (int j = 0; j < 2; j++)
    {
        // step 1: sampling two edges
        // step 3: sample third edge

        total_d_extending_edge = 0;
        vector<uint32_t> open_set_neighbor;
        for (auto it = open_vertex_set.begin(); it != open_vertex_set.end(); it++)
        {
            // uint32_t begin_idx = max(last_edge_idx + 1, it->second.edge_index);
            uint32_t begin_idx = it->second.edge_index;
            uint32_t end_idx = it->second.edge_index + it->second.degree;
            for (uint32_t idx = begin_idx; idx < end_idx; idx++) 
            {
                struct Edge edge = graph.edge_list_[idx];
                if (subgraph.vertex_set.find(edge.v_end) != subgraph.vertex_set.end()) 
                    continue;
                total_d_extending_edge++;
                open_set_neighbor.push_back(idx);
            }
        }
        prob_inverse *= total_d_extending_edge;

        if (total_d_extending_edge == 0) 
            return 0;

        uniform_int_distribution <uint32_t> next_edge_distribution(0, total_d_extending_edge - 1);
        last_edge_idx = open_set_neighbor[next_edge_distribution(rand_generator)];
        struct Edge extending_edge = graph.edge_list_[last_edge_idx];
        subgraph.insert_edge(extending_edge);

        
        if (j == 0)
        {
            node2 = extending_edge.v_start;
            node3 = extending_edge.v_end;
            if (node2 == e1.v_start)
                node1 = e1.v_end;
            else 
                node1 = e1.v_start;
            // step 2: test triangle
            if (!graph.neighbor_test(node1, node3)) 
                return 0;
        } 
        else 
        {
            node3_ = extending_edge.v_start;
            node4 = extending_edge.v_end;
        }

        open_vertex_set.emplace(extending_edge.v_end, graph.v_table_.at(extending_edge.v_end));
        
    }

    // step 4: test closure for 4-clique
    if (node3_ != node1)
        if (!graph.neighbor_test(node1, node4)) return 0;
    if (node3_ != node2)
        if (!graph.neighbor_test(node2, node4)) return 0;
    if (node3_ != node3) 
        if (!graph.neighbor_test(node3, node4)) return 0;    

    return prob_inverse;
}

double ASAP_estimate_4_clique(Graph &graph, uint64_t sampling_times, default_random_engine & rand_generator)
{
    double total_prob_inverse = 0;
    for (uint64_t i = 0; i < sampling_times; i++) 
    {
        total_prob_inverse += ASAP_4_clique_neighbor_sampler(graph, rand_generator);
    }

    return total_prob_inverse;
}

double ASAP_5_house_neighbor_sampler(Graph &graph, default_random_engine & rand_generator)
{  
    
    double prob_inverse = 0;
    Subgraph subgraph;
    map<uint32_t, Vertex> open_vertex_set;
    
    uint32_t last_edge_idx = graph.ASAP_edge_sampling(rand_generator);
    Edge e1 = graph.edge_list_[last_edge_idx];
    subgraph.insert_edge(e1);
    open_vertex_set.emplace(e1.v_start, graph.v_table_.at(e1.v_start));
    open_vertex_set.emplace(e1.v_end, graph.v_table_.at(e1.v_end));

    prob_inverse = graph.total_edges_;

    uint32_t total_d_extending_edge = 0;
    uint32_t node1, node2, node3, node4, node3_, node5, node4_;
    node1 = node2 = node3 = node3_ = node4 = node5 = node4_ = 0;
    vector<uint32_t> chain;

    for (int j = 0; j < 3; j++)
    {

        total_d_extending_edge = 0;
        vector<uint32_t> open_set_neighbor;
        for (auto it = open_vertex_set.begin(); it != open_vertex_set.end(); it++)
        {
            // uint32_t begin_idx = max(last_edge_idx + 1, it->second.edge_index);
            uint32_t begin_idx = it->second.edge_index;
            uint32_t end_idx = it->second.edge_index + it->second.degree;
            for (uint32_t idx = begin_idx; idx < end_idx; idx++) 
            {
                struct Edge edge = graph.edge_list_[idx];
                if (subgraph.vertex_set.find(edge.v_end) != subgraph.vertex_set.end()) 
                    continue;
                total_d_extending_edge++;
                open_set_neighbor.push_back(idx);
            }
        }
        prob_inverse *= total_d_extending_edge;

        if (total_d_extending_edge == 0) 
            return 0;

        uniform_int_distribution <uint32_t> next_edge_distribution(0, total_d_extending_edge - 1);
        last_edge_idx = open_set_neighbor[next_edge_distribution(rand_generator)];
        struct Edge extending_edge = graph.edge_list_[last_edge_idx];
        subgraph.insert_edge(extending_edge);

        
        if (j == 0)
        {
            node2 = extending_edge.v_start;
            node3 = extending_edge.v_end;
            if (node2 == e1.v_start)
                node1 = e1.v_end;
            else 
                node1 = e1.v_start;
            if (!graph.neighbor_test(node1, node3)) return 0;
            chain.push_back(node1);
            chain.push_back(node2);
            chain.push_back(node3);
        } 
        else if (j == 1)
        {
            node3_ = extending_edge.v_start;
            node4 = extending_edge.v_end;
            if (node3_ == chain[0])
                chain.insert(chain.begin(), node4);
            else
                chain.push_back(node4);
        }
        else if (j == 2)
        {
            node4_ = extending_edge.v_start;
            node5 = extending_edge.v_end;
            if (node4_ == chain[0])
                chain.insert(chain.begin(), node5);
            else    
                chain.push_back(node5);
        }

        open_vertex_set.erase(extending_edge.v_start);
        open_vertex_set.emplace(extending_edge.v_end, graph.v_table_.at(extending_edge.v_end));
    }

    // step 4: test closure for last 5_cycle edge
    if (!graph.neighbor_test(chain[0], chain[4])) return 0;

    return prob_inverse;
}

double ASAP_estimate_5_house(Graph &graph, uint64_t sampling_times, default_random_engine & rand_generator)
{
    double total_prob_inverse = 0;
    for (uint64_t i = 0; i < sampling_times; i++) 
    {
        total_prob_inverse += ASAP_5_house_neighbor_sampler(graph, rand_generator);
    }

    return total_prob_inverse;
}


double ASAP_triangle_triangle_neighbor_sampler(Graph &graph, default_random_engine & rand_generator)
{  
    
    double prob_inverse = 0;
    Subgraph subgraph;
    map<uint32_t, Vertex> open_vertex_set;
    
    uint32_t last_edge_idx = graph.ASAP_edge_sampling(rand_generator);
    Edge e1 = graph.edge_list_[last_edge_idx];
    subgraph.insert_edge(e1);
    open_vertex_set.emplace(e1.v_start, graph.v_table_.at(e1.v_start));
    open_vertex_set.emplace(e1.v_end, graph.v_table_.at(e1.v_end));

    prob_inverse = graph.total_edges_;

    uint32_t total_d_extending_edge = 0;
    uint32_t node1, node2, node3, node4, node3_, node5, node4_, node5_, node6;
    node1 = node2 = node3 = node3_ = node4 = node5 = node4_ = node5_ = node6 = 0;
    vector<uint32_t> chain;

    for (int j = 0; j < 4; j++)
    {
        // step 1: sampling two edges
        // step 3: sample third edge

        total_d_extending_edge = 0;
        vector<uint32_t> open_set_neighbor;
        for (auto it = open_vertex_set.begin(); it != open_vertex_set.end(); it++)
        {
            // uint32_t begin_idx = max(last_edge_idx + 1, it->second.edge_index);
            uint32_t begin_idx = it->second.edge_index;
            uint32_t end_idx = it->second.edge_index + it->second.degree;
            for (uint32_t idx = begin_idx; idx < end_idx; idx++) 
            {
                struct Edge edge = graph.edge_list_[idx];
                if (subgraph.vertex_set.find(edge.v_end) != subgraph.vertex_set.end()) 
                    continue;
                total_d_extending_edge++;
                open_set_neighbor.push_back(idx);
            }
        }
        prob_inverse *= total_d_extending_edge;

        if (total_d_extending_edge == 0) 
            return 0;

        uniform_int_distribution <uint32_t> next_edge_distribution(0, total_d_extending_edge - 1);
        last_edge_idx = open_set_neighbor[next_edge_distribution(rand_generator)];
        struct Edge extending_edge = graph.edge_list_[last_edge_idx];
        subgraph.insert_edge(extending_edge);

        
        if (j == 0)
        {
            node2 = extending_edge.v_start;
            node3 = extending_edge.v_end;
            if (node2 == e1.v_start)
                node1 = e1.v_end;
            else 
                node1 = e1.v_start;
            chain.push_back(node1);
            chain.push_back(node2);
            chain.push_back(node3);
        } 
        else if (j == 1)
        {
            node3_ = extending_edge.v_start;
            node4 = extending_edge.v_end;
            if (node3_ == chain[0])
                chain.insert(chain.begin(), node4);
            else
                chain.push_back(node4);
        }
        else if (j == 2)
        {
            node4_ = extending_edge.v_start;
            node5 = extending_edge.v_end;
            if (node4_ == chain[0])
                chain.insert(chain.begin(), node5);
            else    
                chain.push_back(node5);
        }
        else if (j == 3)
        {
            node5_ = extending_edge.v_start;
            node6 = extending_edge.v_end;
            if (node5_ == chain[0])
                chain.insert(chain.begin(), node6);
            else    
                chain.push_back(node6);
        }

        open_vertex_set.erase(extending_edge.v_start);
        open_vertex_set.emplace(extending_edge.v_end, graph.v_table_.at(extending_edge.v_end));
        
    }

    // step 4: test 2 remaining edges for 2 triangles
    if (!graph.neighbor_test(chain[0], chain[2])) return 0;
    if (!graph.neighbor_test(chain[5], chain[3])) return 0;
    return prob_inverse;
}

double ASAP_estimate_triangle_triangle(Graph &graph, uint64_t sampling_times, default_random_engine & rand_generator)
{
    double total_prob_inverse = 0;
    for (uint64_t i = 0; i < sampling_times; i++) 
    {
        total_prob_inverse += ASAP_triangle_triangle_neighbor_sampler(graph, rand_generator);
    }

    return total_prob_inverse;
}

#endif
